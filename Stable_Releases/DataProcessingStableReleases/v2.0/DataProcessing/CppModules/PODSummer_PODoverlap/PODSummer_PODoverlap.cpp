#include <iostream>
#include <fstream>
#include <string>

#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>

#include "../Utilities/LEvtEvent.h"
#include "../Utilities/LEvtFile.h"
#include "../Utilities/LCvtEvent.h"
#include "../Utilities/LCvtFile.h"

#include "PODoverlap.h"

using namespace std;
using boost::property_tree::ptree;

//______________________________________________________________________________
void FindModule(ptree &full_tree, ptree& module, string module_name) {  
  string dp = "data_processing_settings";
  string tag = "module_name";
  BOOST_FOREACH(ptree::value_type &v, full_tree.get_child(dp)) {
    if (v.second.count(tag)==1) {
      if (v.second.get<string>(tag)==module_name){
        module = v.second;
      }
    }
  }  
}
//______________________________________________________________________________
int main(int argc, char* argv[]) {
  
  //parse inputs_______________________________________________________________
  if( argc !=8  ) {
      cerr << "Please specify:\n-evt_filename, evt_dir, rq_filename, rq_dir, module no., dp_xml_full_path, iq_xml_full_path" << endl;
      return -1;
  }
  
  string evt_dir = string(argv[2]);
  string evt_filename = string(argv[1]);
  string rq_dir = string(argv[4]);
  string rq_filename =string(argv[3]);
  string dp_xml_full_path = string(argv[6]);
  string iq_xml_full_path = string(argv[7]);
  
  // Read input xml files into BOOST ptrees____________________________________
  ptree dp_settings, iqs;
  
  ifstream dp_settings_ifs(dp_xml_full_path.c_str());
  read_xml(dp_settings_ifs,dp_settings);
  dp_settings_ifs.close();
  
  ifstream iq_settings_ifs(iq_xml_full_path.c_str());
  read_xml(iq_settings_ifs,iqs);
  iq_settings_ifs.close();
  
  // Fill variables from xml __________________________________________________  
  ptree module_settings;
  FindModule(dp_settings, module_settings, "PODSummer_PODoverlap");
  cout<<module_settings.get<string>("relative_path")<<endl;
  
  
  //Load data into LCvtFile object_____________________________________________
  string cvt_filename(evt_filename);
  cvt_filename.replace(cvt_filename.rfind(".evt"), 4, ".cvt");
  {
    // This code quickly checks if a cvt file exist 
    string name = evt_dir + "/" + cvt_filename;
    ifstream tmp(name.c_str());
    if (!tmp.good()) {
      cerr << "Error: Cannot find cvt file.\n";
      tmp.close();
      return 1;
    }
    tmp.close();
  }  
  LCvtFile* cvt = new LCvtFile(evt_dir, cvt_filename);
  
  //-----EVENT PROCESSING LOOP BEGINS-----
  vector<int> pstart, pend;
  for (size_t e=0; e<cvt->cvt_events.size(); e++) {
    PODoverlapXe(cvt->cvt_events[e], pstart, pend);
    cvt->cvt_events[e]->AddXeSumPodCh(pstart, pend);
    PODoverlapWa(cvt->cvt_events[e], pstart, pend);
    cvt->cvt_events[e]->AddWaSumPodCh(pstart, pend);
  }
  //-----EVENT PROCESSING LOOP ENDS-----
  
  // Write the output 
  cvt->WriteCvtFile(evt_dir, cvt_filename);
  
  //clean up___________________________________________________________________
  delete cvt;

  return 0;
}
