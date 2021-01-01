#include <iostream>
#include <fstream>
#include <string>

#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>

#include "../Utilities/LEvtEvent.h"
#include "../Utilities/LEvtFile.h"
#include "../Utilities/LCvtEvent.h"
#include "../Utilities/LCvtFile.h"

using std::vector;
using std::cout;
using std::endl;
using std::string;
using std::cerr;
using std::ifstream;
using boost::property_tree::ptree;

//______________________________________________________________________________
void load_mVns_per_phe(ptree& pt, vector<double> &mVns_per_phe){
  mVns_per_phe.clear();
  mVns_per_phe.resize(122);
  for(ptree::iterator iter=pt.begin();iter!=pt.end();iter++){
    if (iter->second.count("global")==0) continue;
    if (iter->second.get_child("global").count("iq_type")==0) continue;
    if (iter->second.get_child("global").get<string>("iq_type").compare("pmt_gains")) continue;
    if (iter->second.count("fit")==0) continue;
    
    BOOST_FOREACH(ptree::value_type &v,iter->second.get_child("fit")){
      if(v.first=="channel"){
        int ind = v.second.get<int>("PMT_number")-1;
        mVns_per_phe[ind]=v.second.get<double>("mVns_per_phe");
      }
    }
  }
}
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
void Convert(vector<float>* vec, float baseline, float scale) {
  for (size_t i=0; i<vec->size(); i++) {
    vec->at(i) = (baseline - vec->at(i))*scale;
  }
}
//______________________________________________________________________________
int main(int argc, char* argv[]) {
  
  //parse inputs
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
  
  // Read input xml files into BOOST ptrees
  ptree dp_settings, iqs;
  
  ifstream dp_settings_ifs(dp_xml_full_path.c_str());
  read_xml(dp_settings_ifs,dp_settings);
  dp_settings_ifs.close();
  
  ifstream iq_settings_ifs(iq_xml_full_path.c_str());
  read_xml(iq_settings_ifs,iqs);
  iq_settings_ifs.close();
  
  // Fill variables from xml
  ptree module_settings;
  FindModule(dp_settings, module_settings, "PulseCalibration_MakeCvtFile");
  cout<<module_settings.get<string>("relative_path")<<endl;  
  
  // Fill variables from xml __________________________________________________
  vector<double> mVns_per_phe;
  load_mVns_per_phe(iqs,mVns_per_phe);
  
  // Load evt file data into LCvtFile object
  LEvtFile* ef = new LEvtFile(evt_dir,evt_filename);
  LCvtFile* cvt = new LCvtFile(ef);
  delete ef;
  
  // Convert ADC to phe 
  cvt->Convert(mVns_per_phe);  
  
  //Write cvt file 
  string cvt_filename(evt_filename);
  cvt_filename.replace(cvt_filename.rfind(".evt"), 4, ".cvt");
  cvt->WriteCvtFile(evt_dir, cvt_filename);
	
  //clean up  
  delete cvt;

  return 0;
}
