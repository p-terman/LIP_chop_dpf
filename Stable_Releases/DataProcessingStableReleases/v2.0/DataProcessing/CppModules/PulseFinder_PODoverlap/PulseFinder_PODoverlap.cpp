#include <iostream>
#include <fstream>
#include <string>

#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>

#include "../Utilities/LEvtEvent.h"
#include "../Utilities/LEvtFile.h"
#include "../Utilities/RQFile_IO.hh"

#include "PODoverlap.h"

using std::vector;
using std::cout;
using std::endl;
using std::string;
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
  
#if VERBOSITY>0
  cout<<evt_dir+"/"+evt_filename<<"\n"<<rq_dir+"/"+rq_filename<<"\n"<<dp_xml_full_path<<"\n"<<iq_xml_full_path<<endl;
#endif
  
  // Read input xml files into BOOST ptrees____________________________________
  ptree dp_settings, iqs;
  
  ifstream dp_settings_ifs(dp_xml_full_path.c_str());
  read_xml(dp_settings_ifs,dp_settings);
  dp_settings_ifs.close();
  
  ifstream iq_settings_ifs(iq_xml_full_path.c_str());
  read_xml(iq_settings_ifs,iqs);
  iq_settings_ifs.close();
  
  // Fill variables from xml __________________________________________________
  vector<double> mVns_per_phe;
  load_mVns_per_phe(iqs,mVns_per_phe);
  
  ptree module_settings;
  FindModule(dp_settings, module_settings, "PulseFinder_PODoverlap");
  cout<<module_settings.get<string>("relative_path")<<endl;  
  
  //Load data into LEvtFile object_____________________________________________
  LEvtFile* ef = new LEvtFile(evt_dir,evt_filename);
  ef->CalibrateChannels(mVns_per_phe);
  
  //Load data into RQFileIO object_____________________________________________  
  RQFileIO *rqio = new RQFileIO;
  rqio->ReadFile(rq_dir+"/"+rq_filename);
  
  // Add the RQs that this module is going to generate__________________________
  size_t pulse_dim = dp_settings.get<int>("data_processing_settings.global.max_num_pulses");
  string pulse_dim_str=dp_settings.get<string>("data_processing_settings.global.max_num_pulses");
  if (!rqio->events.Get("pulse_start_samples"))
    rqio->events.AddRQ("pulse_start_samples", "int32", pulse_dim_str);
  if (!rqio->events.Get("pulse_end_samples"))
    rqio->events.AddRQ("pulse_end_samples", "int32", pulse_dim_str);
  if (!rqio->events.Get("pulse_start_pod"))
    rqio->events.AddRQ("pulse_start_pod", "uint32", pulse_dim_str+",122");
  if (!rqio->events.Get("pulse_end_pod"))
    rqio->events.AddRQ("pulse_end_pod", "uint32", pulse_dim_str+",122");  
  if (!rqio->events.Get("pulse_length_samples"))
    rqio->events.AddRQ("pulse_length_samples", "uint32", pulse_dim_str);  
  // grap rqs 
  Rq* pulse_start_samples = rqio->events.Get("pulse_start_samples");
  Rq* pulse_end_samples = rqio->events.Get("pulse_end_samples");
  Rq* pulse_start_pod = rqio->events.Get("pulse_start_pod");
  Rq* pulse_end_pod = rqio->events.Get("pulse_end_pod");
  Rq* pulse_length_samples = rqio->events.Get("pulse_length_samples");
  
  //-----EVENT PROCESSING LOOP BEGINS-----
  vector<int> pstart, pend;
  vector< vector<int> > start_pod, end_pod;
  for (size_t e=0; e<ef->evt_events.size(); e++) {
    PODoverlap(ef->evt_events[e], pulse_dim, pstart, pend);
    WhichPODs(ef->evt_events[e], &pstart, &pend, start_pod, end_pod);
    rqio->GetEvent(e);
    
    for (size_t p=0; p<pulse_dim; p++) {
      if (p<pstart.size()) {
        pulse_start_samples->SetInt(pstart[p], p);
        pulse_end_samples->SetInt(pend[p], p);
        pulse_length_samples->SetInt(pend[p]-pstart[p]+1, p);
        for (size_t c=0; c<122; c++) {
          pulse_start_pod->SetInt(start_pod[p][c], p, c);
          pulse_end_pod->SetInt(end_pod[p][c], p, c);
        }
      }else {
        pulse_start_samples->SetInt(0, p);
        pulse_end_samples->SetInt(0, p);
        pulse_length_samples->SetInt(0, p);
        for (size_t c=0; c<122; c++) {
          pulse_start_pod->SetInt(0, p, c);
          pulse_end_pod->SetInt(0, p, c);
        }
      }
    }
  }
  //-----EVENT PROCESSING LOOP ENDS-----
  
  // Write the output rqs
  rqio->WriteFile(rq_dir+"/"+rq_filename);

  //clean up___________________________________________________________________
  delete ef;
  delete rqio;

  return 0;
}
