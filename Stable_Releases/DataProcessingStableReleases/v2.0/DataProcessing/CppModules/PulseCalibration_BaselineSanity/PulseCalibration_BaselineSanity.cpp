#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>

#include "../Utilities/LEvtEvent.h"
#include "../Utilities/LEvtFile.h"
#include "../Utilities/LCvtEvent.h"
#include "../Utilities/LCvtFile.h"
#include "../Utilities/RQFile_IO.hh"

#include "BaselineSanity.h"

using namespace std;
using boost::property_tree::ptree;

////////////////////////////////////////////////////////////////////////////////
// Supporting Functions
////////////////////////////////////////////////////////////////////////////////
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
float Max(vector<float>* vec) {
  float max = 0;
  if (vec->size()) vec->at(0);
  for (size_t i=1; i<vec->size(); i++) {
    if (max < (*vec)[i]) max = (*vec)[i];
  }
  return max;
}
////////////////////////////////////////////////////////////////////////////////
// Main Function 
////////////////////////////////////////////////////////////////////////////////
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
  FindModule(dp_settings, module_settings, "PulseCalibration_BaselineSanity");
  cout<<module_settings.get<string>("relative_path")<<endl;  
  
  // Fill variables from xml __________________________________________________
  vector<double> mVns_per_phe;
  load_mVns_per_phe(iqs, mVns_per_phe);
  
  // Get some module parameters
  float daq_saturation_mV, pmt_saturation_mV, max_shift_mV;
  size_t min_width;
  try { daq_saturation_mV = module_settings.get<float>("parameters.daq_saturation_threshold_mV"); }
  catch (exception const& ex) { daq_saturation_mV = 1800; }
  try { pmt_saturation_mV = module_settings.get<float>("parameters.pmt_saturation_2pct_threshold_mV");}
  catch (exception const& ex) { pmt_saturation_mV = 1600; }
  try { max_shift_mV = module_settings.get<size_t>("parameters.max_shift_mV"); }
  catch (exception const& ex) { max_shift_mV = 20; }
  try { min_width = module_settings.get<size_t>("parameters.min_window_size_samples"); }
  catch (exception const& ex) { min_width = 20; }
  
  // Load evt file, cvt file and rq file 
  LEvtFile* ef = new LEvtFile(evt_dir,evt_filename);
  
  RQFileIO *rqio = new RQFileIO;
  rqio->ReadFile(rq_dir+"/"+rq_filename);
  
  LCvtFile* cvt = new LCvtFile(ef);
  delete ef;
  cvt->Convert(mVns_per_phe);
  
  // Add new RQ 
  rqio->events.AddRQ("max_baseline_shift_mV", "float", "122");
  rqio->events.AddRQ("daq_saturation_flag", "bool", "122");
  rqio->events.AddRQ("pmt_2pct_saturation_flag", "bool", "122");
  
  // Grab rqs 
  Rq* shift = rqio->events.Get("max_baseline_shift_mV");
  Rq* daq_sat = rqio->events.Get("daq_saturation_flag");
  Rq* pmt_sat = rqio->events.Get("pmt_2pct_saturation_flag");
  
  //-----EVENT PROCESSING LOOP BEGINS-----
  LCvtChannel *ch;
  double local_max_diff, max_diff, max_shift_phe, local_max_height, max_height; 
  for (size_t e=0; e<cvt->cvt_events.size(); e++) {
    rqio->GetEvent(e);
    
    // loop through the channels
    for (unsigned int cc=0; cc<122; cc++) {
      ch = cvt->cvt_events[e]->channels[cc];
      max_diff = 0;
      max_shift_phe = max_shift_mV*10.0/mVns_per_phe[cc];
      max_height = 0;
      
      // loop through the pods 
      for (unsigned int pp=0; pp<ch->number_of_pods; pp++) {
        
        local_max_height = Max(&ch->pod_data_phe[pp]);
        if (max_height < local_max_height) max_height = local_max_height;
        
        local_max_diff = BaselineSanity(&ch->pod_data_phe[pp], max_shift_phe, min_width);
        if (max_diff < local_max_diff) max_diff = local_max_diff;
        
      } // end of pod loop 
      
      // Store the rqs  
      max_height *= mVns_per_phe[cc]/10.0;
      if (max_height >= daq_saturation_mV) daq_sat->SetInt(1, cc);
      if (max_height >= pmt_saturation_mV) pmt_sat->SetInt(1, cc);
      
      max_diff *= mVns_per_phe[cc]/10.0;
      shift->SetDouble(max_diff, cc);
    }  
  }
  //-----EVENT PROCESSING LOOP ENDS----- 
  
  // Write cvt file and rq file 
  string cvt_filename(evt_filename);
  cvt_filename.replace(cvt_filename.rfind(".evt"), 4, ".cvt");
  cvt->WriteCvtFile(evt_dir, cvt_filename);
	rqio->WriteFile(rq_dir+"/"+rq_filename);
	
  //clean up  
  delete cvt;  
  delete rqio;
  
  return 0;
}
