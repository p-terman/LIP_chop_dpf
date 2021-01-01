#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <exception>

#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>

#include "../Utilities/LCvtChannel.h"
#include "../Utilities/LCvtEvent.h"
#include "../Utilities/LCvtFile.h"
#include "../Utilities/RQFile_IO.hh"

#include "Timing.h"

using namespace std;
using boost::property_tree::ptree;

typedef vector<float> FVec;

////////////////////////////////////////////////////////////////////////////////
// Supporting Functions
////////////////////////////////////////////////////////////////////////////////
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
double FindPeakStd(LCvtChannel *ch, int pstart, int pend) {
  int pod_start, pod_end, npods = 0;
  double pre_std, pos_std, avg_pod_std = 0;
  for (size_t pp=0; pp<ch->pod_data_phe.size(); pp++) {
    pod_start = ch->pod_starts[pp];
    pod_end = pod_start+int(ch->pod_lengths[pp])-1;
    if (pod_start > pend or pod_end < pstart) continue;
    // find the standard deviation of a pod 
    pre_std = 0;
    for (size_t i=0; i<24; i++) 
      pre_std += ch->pod_data_phe[pp][i]*ch->pod_data_phe[pp][i]/24.0;
    pos_std = 0;
    for (int i=pod_end-pod_start; i>pod_end-pod_start-31; i--)
      pos_std += ch->pod_data_phe[pp][i]*ch->pod_data_phe[pp][i]/31.0;
    // add the bigger std 
    avg_pod_std += pre_std>pos_std ? pre_std : pos_std;    
    npods++;
  }
  if (npods) avg_pod_std = sqrt(avg_pod_std/npods);
  return avg_pod_std;
}
////////////////////////////////////////////////////////////////////////////////
// Main Function
////////////////////////////////////////////////////////////////////////////////
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
  ptree dp_settings;
  
  ifstream dp_settings_ifs(dp_xml_full_path.c_str());
  read_xml(dp_settings_ifs,dp_settings);
  dp_settings_ifs.close();
  
  // Fill variables from xml __________________________________________________
  ptree module_settings;
  FindModule(dp_settings, module_settings, "PulseTiming_HeightTiming");
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
  
  //Load data into RQFileIO object_____________________________________________  
  RQFileIO *rqio = new RQFileIO;
  rqio->ReadFile(rq_dir+"/"+rq_filename);
  
  // Add the RQs that this module is going to generate__________________________
  size_t pulse_dim=dp_settings.get<size_t>("data_processing_settings.global.max_num_pulses");
  string pulse_dim_str=dp_settings.get<string>("data_processing_settings.global.max_num_pulses");
  
  rqio->events.AddRQ("hft_t0_samples", "int32", pulse_dim_str, -999999);
  rqio->events.AddRQ("hft_t10l_samples", "int32", pulse_dim_str, -999999);
  rqio->events.AddRQ("hft_t50l_samples", "int32", pulse_dim_str, -999999);
  rqio->events.AddRQ("hft_t1_samples", "int32", pulse_dim_str, -999999);
  rqio->events.AddRQ("hft_t50r_samples", "int32", pulse_dim_str, -999999);
  rqio->events.AddRQ("hft_t10r_samples", "int32", pulse_dim_str, -999999);
  rqio->events.AddRQ("hft_t2_samples", "int32", pulse_dim_str, -999999);
  rqio->events.AddRQ("pulse_std_phe_per_sample", "float", pulse_dim_str);
  rqio->events.AddRQ("aft_t0_samples", "int32", pulse_dim_str, -999999);
  rqio->events.AddRQ("aft_t05_samples", "int32", pulse_dim_str, -999999);
  rqio->events.AddRQ("aft_tlx_samples", "int32", pulse_dim_str, -999999);
  rqio->events.AddRQ("aft_t25_samples", "int32", pulse_dim_str, -999999);
  rqio->events.AddRQ("aft_t1_samples", "int32", pulse_dim_str, -999999);
  rqio->events.AddRQ("aft_t75_samples", "int32", pulse_dim_str, -999999);
  rqio->events.AddRQ("aft_trx_samples", "int32", pulse_dim_str, -999999);
  rqio->events.AddRQ("aft_t95_samples", "int32", pulse_dim_str, -999999);
  rqio->events.AddRQ("aft_t2_samples", "int32", pulse_dim_str, -999999);
  
  // grap rqs 
  Rq* start = rqio->events.Get("pulse_start_samples");
  Rq* end = rqio->events.Get("pulse_end_samples");  
  Rq* t0 = rqio->events.Get("hft_t0_samples");
  Rq* t10l = rqio->events.Get("hft_t10l_samples");
  Rq* t50l = rqio->events.Get("hft_t50l_samples");
  Rq* t1 = rqio->events.Get("hft_t1_samples");
  Rq* t50r = rqio->events.Get("hft_t50r_samples");
  Rq* t10r = rqio->events.Get("hft_t10r_samples");
  Rq* t2 = rqio->events.Get("hft_t2_samples");
  Rq* pulse_std = rqio->events.Get("pulse_std_phe_per_sample");
  Rq* t01 = rqio->events.Get("aft_t0_samples");
  Rq* t05 = rqio->events.Get("aft_t05_samples");
  Rq* t10 = rqio->events.Get("aft_tlx_samples");
  Rq* t25 = rqio->events.Get("aft_t25_samples");
  Rq* t50 = rqio->events.Get("aft_t1_samples");
  Rq* t75 = rqio->events.Get("aft_t75_samples");
  Rq* t90 = rqio->events.Get("aft_trx_samples");
  Rq* t95 = rqio->events.Get("aft_t95_samples");
  Rq* t99 = rqio->events.Get("aft_t2_samples");
  
  
  //-----EVENT PROCESSING LOOP BEGINS-----
  double times[10], threshold, pstart, pend, tmp;
  LCvtChannel *ch;
  vector<float> pulse_waveform;
  for (size_t e=0; e<cvt->cvt_events.size(); e++) {
    rqio->GetEvent(e);
    
    for (size_t p=0; p<pulse_dim; p++) {
      if (end->GetInt(p)-start->GetInt(p) <= 0) break;

      // find the timing information 
      pstart = start->GetInt(p);
      pend = end->GetInt(p); 
      pulse_waveform = cvt->cvt_events[e]->channels[136]->GetPeak(pstart, pend);      
      
      // Calculate the pulse std 
      threshold = 0;
      for (size_t c=0; c<122; c++) {
        ch = cvt->cvt_events[e]->channels[c];
        tmp = FindPeakStd(ch, pstart, pend);
        threshold += tmp*tmp;
      }
      threshold = sqrt(threshold);
      HeightTiming(&pulse_waveform, threshold, times);
      
      // store the rq values 
      t0->SetDouble(times[0]+pstart, p);
      t10l->SetDouble(times[1]+pstart, p);
      t50l->SetDouble(times[2]+pstart, p);
      t1->SetDouble(times[3]+pstart, p);
      t50r->SetDouble(times[4]+pstart, p);
      t10r->SetDouble(times[5]+pstart, p);
      t2->SetDouble(times[6]+pstart, p);
      pulse_std->SetDouble(threshold, p);
      
      AreaTiming(&pulse_waveform, times);
      t01->SetDouble(times[0]+pstart, p);
      t05->SetDouble(times[1]+pstart, p);
      t10->SetDouble(times[2]+pstart, p);
      t25->SetDouble(times[3]+pstart, p);
      t50->SetDouble(times[4]+pstart, p);
      t75->SetDouble(times[5]+pstart, p);
      t90->SetDouble(times[6]+pstart, p);
      t95->SetDouble(times[7]+pstart, p);
      t99->SetDouble(times[8]+pstart, p);
      
    }
  }
  //-----EVENT PROCESSING LOOP ENDS-----
  
  // Write the output rqs
  rqio->WriteFile(rq_dir+"/"+rq_filename);

  //clean up___________________________________________________________________
  delete cvt;
  delete rqio;

  return 0;
}
