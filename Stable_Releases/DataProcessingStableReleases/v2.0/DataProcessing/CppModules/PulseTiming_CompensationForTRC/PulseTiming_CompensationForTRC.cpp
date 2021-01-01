#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <exception>
#include <numeric>

#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>

#include "../Utilities/LCvtChannel.h"
#include "../Utilities/LCvtEvent.h"
#include "../Utilities/LCvtFile.h"
#include "../Utilities/RQFile_IO.hh"

#include "TRCAlgorithms.h"

using namespace std;
using boost::property_tree::ptree;

typedef vector<float> FVec;

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
  ptree dp_settings;
  
  ifstream dp_settings_ifs(dp_xml_full_path.c_str());
  read_xml(dp_settings_ifs,dp_settings);
  dp_settings_ifs.close();
  
  // Fill variables from xml __________________________________________________
  ptree module_settings;
  FindModule(dp_settings, module_settings, "PulseTiming_CompensationForTRC");
  cout<<module_settings.get<string>("relative_path")<<endl;  
  
  double threshold, ratio;
  int width;
  bool change;  
  try {	
  	ratio = module_settings.get<double>("parameters.txFraction");
  } catch (exception const& ex) { ratio = 0.1; }
  try {	
  	width = module_settings.get<int>("parameters.skinnyBoxSamples");
  } catch (exception const& ex) { width = 10; }
  try {	
  	threshold = module_settings.get<double>("parameters.noiseThre");
  } catch (exception const& ex) { threshold = 0.15; }
  try {	
  	change = module_settings.get<bool>("parameters.changePulseStartEnd");
  } catch (exception const& ex) { change = true; }
  
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
  
  if (!rqio->events.Get("skinny_pulse_start_samples"))
    rqio->events.AddRQ("skinny_pulse_start_samples", "double", pulse_dim_str);
  if (!rqio->events.Get("skinny_pulse_end_samples"))
    rqio->events.AddRQ("skinny_pulse_end_samples", "double", pulse_dim_str);
  if (!rqio->events.Get("aft_tlx_samples"))
    rqio->events.AddRQ("aft_tlx_samples", "double", pulse_dim_str);
  if (!rqio->events.Get("aft_trx_samples"))
    rqio->events.AddRQ("aft_trx_samples", "double", pulse_dim_str);
  if (!rqio->events.Get("n_samples_in_evt"))
    rqio->events.AddRQ("n_samples_in_evt", "double", "1");
  if (!rqio->events.Get("full_evt_area_phe"))
    rqio->events.AddRQ("full_evt_area_phe", "double", "1");
  
  // grap rqs 
  Rq* start = rqio->events.Get("pulse_start_samples");
  Rq* end = rqio->events.Get("pulse_end_samples");
  Rq* t0 = rqio->events.Get("hft_t0_samples");
  Rq* t2 = rqio->events.Get("hft_t2_samples");
  Rq* tlx = rqio->events.Get("aft_tlx_samples");
  Rq* trx = rqio->events.Get("aft_trx_samples");
  Rq* skinny_start = rqio->events.Get("skinny_pulse_start_samples");
  Rq* skinny_end = rqio->events.Get("skinny_pulse_end_samples");
  Rq* nsamples = rqio->events.Get("n_samples_in_evt");
  Rq* eventArea = rqio->events.Get("full_evt_area_phe");
  
  
  //-----EVENT PROCESSING LOOP BEGINS-----
  double pstart, pend, tmp[2];
  LCvtChannel *ch;
  vector<float> pulse_waveform;
  for (size_t e=0; e<cvt->cvt_events.size(); e++) {
    rqio->GetEvent(e);     
    ch = cvt->cvt_events[e]->channels[136];
    
    for (size_t p=0; p<pulse_dim; p++) {
      if (end->GetInt(p)-start->GetInt(p) <= 0) break;
      // find the timing information 
      pstart = start->GetDouble(p);//t0->GetInt(p);
      pend = end->GetDouble(p);//t2->GetInt(p);
      
      pulse_waveform = ch->GetPeak(pstart, pend); 
      
      // change all the pulse start and pulse end rqs to the hft t0 and t2 
      // which is a more accurate algorithm. This is important to get consistent
      // results between FastMinimumSet and MinimumSet modules 
      /*if (change) {
      	start->SetDouble(pstart, p);
      	end->SetDouble(pend, p);
      }*/    
      
      // Calculate the skinny box rqs 
      tmp[0] = tmp[1] = 0;
      SkinnyPulseTiming(&pulse_waveform, threshold, width, tmp[0], tmp[1]);
      skinny_start->SetDouble(tmp[0]+pstart, p);
      skinny_end->SetDouble(tmp[1]+pstart, p);
      
      // Calculate the tlx and trx  rqs 
      tmp[0] = tmp[1] = 0;
      AftTxTiming(&pulse_waveform, threshold, ratio, tmp[0], tmp[1]);
      if (tmp[0] == -1 || tmp[1] == -1) {
      	tlx->SetDouble(pstart, p);
      	trx->SetDouble(pend, p);
      }else {
      	tlx->SetDouble(tmp[0]+pstart, p);
      	trx->SetDouble(tmp[1]+pstart, p);
      }
      
    }
    
    // Calculate the event level rqs 
    tmp[0] = tmp[1] = 0;
    for (size_t sp=0; sp<ch->pod_starts.size(); sp++) {
    	tmp[0] += ch->pod_data_phe[sp].size();    	
    	tmp[1] = accumulate(ch->pod_data_phe[sp].begin(), ch->pod_data_phe[sp].end(), tmp[1]);
    }
    nsamples->SetDouble(tmp[0]);
    eventArea->SetDouble(tmp[1]);
    
  }
  //-----EVENT PROCESSING LOOP ENDS-----
  
  // Write the output rqs
  rqio->WriteFile(rq_dir+"/"+rq_filename);

  //clean up___________________________________________________________________
  delete cvt;
  delete rqio;

  return 0;
}
