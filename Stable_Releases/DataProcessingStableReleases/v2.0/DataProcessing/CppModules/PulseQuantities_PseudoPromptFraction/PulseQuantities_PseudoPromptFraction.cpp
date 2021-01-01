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

#include "PseudoPromptFraction.h"

using namespace std;
using boost::property_tree::ptree;

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
  double pod_std, avg_pod_std = 0;
  for (size_t pp=0; pp<ch->pod_data_phe.size(); pp++) {
    pod_start = ch->pod_starts[pp];
    pod_end = pod_start+int(ch->pod_lengths[pp])-1;
    if (pod_start > pend or pod_end < pstart) continue;
    // find the standard deviation of a pod 
    pod_std = 0;
    for (size_t i=0; i<24; i++) 
      pod_std += ch->pod_data_phe[pp][i]*ch->pod_data_phe[pp][i];
    avg_pod_std += pod_std/24;
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
  FindModule(dp_settings, module_settings, "PulseQuantities_PseudoPromptFraction");
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
  
  // Get some module parameters
  size_t nsamples;
  try { nsamples = module_settings.get<size_t>("parameters.max_n_samples");}
  catch (exception const &ex) { nsamples = 9;}
  
  //Load data into RQFileIO object_____________________________________________  
  RQFileIO *rqio = new RQFileIO;
  rqio->ReadFile(rq_dir+"/"+rq_filename);
  
  // Add the RQs that this module is going to generate__________________________
  size_t pulse_dim=dp_settings.get<size_t>("data_processing_settings.global.max_num_pulses");
  string pulse_dim_str=dp_settings.get<string>("data_processing_settings.global.max_num_pulses");
  rqio->events.AddRQ("pseudo_prompt_fraction", "float", pulse_dim_str);
  rqio->events.AddRQ("pulse_std_phe_per_sample_ppf", "float", pulse_dim_str);
  
  // grap rqs 
  Rq* pstart = rqio->events.Get("pulse_start_samples");
  Rq* pend = rqio->events.Get("pulse_end_samples");  
  Rq* t0 = rqio->events.Get("hft_t0_samples");
  Rq* t10l = rqio->events.Get("hft_t10l_samples");
  Rq* t2 = rqio->events.Get("hft_t2_samples");
  Rq* spf = rqio->events.Get("pseudo_prompt_fraction");
  Rq* pulse_std = rqio->events.Get("pulse_std_phe_per_sample_ppf");  
  
  //rqio->events.PrintRQs();
  
  //-----EVENT PROCESSING LOOP BEGINS-----
  double results[3];
  int start, end, offset;
  LCvtChannel *ch;
  vector<float> waveform;
  for (size_t e=0; e<cvt->cvt_events.size(); e++) {
    rqio->GetEvent(e);
    //cout << "Entry #" << e+1 << endl;
    // pulse level rqs 
    for (size_t p=0; p<pulse_dim; p++) {
      /*
      cout << "Pulse " << p+1 
           << " pulse_start=" << pstart->GetDouble(p) 
           << " pulse_end=" << pend->GetDouble(p) 
           << " hft_t0=" << t0->GetDouble(p) 
           << " hft_t10l=" << t10l->GetDouble(p) 
           << " hft_t2=" << t2->GetDouble(p) << endl;
      */
      
      if (pend->GetDouble(p)-pstart->GetDouble(p) <= 0) break;
      if (!isfinite(pstart->GetDouble(p))) continue;
      
      start = t0->GetDouble(p) - pstart->GetDouble(p);
      end = t2->GetDouble(p) - pstart->GetDouble(p);
      offset = t10l->GetDouble(p) - pstart->GetDouble(p);
      ch = cvt->cvt_events[e]->channels[136];
      waveform = ch->GetPeak(pstart->GetDouble(p), pend->GetDouble(p));
      /*
      cout << "start=" << start << " offset=" << offset 
           << " end=" << end << " waveform size=" << waveform.size() << endl;
      */
      // Calculate the pulse std 
      results[0] = 0;
      for (size_t c=0; c<122; c++) {
        ch = cvt->cvt_events[e]->channels[c];
        results[2] = FindPeakStd(ch, pstart->GetDouble(p), pend->GetDouble(p));
        results[0] += results[2]*results[2];
      }
      results[0] = sqrt(results[0]);
      
      // Calculate the pseudo_prompt_fraction
      results[1] = PseudoPromptFraction(&waveform, start, end, offset, nsamples, results[0]);
      
      // Store the new rqs 
      pulse_std->SetDouble(results[0], p);
      spf->SetDouble(results[1], p);
      
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
