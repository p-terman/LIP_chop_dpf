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

#include "Fits.h"
#include "Quantities.h"

using namespace std;
using boost::property_tree::ptree;
double small_number = 1e-10;
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
//______________________________________________________________________________
void PeakBaselineStats(LCvtChannel *ch, int pstart, int pend, double *results) {
  // Results are 
  // [0] average pod predetect mean 
  // [1] average pod predetect std 
  // [2] average pod postdetect mean 
  // [3] average pod postdetect std 
  
  int pod_start, pod_end, npods = 0;
  double tmp[4];
  results[0] = results[1] = results[2] = results[3] = results[4] = 0;
  for (size_t pp=0; pp<ch->pod_data_phe.size(); pp++) {
    pod_start = ch->pod_starts[pp];
    pod_end = pod_start+int(ch->pod_lengths[pp])-1;
    if (pod_start > pend or pod_end < pstart) continue;
    
    PODBaselineStats(&ch->pod_data_phe[pp], tmp);
    results[0] += tmp[0];
    results[1] += tmp[1]*tmp[1];
    results[2] += tmp[2];
    results[3] += tmp[3]*tmp[3];
    results[4] += ch->pod_baselines[pp];
    npods++;
  }
  if (npods) {
    results[0] /= npods;
    results[1] = sqrt(results[1]/npods);
    results[2] /= npods;
    results[3] = sqrt(results[3]/npods);
    results[4] /= npods;
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
  vector<double> mVns_per_phe;
  load_mVns_per_phe(iqs,mVns_per_phe);
   
  ptree module_settings;
  FindModule(dp_settings, module_settings, "PulseQuantities_FastMinimumSet");
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
  //cout << "About to load cvt..." << endl;
  LCvtFile* cvt = new LCvtFile(evt_dir, cvt_filename);
  //cout << "done." << endl;
  
  // Get some cvt file settings 
  float preamp = cvt->settings_pt.get<float>("daq_settings.global.preamp");
  float postamp = cvt->settings_pt.get<float>("daq_settings.global.postamp");
  
  // Get some module parameters
  size_t prebins, window, nsamples, S1_width, S2_width;
  float daq_saturation_mV, pmt_saturation_mV;
  try { prebins = module_settings.get<size_t>("parameters.prompt_fraction.preBins_samples");}
  catch (exception const &ex) { prebins = 2;}
  try { window = module_settings.get<size_t>("parameters.prompt_fraction.windowBins_samples");}
  catch (exception const &ex) { window = 10;}
  try { nsamples = module_settings.get<size_t>
                           ("parameters.pseudo_prompt_fraction.max_n_samples");}
  catch (exception const &ex) { nsamples = 9;}
  try { daq_saturation_mV = module_settings.get<float>("parameters.saturation_flags.daq_saturation_threshold_mV");}
  catch (exception const &ex) { daq_saturation_mV = 1800;}
  try { pmt_saturation_mV = module_settings.get<float>("parameters.saturation_flags.pmt_saturation_2pct_threshold_mV");}
  catch (exception const &ex) { pmt_saturation_mV = 1600;}
  try { S1_width = module_settings.get<size_t>("parameters.s2filter.s1window_samples");}
  catch (exception const &ex) { S1_width = 50;}
  try { S2_width = module_settings.get<size_t>("parameters.s2filter.s2window_samples");}
  catch (exception const &ex) { S2_width = 200;}
  
  //Load data into RQFileIO object_____________________________________________  
  //cout << "About to load RQ file." << endl;
  RQFileIO *rqio = new RQFileIO;
  rqio->ReadFile(rq_dir+"/"+rq_filename);
  //cout << "done." << endl;
  
  // Add the RQs that this module is going to generate__________________________
  size_t pulse_dim=dp_settings.get<size_t>("data_processing_settings.global.max_num_pulses");
  string pulse_dim_str=dp_settings.get<string>("data_processing_settings.global.max_num_pulses");

  rqio->events.AddRQ("pulse_area_phe", "float", pulse_dim_str);
  rqio->events.AddRQ("pulse_height_phe_per_sample", "float", pulse_dim_str);
  rqio->events.AddRQ("prompt_fraction", "float", pulse_dim_str);
  rqio->events.AddRQ("pseudo_prompt_fraction", "float", pulse_dim_str);
  rqio->events.AddRQ("pulse_std_phe_per_sample", "float", pulse_dim_str);  
  rqio->events.AddRQ("exp_fit_amplitude_phe_per_sample", "float", pulse_dim_str);
  rqio->events.AddRQ("exp_fit_tau_fall_samples", "float", pulse_dim_str);
  rqio->events.AddRQ("exp_fit_time_offset_samples", "float", pulse_dim_str);
  rqio->events.AddRQ("exp_fit_tau_rise_samples", "float", pulse_dim_str);
  rqio->events.AddRQ("exp_fit_chisq", "float", pulse_dim_str);
  rqio->events.AddRQ("exp_fit_dof", "float", pulse_dim_str);
  rqio->events.AddRQ("gaus_fit_amplitude_phe_per_sample", "float", pulse_dim_str);
  rqio->events.AddRQ("gaus_fit_mu_samples", "float", pulse_dim_str);
  rqio->events.AddRQ("gaus_fit_sigma_samples", "float", pulse_dim_str);
  rqio->events.AddRQ("gaus_fit_chisq", "float", pulse_dim_str);
  rqio->events.AddRQ("gaus_fit_dof", "float", pulse_dim_str);
  rqio->events.AddRQ("peak_area_phe", "float", pulse_dim_str+",122");
  rqio->events.AddRQ("peak_height_mV", "float", pulse_dim_str+",122");
  rqio->events.AddRQ("peak_height_phe_per_sample", "float", pulse_dim_str+",122");
  rqio->events.AddRQ("mean_first_last_pts_phe_per_sample", "float", pulse_dim_str+",122");
  rqio->events.AddRQ("pod_detect_pre_pod_mean_phe", "float", pulse_dim_str+",122");
  rqio->events.AddRQ("pod_detect_pos_pod_mean_phe", "float", pulse_dim_str+",122");
  rqio->events.AddRQ("pod_detect_pre_pod_std_phe", "float", pulse_dim_str+",122");
  rqio->events.AddRQ("pod_detect_pos_pod_std_phe", "float", pulse_dim_str+",122");
  rqio->events.AddRQ("baseline_daq_mV", "float", pulse_dim_str+",122");
  rqio->events.AddRQ("daq_saturation_flag", "bool", "122");
  rqio->events.AddRQ("pmt_2pct_saturation_flag", "bool", "122");
  rqio->events.AddRQ("top_bottom_ratio", "float", pulse_dim_str, -99);
  rqio->events.AddRQ("top_bottom_asymmetry", "float", pulse_dim_str, -99);
  rqio->events.AddRQ("pulse_area_positive_phe", "float", pulse_dim_str);
  rqio->events.AddRQ("pulse_area_negative_phe", "float", pulse_dim_str);
  rqio->events.AddRQ("pulse_top_array_area_phe", "float", pulse_dim_str);
  rqio->events.AddRQ("pulse_bot_array_area_phe", "float", pulse_dim_str);  
  rqio->events.AddRQ("s2filter_max_s1_area", "float", pulse_dim_str);
  rqio->events.AddRQ("s2filter_max_s2_area", "float", pulse_dim_str);
  rqio->events.AddRQ("s2filter_max_area_diff", "float", pulse_dim_str);
  rqio->events.AddRQ("max_peak_area_phe", "float", pulse_dim_str);
  
  // grap rqs 
  Rq* pstart = rqio->events.Get("pulse_start_samples");
  Rq* pend = rqio->events.Get("pulse_end_samples");  
  Rq* t0 = rqio->events.Get("hft_t0_samples");
  Rq* t10l = rqio->events.Get("hft_t10l_samples");
  Rq* t2 = rqio->events.Get("hft_t2_samples");
  Rq* pulse_std = rqio->events.Get("pulse_std_phe_per_sample");
  Rq* parea = rqio->events.Get("pulse_area_phe");
  Rq* pheight = rqio->events.Get("pulse_height_phe_per_sample");
  Rq* pf = rqio->events.Get("prompt_fraction");
  Rq* spf = rqio->events.Get("pseudo_prompt_fraction");
  Rq* exp_amp = rqio->events.Get("exp_fit_amplitude_phe_per_sample");
  Rq* exp_fall = rqio->events.Get("exp_fit_tau_fall_samples");
  Rq* exp_t0 = rqio->events.Get("exp_fit_time_offset_samples");
  Rq* exp_rise = rqio->events.Get("exp_fit_tau_rise_samples");
  Rq* exp_chisq = rqio->events.Get("exp_fit_chisq");
  Rq* exp_dof = rqio->events.Get("exp_fit_dof");
  Rq* gaus_amp = rqio->events.Get("gaus_fit_amplitude_phe_per_sample");
  Rq* gaus_mu = rqio->events.Get("gaus_fit_mu_samples");
  Rq* gaus_sigma = rqio->events.Get("gaus_fit_sigma_samples");
  Rq* gaus_chisq = rqio->events.Get("gaus_fit_chisq");
  Rq* gaus_dof = rqio->events.Get("gaus_fit_dof");
  Rq* pparea = rqio->events.Get("peak_area_phe");
  Rq* ppheight_mV = rqio->events.Get("peak_height_mV");
  Rq* ppheight = rqio->events.Get("peak_height_phe_per_sample");
  Rq* f_l_points = rqio->events.Get("mean_first_last_pts_phe_per_sample");
  Rq* pre_mean = rqio->events.Get("pod_detect_pre_pod_mean_phe");
  Rq* pos_mean = rqio->events.Get("pod_detect_pos_pod_mean_phe");
  Rq* pre_std = rqio->events.Get("pod_detect_pre_pod_std_phe");
  Rq* pos_std = rqio->events.Get("pod_detect_pos_pod_std_phe");
  Rq* baseline = rqio->events.Get("baseline_daq_mV");
  Rq* daq_sat = rqio->events.Get("daq_saturation_flag");
  Rq* pmt_sat = rqio->events.Get("pmt_2pct_saturation_flag");
  Rq* tb_ratio = rqio->events.Get("top_bottom_ratio");
  Rq* tb_asym = rqio->events.Get("top_bottom_asymmetry");
  Rq* pulse_pos_area = rqio->events.Get("pulse_area_positive_phe");
  Rq* pulse_neg_area = rqio->events.Get("pulse_area_negative_phe");
  Rq* top_area = rqio->events.Get("pulse_top_array_area_phe");
  Rq* bot_area = rqio->events.Get("pulse_bot_array_area_phe");
  Rq* s2filter_s1 = rqio->events.Get("s2filter_max_s1_area");
  Rq* s2filter_s2 = rqio->events.Get("s2filter_max_s2_area");
  Rq* s2filter_diff = rqio->events.Get("s2filter_max_area_diff");
  Rq* max_peak_area = rqio->events.Get("max_peak_area_phe");
  
  
  //-----EVENT PROCESSING LOOP BEGINS-----
  double results[12], gain, threshold;
  size_t start, end, offset;
  LCvtChannel *ch;
  vector<float> waveform;
  //cout << "About to loop" << endl;
  for (size_t e=0; e<cvt->cvt_events.size(); e++) {
    rqio->GetEvent(e);
    //cout << "Working on event " << e << " of " << cvt->cvt_events.size() << endl;
    
    // pulse level rqs 
    for (size_t p=0; p<pulse_dim; p++) {
      if (pend->GetInt(p)-pstart->GetInt(p) <= 0) break;
      
      start = t0->GetDouble(p) - pstart->GetInt(p);
      end = t2->GetDouble(p) - pstart->GetInt(p);
      offset = int(t10l->GetDouble(p)) - int(pstart->GetInt(p));
      ch = cvt->cvt_events[e]->channels[136];
      waveform = ch->GetPeak(pstart->GetInt(p), pend->GetInt(p));
      
      // calculate areas and height       
      //PulseBasicQuantites(&waveform, start, end, results);
      PulseBasicQuantites(&waveform, 0, waveform.size()-1, results);
      parea->SetDouble(results[0], p);
      pheight->SetDouble(results[1], p);
      pulse_pos_area->SetDouble(results[2], p);
      pulse_neg_area->SetDouble(results[3], p);
      
      // if pulse_std isn't defined, calculate it. 
      if (!pulse_std->GetDouble(p)) {
        threshold = 0;
        for (size_t c=0; c<122; c++) {
          ch = cvt->cvt_events[e]->channels[c];
          results[5] = FindPeakStd(ch, pstart->GetInt(p), pend->GetInt(p));
          threshold += results[5]*results[5];
        }
        threshold = sqrt(threshold);
        pulse_std->SetDouble(threshold, p);
      }else {
        threshold = pulse_std->GetDouble(p);
      }
      
      // calculate prompt fractions
      results[5] = PromptFraction(&waveform, results[0], offset, prebins, window);
      pf->SetDouble(results[5], p);
      results[6] = PseudoPromptFraction(&waveform, start, end, offset, nsamples, threshold);
      spf->SetDouble(results[6], p);
      
      // Do the two exponential fits and gaussian fits 
      Fit(&waveform, start, end, threshold, results);
      exp_amp->SetDouble(results[0], p);
      exp_t0->SetDouble(results[1]+t0->GetDouble(p), p);
      exp_rise->SetDouble(results[2], p);
      exp_fall->SetDouble(results[3], p);
      exp_chisq->SetDouble(results[4], p);
      exp_dof->SetDouble(results[5], p);
      gaus_amp->SetDouble(results[6], p);
      gaus_mu->SetDouble(results[7]+t0->GetDouble(p), p);
      gaus_sigma->SetDouble(results[8], p);
      gaus_chisq->SetDouble(results[9], p);
      gaus_dof->SetDouble(results[10], p);
      
      // Do the S2Filter 
      S2Filter(&waveform, S1_width, S2_width, results[0], results[1]);
      s2filter_s1->SetDouble(results[0], p);
      s2filter_s2->SetDouble(results[1], p);
      s2filter_diff->SetDouble(fabs(results[1]-results[0]), p);
      
      // peak rqs 
      results[9] = results[10] = results[11] = 0;
      for (size_t pp=0; pp<122; pp++) {
        ch = cvt->cvt_events[e]->channels[pp];
        PeakBaselineStats(ch, pstart->GetInt(p), pend->GetInt(p), results);
        if (results[1] == 0 and results[3] == 0) continue;
        
        // store peak baseline stats 
        pre_mean->SetDouble(results[0], p, pp);
        pre_std->SetDouble(results[1], p, pp);
        pos_mean->SetDouble(results[2], p, pp);
        pos_std->SetDouble(results[3], p, pp);
        baseline->SetDouble(results[4], p, pp);
        
        // calculate peak area and heights         
        gain = mVns_per_phe[pp]*preamp*postamp;
        waveform = ch->GetPeak(pstart->GetInt(p), pend->GetInt(p));
        PeakQuantities(&waveform, start, end, results);
        pparea->SetDouble(results[0], p, pp);
        ppheight->SetDouble(results[1], p, pp);
        ppheight_mV->SetDouble(results[1]*gain/10, p, pp);
        f_l_points->SetDouble(results[2], p, pp);
        
        // calculate saturation flags 
        if (results[1]*gain/10 >= daq_saturation_mV) {
          daq_sat->SetInt(1, pp);
        }
        if (results[1]*gain/10 >= pmt_saturation_mV) {
          pmt_sat->SetInt(1, pp);
        }
        
        // calculate quantities for top and bottom array 
        if (pp+1<61 or pp==121) {
          results[10] += results[0];
        }else {
          results[11] += results[0];
        }
        if (results[9] < results[0]) results[9] = results[0];
      }
      
      // calculate top and bottom array 
      top_area->SetDouble(results[10], p);
      bot_area->SetDouble(results[11], p);
      tb_ratio->SetDouble(results[10]/(results[11]+small_number), p);
      tb_asym->SetDouble((results[10]-results[11])/(results[10]+results[11]+small_number), p);
      max_peak_area->SetDouble(results[9], p);
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
