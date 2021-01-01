#ifndef __RQS_H__
#define __RQS_H__ 1
/*
This Header has the MinimumSet Class which allows for quick loading of the rqs 
and their waveforms.

Author: Sergey Uvarov
Version 0.9 created at 2/27/2013

*/
#include "RQFile_IO.hh"
#include "LEvtEvent.h"
#include <vector>
#include <iostream>
#include <cmath>

using std::vector;
using std::cerr;
using std::sqrt;
using std::abs;
using std::cout;
using std::endl;

typedef std::vector<double> DVec;
typedef std::vector< vector<double> > D2Vec;
typedef std::vector< vector< vector<double> > > D3Vec;

inline void PODVariance(DVec *pod, double& pre_mean, double& pre_std, 
                  double& pos_mean, double& pos_std) {
  // 
  pre_mean = pre_std = pos_mean = pos_std = 0;
  int last = int(pod->size())-1;
  for (int i=0; i<24; i++) {
    pre_mean += pod->at(i);
    pre_std  += pod->at(i)*pod->at(i);
  }
  for (int i=0; i<31; i++) {
    pos_mean += pod->at(last-i);
    pos_std  += pod->at(last-i)*pod->at(last-i);
  }
  pre_std = sqrt(abs(pre_std-pre_mean*pre_mean/24)/23);
  pos_std = sqrt(abs(pos_std-pos_mean*pos_mean/31)/30);
  pre_mean /= 24;
  pos_mean /= 31;  
}
//______________________________________________________________________________
class RQSet {
  public:
    ////////////////////////////
    //////// Members ///////////
    ////////////////////////////    
    Rq* event_number;
    // pulse finder
    Rq* pulse_start_samples;
    Rq* pulse_end_samples;
    Rq* pulse_start_pod;
    Rq* pulse_end_pod;
    Rq* pulse_length_samples;
    // pulse timing
    Rq* t0_area_samples;
    Rq* t05_area_samples;
    Rq* t25_area_samples;
    Rq* t50_area_samples;
    Rq* t75_area_samples;
    Rq* t95_area_samples;
    Rq* t2_area_samples;
    Rq* t0_height_samples;
    Rq* t10l_height_samples;
    Rq* t50l_height_samples;
    Rq* t1_height_samples;
    Rq* t50r_height_samples;
    Rq* t10r_height_samples;
    Rq* t2_height_samples;
    // pulse quantities
    Rq* pulse_area_phe;        
    Rq* pulse_height_phe_per_sample; 
    Rq* prompt_fraction;
    Rq* pseudo_prompt_fraction;
    Rq* exp_fit_amplitude_phe_per_sample;
    Rq* exp_fit_tau_fall_samples;
    Rq* exp_fit_time_offset_samples;
    Rq* exp_fit_tau_rise_samples;
    Rq* exp_fit_chisq;
    Rq* exp_fit_dof;
    Rq* gaus_fit_amplitude_phe_per_sample;
    Rq* gaus_fit_mu_samples;
    Rq* gaus_fit_sigma_samples; 
    Rq* gaus_fit_chisq;
    Rq* gaus_fit_dof;
    Rq* peak_area_phe;
    Rq* peak_height_mV;
    Rq* peak_height_phe_per_sample;
    Rq* mean_first_last_pts_phe_per_sample;
    Rq* pod_detect_pre_pod_mean_phe;
    Rq* pod_detect_pos_pod_mean_phe;
    Rq* pod_detect_pre_pod_std_phe;
    Rq* pod_detect_pos_pod_std_phe;    
    Rq* baseline_daq_mV;
    Rq* daq_saturation_flag;
    Rq* pmt_2pct_saturation_flag;
    Rq* top_to_bottom_ratio;
    Rq* top_bottom_asymmetry;
    Rq* pulse_area_positive_phe;
    Rq* pulse_area_negative_phe;
    // pulse quality check
    Rq* pulse_quality_flag;
    // pulse classifier
    Rq* pulse_classification;
    // position reconstruction
    Rq* x_cm;
    Rq* y_cm;
    // waveforms
    D2Vec pulse_waveforms;
    D3Vec peak_waveforms;
    // extra 
    DVec pulse_thresholds;
    D2Vec peak_pre_mean_phe;
    D2Vec peak_pos_mean_phe;
    D2Vec peak_pre_std_phe;
    D2Vec peak_pos_std_phe;
    D2Vec baseline_daq;
    
    ////////////////////////////
    //////// Methods ///////////
    ////////////////////////////
    inline void GetRQs(RQBlock *block);
    inline void BuildWaveforms(LEvtEvent *event);
};
//______________________________________________________________________________
inline void RQSet::GetRQs(RQBlock *block) {
  // Assigns the RQs from the rq file to the RQSet members. If the rq does
  // not exist in the rq file, the variable will be a null pointer. Trying to 
  // access anything from a variable with a null pointer will cause a seg fault. 
  
  event_number = block->Get("event_number");
  // pulse finder 
  pulse_start_samples = block->Get("pulse_start_samples");
  pulse_end_samples = block->Get("pulse_end_samples");
  pulse_start_pod = block->Get("pulse_start_pod");
  pulse_end_pod = block->Get("pulse_end_pod");
  pulse_length_samples = block->Get("pulse_length_samples");
  // pulse timing 
  t0_area_samples = block->Get("t0_area_samples");
  t05_area_samples = block->Get("t05_area_samples");
  t25_area_samples = block->Get("t25_area_samples");
  t50_area_samples = block->Get("t50_area_samples");
  t75_area_samples = block->Get("t75_area_samples");
  t95_area_samples = block->Get("t95_area_samples");
  t2_area_samples = block->Get("t2_area_samples");
  t0_height_samples = block->Get("t0_height_samples");
  t10l_height_samples = block->Get("t10l_height_samples");
  t50l_height_samples = block->Get("t50l_height_samples");
  t1_height_samples = block->Get("t1_height_samples");
  t50r_height_samples = block->Get("t50r_height_samples");
  t10r_height_samples = block->Get("t10r_height_samples");
  t2_height_samples = block->Get("t2_height_samples");
  // pulse quantities
  pulse_area_phe = block->Get("pulse_area_phe");  
  pulse_height_phe_per_sample = block->Get("pulse_height_phe_per_sample");
  prompt_fraction = block->Get("prompt_fraction");
  pseudo_prompt_fraction = block->Get("pseudo_prompt_fraction");   
  exp_fit_amplitude_phe_per_sample = block->Get("exp_fit_amplitude_phe_per_sample");
  exp_fit_tau_fall_samples = block->Get("exp_fit_tau_fall_samples");
  exp_fit_time_offset_samples = block->Get("exp_fit_time_offset_samples");
  exp_fit_tau_rise_samples = block->Get("exp_fit_tau_rise_samples");
  exp_fit_chisq = block->Get("exp_fit_chisq");
  exp_fit_dof = block->Get("exp_fit_dof");
  gaus_fit_amplitude_phe_per_sample = block->Get("gaus_fit_amplitude_phe_per_sample");
  gaus_fit_mu_samples = block->Get("gaus_fit_mu_samples");
  gaus_fit_sigma_samples = block->Get("gaus_fit_sigma_samples");
  gaus_fit_chisq = block->Get("gaus_fit_chisq");
  gaus_fit_dof = block->Get("gaus_fit_dof");
  peak_area_phe = block->Get("peak_area_phe");
  peak_height_mV = block->Get("peak_height_mV");
  peak_height_phe_per_sample = block->Get("peak_height_phe_per_sample");
  mean_first_last_pts_phe_per_sample = block->Get("mean_first_last_pts_phe_per_sample");
  pod_detect_pre_pod_mean_phe = block->Get("pod_detect_pre_pod_mean_phe");
  pod_detect_pos_pod_mean_phe = block->Get("pod_detect_pos_pod_mean_phe");
  pod_detect_pre_pod_std_phe = block->Get("pod_detect_pre_pod_std_phe");
  pod_detect_pos_pod_std_phe = block->Get("pod_detect_pos_pod_std_phe");
  baseline_daq_mV = block->Get("baseline_daq_mV");
  daq_saturation_flag = block->Get("daq_saturation_flag");
  pmt_2pct_saturation_flag = block->Get("pmt_2pct_saturation_flag");
  top_to_bottom_ratio = block->Get("top_to_bottom_ratio");
  top_bottom_asymmetry = block->Get("top_bottom_asymmetry"); 
  pulse_area_positive_phe = block->Get("pulse_area_positive_phe");
  pulse_area_negative_phe = block->Get("pulse_area_negative_phe");
  // pulse quality check 
  pulse_quality_flag = block->Get("pulse_quality_flag");
  // pulse classifier
  pulse_classification = block->Get("pulse_classification");
  // position reconstruction 
  x_cm = block->Get("x_cm");
  y_cm = block->Get("y_cm");
}
//______________________________________________________________________________
inline void RQSet::BuildWaveforms(LEvtEvent *event) {
  // Builds the pulse waveforms for a given event. 
  // ASSUMTIONS: 
  //  "pulse_start_samples" and "pulse_end_samples" exist and already loaded
  //  "pulse_start_samples" and "pulse_end_samples" are the same dimension
  //  "pulse_start_pod" and "pulse_end_pod" exist and already loaded
  //  "pulse_start_pod" and "pulse_end_pod" are the same dimension 
  //  "pulse_start_pod" and "pulse_end_pod" indicies start at 1 not 0 
  //  Xe PMTs are mapped 1-1 in data index channels 0:121
  const size_t nch = 122; // number of Xe pmts
  pulse_waveforms.clear();
  peak_waveforms.clear();
  pulse_thresholds.clear();
  peak_pre_mean_phe.clear();
  peak_pos_mean_phe.clear();
  peak_pre_std_phe.clear();
  peak_pos_std_phe.clear();
  baseline_daq.clear();
  peak_pre_mean_phe.resize(pulse_start_samples->GetDim1(), DVec(nch,0));
  peak_pos_mean_phe.resize(pulse_start_samples->GetDim1(), DVec(nch,0));
  peak_pre_std_phe.resize(pulse_start_samples->GetDim1(), DVec(nch,0));
  peak_pos_std_phe.resize(pulse_start_samples->GetDim1(), DVec(nch,0));
  baseline_daq.resize(pulse_start_samples->GetDim1(), DVec(nch,0));
  pulse_thresholds.resize(pulse_start_samples->GetDim1(), 0);
  
  if (!pulse_start_samples or !pulse_end_samples) {
    cerr << "Error: Cannot Build Waveforms, pulse_start_samples or "
         << "pulse_end_samples does not exist\n";
    return;
  }
  if (!pulse_start_pod or !pulse_end_pod) {
    cerr << "Error: Cannot Build Waveforms, pulse_start_pod or "
         << "pulse_end_pod does not exist\n";
    return;
  }
  
  int start_time, end_time, pod_start, pod_end, pod_start_time, pod_end_time;
  DVec p_wave, pp_wave, *pod_waveform;
  double pre_mean, pre_std, pos_mean, pos_std;
  peak_waveforms.resize(pulse_start_samples->GetDim1(), D2Vec(nch));
  
  for (size_t i=0; i<pulse_start_samples->GetDim1(); i++) {    
    start_time = pulse_start_samples->GetInt(i);
    end_time = pulse_end_samples->GetInt(i);
    if (end_time - start_time <= 0.5) continue; // skip empty pulses    
    p_wave.resize(end_time - start_time + 1, 0);
    
    // loop throught the Xe PMTs only
    for (size_t c=0; c<nch; c++) {
      pod_start = pulse_start_pod->GetInt(i, c)-1;
      pod_end = pulse_end_pod->GetInt(i, c)-1;
      if (pod_start < 0) continue;
      
      pp_wave.resize(end_time - start_time + 1, 0);
      for (int p=pod_start; p<=pod_end; p++) {
        pod_start_time = event->channels[c]->pod_starts[p];
        pod_end_time = int(event->channels[c]->pod_data[p].size())+pod_start_time-1;
        if (pod_start_time > end_time or pod_end_time < start_time) continue;
        pod_waveform = &event->channels[c]->pod_data_phe[p];
        
        for (int j=pod_start_time; j<=pod_end_time and j<=end_time; j++) {
          if (j-start_time < 0) continue;
          p_wave[j-start_time] += pod_waveform->at(j-pod_start_time);
          pp_wave[j-start_time] = pod_waveform->at(j-pod_start_time);
        }
        PODVariance(pod_waveform, pre_mean, pre_std, pos_mean, pos_std);
        peak_pre_mean_phe[i][c] += pre_mean;
        peak_pos_mean_phe[i][c] += pos_mean;
        peak_pre_std_phe[i][c] += pre_std;
        peak_pos_std_phe[i][c] += pos_std;
        baseline_daq[i][c] += event->channels[c]->pod_baselines[p]*2000./16384.;
      }
      peak_pre_mean_phe[i][c] /= pod_end-pod_start+1;
      peak_pos_mean_phe[i][c] /= pod_end-pod_start+1;
      peak_pre_std_phe[i][c] /= pod_end-pod_start+1;
      peak_pos_std_phe[i][c] /= pod_end-pod_start+1;
      baseline_daq[i][c] /= pod_end-pod_start+1;
      peak_waveforms[i][c] = pp_wave;
      pp_wave.clear();
      
      pulse_thresholds[i] += peak_pre_std_phe[i][c]*peak_pre_std_phe[i][c];
    }
    
    pulse_thresholds[i] = sqrt(pulse_thresholds[i]);
    pulse_waveforms.push_back(p_wave);
    p_wave.clear();
  }
}
#endif
