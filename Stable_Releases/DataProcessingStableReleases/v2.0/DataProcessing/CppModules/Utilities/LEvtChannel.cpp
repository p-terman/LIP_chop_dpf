#include <algorithm>
#include <functional>
#include <iostream>

#include "LEvtChannel.h"
#include "LEvtUtilities.h"

LEvtChannel::LEvtChannel() : 
ch_binary_datatype(0),
voltage_resolution(1),
voltage_offset(0),
time_resolution(1),
pretrigger(0),
event_size(0),
pulse_detect_pretrigger(0),
pulse_detect_posttrigger(0),
number_of_pods(0)
{}

//Assuming that &input is an open ifstream, instantiate the next channel from it
LEvtChannel::LEvtChannel(std::ifstream &input, bool byteswap){
  input.readval(ch_binary_datatype);
  //std::cout<<"ch_binary_datatype:"<<ch_binary_datatype<<std::endl;
  input.readval(voltage_resolution);
  input.readval(voltage_offset);
  input.readval(time_resolution);
  input.readval(pretrigger);
  input.readval(event_size);
  input.readval(pulse_detect_pretrigger);
  input.readval(pulse_detect_posttrigger);
  //std::cout<<"pulse_detect_posttrigger:"<<pulse_detect_posttrigger<<std::endl;
  input.readval(number_of_pods);
  //std::cout<<"number_of_pods:"<<number_of_pods<<std::endl;
  
  pod_starts.resize(number_of_pods);
  input.readvec(pod_starts);
  //std::cout<<"pod_start:"<<pod_starts.back()<<std::endl;
  pod_lengths.resize(number_of_pods);
  input.readvec(pod_lengths);
  //std::cout<<"pod_length:"<<pod_lengths.back()<<std::endl;
  pod_baselines.resize(number_of_pods);
  input.readvec(pod_baselines);
  pod_data.resize(number_of_pods);
  
  for(unsigned short pp=0;pp<number_of_pods;pp++){
    pod_data.at(pp).resize(pod_lengths.at(pp));
    char* addr= reinterpret_cast<char*>(&pod_data.at(pp).front());
    input.read(addr,sizeof(pod_data.at(pp).front())*pod_data.at(pp).size());
  }
}

int LEvtChannel::Read(std::ifstream &input, bool byteswap, size_t file_end_position){
  input.readval(ch_binary_datatype);
  //std::cout<<"ch_binary_datatype:"<<ch_binary_datatype<<std::endl;
  input.readval(voltage_resolution);
  input.readval(voltage_offset);
  input.readval(time_resolution);
  input.readval(pretrigger);
  input.readval(event_size);
  input.readval(pulse_detect_pretrigger);
  input.readval(pulse_detect_posttrigger);
  //std::cout<<"pulse_detect_posttrigger:"<<pulse_detect_posttrigger<<std::endl;
  input.readval(number_of_pods);
  //std::cout<<"number_of_pods:"<<number_of_pods<<std::endl;
  
  pod_starts.resize(number_of_pods);
  input.readvec(pod_starts);
  //std::cout<<"pod_start:"<<pod_starts.back()<<std::endl;
  pod_lengths.resize(number_of_pods);
  input.readvec(pod_lengths);
  //std::cout<<"pod_length:"<<pod_lengths.back()<<std::endl;
  pod_baselines.resize(number_of_pods);
  input.readvec(pod_baselines);
  pod_data.resize(number_of_pods);
  
  size_t current_position;
  for(unsigned short pp=0;pp<number_of_pods;pp++){
    current_position = input.tellg();
    if (file_end_position < current_position+pod_lengths.at(pp)) return -1;
    
    pod_data.at(pp).resize(pod_lengths.at(pp));
    char* addr= reinterpret_cast<char*>(&pod_data.at(pp).front());
    input.read(addr,sizeof(pod_data.at(pp).front())*pod_data.at(pp).size());
  }
  return 0;
}

void LEvtChannel::CalibratePods(double mVns_per_phe,double electronics_gain, double sample_interval_ns, double full_scale_mV, double adc_bins){
  double phe_per_mVsample=-1.0*sample_interval_ns*full_scale_mV/(electronics_gain*mVns_per_phe*adc_bins);
  pod_data_phe.resize(pod_data.size());
  for(unsigned int pp=0;pp<number_of_pods;pp++){
    pod_data_phe.at(pp).resize(pod_data.at(pp).size());
    std::transform(pod_data[pp].begin(),pod_data[pp].end(),pod_data_phe[pp].begin(),std::bind2nd(std::minus<double>(),pod_baselines[pp]));
    std::transform(pod_data_phe[pp].begin(),pod_data_phe[pp].end(),pod_data_phe[pp].begin(),std::bind2nd(std::multiplies<double>(),phe_per_mVsample));
  }
}
