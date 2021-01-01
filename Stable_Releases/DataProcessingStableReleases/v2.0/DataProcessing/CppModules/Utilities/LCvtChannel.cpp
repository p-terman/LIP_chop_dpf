#include <algorithm>
#include <functional>
#include <iostream>

#include "LCvtChannel.h"
#include "LEvtUtilities.h"

LCvtChannel::LCvtChannel() :
  ch_binary_datatype(0),
  voltage_resolution(2000.0/16384.0),
  voltage_offset(0),
  time_resolution(1e-8),
  pretrigger(0),
  event_size(0),
  pulse_detect_pretrigger(24),
  pulse_detect_posttrigger(31),
  number_of_pods(0)
{;}

//Assuming that &input is an open ifstream, instantiate the next channel from it
LCvtChannel::LCvtChannel(std::ifstream &input, bool byteswap){
  input.readval(ch_binary_datatype);
//  std::cout<<"ch_binary_datatype:"<<ch_binary_datatype<<std::endl;
  input.readval(voltage_resolution);
  input.readval(voltage_offset);
  input.readval(time_resolution);
//  std::cout<<"time_resolution:"<<time_resolution<<std::endl;
  input.readval(pretrigger);
  input.readval(event_size);
  input.readval(pulse_detect_pretrigger);
  input.readval(pulse_detect_posttrigger);
//  std::cout<<"pulse_detect_posttrigger:"<<pulse_detect_posttrigger<<std::endl;
  input.readval(number_of_pods);
//  std::cout<<"number_of_pods:"<<number_of_pods<<std::endl;
  if (number_of_pods==0) return;
  pod_starts.resize(number_of_pods);
  input.readvec(pod_starts);
//  std::cout<<"pod_start:"<<pod_starts.back()<<std::endl;
  pod_lengths.resize(number_of_pods);
  input.readvec(pod_lengths);
//  std::cout<<"pod_length:"<<pod_lengths.back()<<std::endl;
  pod_baselines.resize(number_of_pods);
  input.readvec(pod_baselines);
  
  // set all pod_baselines to 0
  pod_baselines.clear();
  pod_baselines.resize(number_of_pods, 0);
  
  pod_data_phe.resize(number_of_pods);
  for(unsigned short pp=0;pp<number_of_pods;pp++){
    pod_data_phe.at(pp).resize(pod_lengths.at(pp));
    char* addr= reinterpret_cast<char*>(&pod_data_phe.at(pp).front());
    input.read(addr,sizeof(pod_data_phe.at(pp).front())*pod_data_phe.at(pp).size());
  }
}

LCvtChannel::LCvtChannel(LEvtChannel *ch) :   
  ch_binary_datatype(ch->ch_binary_datatype),
  voltage_resolution(ch->voltage_resolution),
  voltage_offset(ch->voltage_offset),
  time_resolution(ch->time_resolution),
  pretrigger(ch->pretrigger),
  event_size(ch->event_size),
  pulse_detect_pretrigger(ch->pulse_detect_pretrigger),
  pulse_detect_posttrigger(ch->pulse_detect_posttrigger),
  number_of_pods(ch->number_of_pods),
  pod_starts(ch->pod_starts),
  pod_lengths(ch->pod_lengths)
{
  // pod data has to converted the longer way because ch->pod_data_phe is a 
  // 2D short array while pod_data_phe is a 2D float array.
  pod_data_phe.resize(number_of_pods);
  pod_baselines.resize(number_of_pods);
  for (unsigned int pp=0; pp<number_of_pods; pp++) {
    pod_baselines[pp] = ch->pod_baselines[pp];
    pod_data_phe[pp].resize(pod_lengths[pp]);
    for (unsigned int i=0; i<pod_lengths[pp]; i++)
      pod_data_phe[pp][i] = ch->pod_data[pp][i];
  }
}
