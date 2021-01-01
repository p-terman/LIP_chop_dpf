#ifndef LCVTCHANNEL_H
#define LCVTCHANNEL_H
#include "LEvtUtilities.h"
#include "LEvtChannel.h"
#include <vector>
#include <fstream>

class LCvtChannel {
  public:
    ////////////////////////////
    //////// Members ///////////
    ////////////////////////////
    char ch_binary_datatype;
    double voltage_resolution;
    double voltage_offset;
    double time_resolution;
    int pretrigger;
    unsigned int event_size;
    unsigned int pulse_detect_pretrigger;
    unsigned int pulse_detect_posttrigger;
    unsigned int number_of_pods;

    std::vector<int> pod_starts;
    std::vector<unsigned int> pod_lengths;
    std::vector<float> pod_baselines;
    std::vector<std::vector<float> > pod_data_phe;
    
    ////////////////////////////
    //////// Methods ///////////
    ////////////////////////////
    
    LCvtChannel();
    LCvtChannel(std::ifstream &input, bool byteswap);
    LCvtChannel(LEvtChannel *ch);
    inline size_t GetSize();
    inline void Convert(float amp);
    inline std::vector<float> GetPeak(int start, int end);
};

inline size_t LCvtChannel::GetSize() {
  size_t size = 0;
  size += sizeof(ch_binary_datatype);
  size += sizeof(voltage_resolution);
  size += sizeof(voltage_offset);
  size += sizeof(time_resolution);
  size += sizeof(pretrigger);
  size += sizeof(event_size);
  size += sizeof(pulse_detect_pretrigger);
  size += sizeof(pulse_detect_posttrigger);
  size += sizeof(number_of_pods);
  size += sizeof(int)*number_of_pods;//sizeof(pod_starts[0])*number_of_pods;
  size += sizeof(unsigned int)*number_of_pods;//sizeof(pod_lengths[0])*number_of_pods;
  size += sizeof(float)*number_of_pods;//sizeof(pod_baselines[0])*number_of_pods;
  for (unsigned int pp=0; pp<number_of_pods; pp++) {
   size += sizeof(float)*pod_lengths[pp];//sizeof(pod_data_phe[pp][0])*number_of_pods;
  }
  return size;
}

inline void LCvtChannel::Convert(float amp) {
  float ADC_to_mV = 2000.0/16384.0;
  float time_resolution = 10.0;  
  for (size_t pp=0; pp<pod_data_phe.size(); pp++) {
    for (size_t i=0; i<pod_data_phe[pp].size(); i++) {
      pod_data_phe[pp][i] = (pod_baselines[pp] - pod_data_phe[pp][i])
                          * time_resolution * ADC_to_mV / amp;
    }
    //pod_baselines[pp] *= ADC_to_mV;
    pod_baselines[pp] = 0;
  }
}

inline std::vector<float> LCvtChannel::GetPeak(int pstart, int pend) {
  std::vector<float> peak(pend-pstart+1, 0);
  int pod_start, pod_end;
  for (size_t pp=0; pp<pod_data_phe.size(); pp++) {
    // select sumpods that overlap pulse times 
    pod_start = pod_starts[pp];
    pod_end = pod_start+int(pod_lengths[pp])-1;
    if (pod_start > pend or pod_end < pstart) continue;
    // assign sumpod values to pulse 
    for (int i=pod_start; i<=pod_end; i++) {
      if (i-pstart < 0) continue;
      if (i > pend) break;
      peak[i-pstart] = pod_data_phe[pp][i-pod_start];
    }    
  }  
  return peak;
}


#endif
