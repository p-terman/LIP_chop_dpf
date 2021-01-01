#ifndef LEVTCHANNEL_H
#define LEVTCHANNEL_H
#include "LEvtUtilities.h"
#include <vector>
#include <fstream>

class LEvtChannel {
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
    std::vector<unsigned int> pod_baselines;
    std::vector<std::vector<short> > pod_data;
    std::vector<std::vector<double> > pod_data_phe;
    
    ////////////////////////////
    //////// Methods ///////////
    ////////////////////////////
    
    LEvtChannel();
    LEvtChannel(std::ifstream &input, bool byteswap);
    int Read(std::ifstream &input, bool byteswap, size_t file_end_position);
    void CalibratePods(double mVns_per_phe,double amps_gain=7.5,double sample_interval_ns=10.0,double full_scale_mV=2000.0, double adc_bins=16384.0);
};

#endif
