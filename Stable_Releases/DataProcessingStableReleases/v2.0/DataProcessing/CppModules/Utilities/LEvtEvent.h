#ifndef LEVTEVENT_H
#define LEVTEVENT_H
#include <vector>
#include <fstream>
#include <boost/property_tree/ptree.hpp>
#include "LEvtUtilities.h"
#include "LEvtChannel.h"

class LEvtEvent {
  public:
    ////////////////////////////
    //////// Members ///////////
    ////////////////////////////
    
    unsigned int gid_date_time;
    unsigned int gid_location;
    unsigned int gid_event_num;
    unsigned int gid_num_chans;
    unsigned int gid_event_size;
    unsigned long long ddc_trigger_timestamp;
    unsigned int ddc_sequence_number;
    unsigned int ddc_max_filter_resp;
    char ddc_max_chan_ID;
    unsigned short num_ddc_boards;
    std::vector<char> ddc_s1;
    std::vector<char> ddc_s2;
    char ddc_check_byte;
    
    unsigned int record_format;
    unsigned int record_size;
    
    unsigned long long event_trigger_timestamp;
    
    std::vector<LEvtChannel*> channels;
    
    int first_sample;
    int last_sample;
    std::vector<double> sum_phe;
    std::vector<std::vector<double> > sumpod_data_phe; // to be implemented

    ////////////////////////////
    //////// Methods ///////////
    ////////////////////////////

    LEvtEvent();
    LEvtEvent(size_t nchannels);
    
    //Assuming that &input is an open ifstream, instantiate the next event from it
    LEvtEvent(std::ifstream &input, bool byteswap, boost::property_tree::ptree &settings_pt);
    int Read(std::ifstream &input, bool byteswap, boost::property_tree::ptree &settings_pt, size_t file_end_position);
    ~LEvtEvent();
        
};

#endif
