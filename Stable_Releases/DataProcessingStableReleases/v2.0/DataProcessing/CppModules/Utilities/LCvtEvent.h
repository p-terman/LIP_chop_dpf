#ifndef LCVTEVENT_H
#define LCVTEVENT_H
#include <vector>
#include <fstream>
#include <boost/property_tree/ptree.hpp>
#include "LEvtUtilities.h"
#include "LEvtEvent.h"
#include "LCvtChannel.h"

using std::vector;

class LCvtEvent {
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
    
    std::vector<LCvtChannel*> channels;
    

    ////////////////////////////
    //////// Methods ///////////
    ////////////////////////////

    LCvtEvent();
    
    //Assuming that &input is an open ifstream, instantiate the next event from it
    LCvtEvent(std::ifstream &input, bool byteswap, boost::property_tree::ptree &settings_pt);
    LCvtEvent(LEvtEvent *ev);
    ~LCvtEvent();
    inline size_t GetSize();
    inline void AddXeSumPodCh(vector<int> start, vector<int> end);
    inline void AddWaSumPodCh(vector<int> start, vector<int> end);
    
};

inline size_t LCvtEvent::GetSize() {
  size_t size = 0;
  size += sizeof(gid_date_time);
  size += sizeof(gid_location);
  size += sizeof(gid_event_num);
  size += sizeof(gid_num_chans);
  size += sizeof(gid_event_size);
  /*
  // The xlm trigger information is not transfered over to the cvt file 
  size += sizeof(ddc_trigger_timestamp);
  size += sizeof(ddc_sequence_number);
  size += sizeof(ddc_max_filter_resp);
  size += sizeof(ddc_max_chan_ID);
  size += sizeof(num_ddc_boards);
  size += sizeof(char)*num_ddc_boards;//sizeof(ddc_s1[0])*num_ddc_boards;
  size += sizeof(char)*num_ddc_boards;//sizeof(ddc_s2[0])*num_ddc_boards;
  size += sizeof(ddc_check_byte);
  */
  size += sizeof(record_format);
  size += sizeof(record_size);
  size += sizeof(event_trigger_timestamp);
  
  for(unsigned int cc=0;cc<gid_num_chans; cc++){
    size += channels[cc]->GetSize();
  }
  return size;
}

inline void LCvtEvent::AddXeSumPodCh(vector<int> start, vector<int> end) {
  // containers
  LCvtChannel *ch;
  int podstart, podend;
  // check if sumpod channel exist 
  if (channels.size() < 137) {
    channels.resize(137, NULL);
    for (size_t i=0; i<137; i++) {
      if (!channels[i]) channels[i] = new LCvtChannel();
    }    
    gid_num_chans = 137;
  }
  ch = channels[136];
  // assign values 
  ch->ch_binary_datatype = 0;
  ch->voltage_resolution = 2000.0/16384.0;
  ch->voltage_offset = 0;
  ch->time_resolution = 1e-8;
  ch->pretrigger = 0;
  ch->event_size = 0;
  ch->pulse_detect_pretrigger = 24;
  ch->pulse_detect_posttrigger = 31;
  ch->number_of_pods = start.size();  //  <-- very important
  ch->pod_starts.resize(ch->number_of_pods);
  ch->pod_lengths.resize(ch->number_of_pods);
  ch->pod_baselines.resize(ch->number_of_pods);
  ch->pod_data_phe.resize(ch->number_of_pods);
  // trasfer the start and end data 
  for (unsigned int i=0; i<ch->number_of_pods; i++) {
    ch->pod_starts[i] = start[i];
    ch->pod_lengths[i] = end[i]-start[i]+1;
    ch->pod_baselines[i] = 0;
    ch->pod_data_phe[i].clear(); // remove any old data 
    ch->pod_data_phe[i].resize(ch->pod_lengths[i], 0); 
    // loop through Xe channels
    for (size_t cc=0; cc<122; cc++) {
      for (size_t pp=0; pp<channels[cc]->pod_starts.size(); pp++) {
        // select proper pods 
        podstart = channels[cc]->pod_starts[pp];
        podend = podstart+int(channels[cc]->pod_lengths[pp])-1;
        if (podstart > end[i] or podend < start[i]) continue;
        
        // add pod to sumpod
        for (int j=podstart; j<=podend; j++) {
          if (j < start[i]) continue;
          if (j > end[i]) break; 
          ch->pod_data_phe[i][j-start[i]] += channels[cc]->pod_data_phe[pp][j-podstart];
        }
      }      
    } // end of Xe channel loop 
  }
}

inline void LCvtEvent::AddWaSumPodCh(vector<int> start, vector<int> end) {
  // containers
  LCvtChannel *ch;
  int podstart, podend;
  // check if sumpod channel exist 
  if (channels.size() < 138) {
    channels.resize(138, NULL);
    for (size_t i=0; i<138; i++) {
      if (!channels[i]) channels[i] = new LCvtChannel();
    } 
    gid_num_chans = 138;
  }
  ch = channels[137];
  // assign values 
  ch->ch_binary_datatype = 0;
  ch->voltage_resolution = 2000.0/16384.0;
  ch->voltage_offset = 0;
  ch->time_resolution = 1e-8;
  ch->pretrigger = 0;
  ch->event_size = 0;
  ch->pulse_detect_pretrigger = 24;
  ch->pulse_detect_posttrigger = 31;
  ch->number_of_pods = start.size();  //  <-- very important
  ch->pod_starts.resize(ch->number_of_pods);
  ch->pod_lengths.resize(ch->number_of_pods);
  ch->pod_baselines.resize(ch->number_of_pods);
  ch->pod_data_phe.resize(ch->number_of_pods);
  // trasfer the start and end data 
  for (unsigned int i=0; i<ch->number_of_pods; i++) {
    ch->pod_starts[i] = start[i];
    ch->pod_lengths[i] = end[i]-start[i]+1;
    ch->pod_baselines[i] = 0;
    ch->pod_data_phe[i].clear(); // remove any old data 
    ch->pod_data_phe[i].resize(ch->pod_lengths[i], 0); 
    // loop through Water channels
    for (size_t cc=128; cc<136; cc++) {
      for (size_t pp=0; pp<channels[cc]->pod_starts.size(); pp++) {
        // select proper pods 
        podstart = channels[cc]->pod_starts[pp];
        podend = podstart+int(channels[cc]->pod_lengths[pp])-1;
        if (podstart > end[i] or podend < start[i]) continue;
        
        // add pod to sumpod
        for (int j=podstart; j<=podend; j++) {
          if (j < start[i]) continue;
          if (j > end[i]) break; 
          ch->pod_data_phe[i][j-start[i]] += channels[cc]->pod_data_phe[pp][j-podstart];
        }
      }      
    } // end of Water channel loop 
  }
}


#endif
