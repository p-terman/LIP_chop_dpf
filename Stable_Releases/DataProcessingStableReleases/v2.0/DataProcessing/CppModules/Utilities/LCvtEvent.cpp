#include "LCvtEvent.h"
#include "LEvtUtilities.h"
#include <iostream>
#include <algorithm>
#include <numeric>
#include <boost/property_tree/ptree.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>

#define VERBOSE 0

using namespace std;

LCvtEvent::LCvtEvent(){;}

//Assuming that &input is an open ifstream, instantiate the next event from it
LCvtEvent::LCvtEvent(std::ifstream &input, bool byteswap,boost::property_tree::ptree &settings_pt){
  input.readval(gid_date_time);
  input.readval(gid_location);
  input.readval(gid_event_num);
  input.readval(gid_num_chans);
  input.readval(gid_event_size);
//   if (settings_pt.get<bool>("daq_settings.sis3301.global.read_xlm")) {
// 		input.readval(ddc_trigger_timestamp);
// 		if(settings_pt.get<float>("daq_settings.global.daq_version")>7.0){
// 		  input.readval(ddc_sequence_number);
// 		}
// 		input.readval(ddc_max_filter_resp);
// 		input.readval(ddc_max_chan_ID);
// 		input.readval(num_ddc_boards);
// 		ddc_s1.resize(num_ddc_boards);
// 		ddc_s2.resize(num_ddc_boards);
// 		for(unsigned short bb=0;bb<num_ddc_boards;bb++){
// 		  input.readval(ddc_s1[bb]);
// 		  input.readval(ddc_s2[bb]);
// 		}
// 		if(settings_pt.get<float>("daq_settings.global.daq_version")>7.0){
// 		  input.readval(ddc_check_byte);
// 		}
// #if VERBOSE>0
//   cout<<num_ddc_boards<<endl;
// #endif
//   }
  input.readval(record_format);
  input.readval(record_size);
  
  input.readval(event_trigger_timestamp);
#if VERBOSE>0
  cout<<"trigger_timestamp"<<event_trigger_timestamp<<endl;
#endif
  
  for(unsigned int cc=0;cc<gid_num_chans; cc++){
    channels.push_back(new LCvtChannel(input, byteswap));
  }
}


LCvtEvent::LCvtEvent(LEvtEvent *ev) :
  gid_date_time(ev->gid_date_time),
  gid_location(ev->gid_location),
  gid_event_num(ev->gid_event_num),
  gid_num_chans(ev->gid_num_chans),
  gid_event_size(ev->gid_event_size),
  /*
  ddc_trigger_timestamp(ev->ddc_trigger_timestamp),
  ddc_trigger_timestamp(ev->ddc_sequence_number),
  ddc_max_filter_resp(ev->ddc_max_filter_resp),
  ddc_max_chan_ID(ev->ddc_max_chan_ID),
  num_ddc_boards(ev->num_ddc_boards),
  ddc_s1(ev->ddc_s1),
  ddc_s2(ev->ddc_s2),
  ddc_check_byte(ev->ddc_check_byte),
  */
  record_format(ev->record_format),
  record_size(ev->record_size),  
  event_trigger_timestamp(ev->event_trigger_timestamp)
{
  for(unsigned int cc=0;cc<gid_num_chans; cc++){
    channels.push_back(new LCvtChannel(ev->channels[cc]));
  }
}


LCvtEvent::~LCvtEvent(){
  for (std::vector<LCvtChannel*>::iterator iter = channels.begin();
    iter != channels.end(); ++iter) {
    delete *iter;
  }  
}

