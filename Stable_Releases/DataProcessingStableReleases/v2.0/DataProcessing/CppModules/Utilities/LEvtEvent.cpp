#include "LEvtEvent.h"
#include "LEvtUtilities.h"
#include <iostream>
#include <algorithm>
#include <numeric>
#include <exception>
#include <boost/property_tree/ptree.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>

#define VERBOSE 0

using boost::numeric::ublas::mapped_vector;
using namespace std;

LEvtEvent::LEvtEvent(){;}


LEvtEvent::LEvtEvent(size_t nchannels) :
gid_date_time(0),
gid_location(0),
gid_event_num(0),
gid_num_chans(nchannels),
gid_event_size(0),
ddc_trigger_timestamp(0),
ddc_sequence_number(0),
ddc_max_filter_resp(0),
ddc_max_chan_ID(0),
num_ddc_boards(0),
ddc_check_byte(0),
record_format(0),
record_size(0),
event_trigger_timestamp(0)
{
	for(unsigned int cc=0;cc</*1*/gid_num_chans; cc++){
    channels.push_back(new LEvtChannel());
  }
}


//Assuming that &input is an open ifstream, instantiate the next event from it
LEvtEvent::LEvtEvent(std::ifstream &input, bool byteswap,boost::property_tree::ptree &settings_pt){
  input.readval(gid_date_time);
  input.readval(gid_location);
  input.readval(gid_event_num);
  input.readval(gid_num_chans);
  input.readval(gid_event_size);
  
  // 20130528 - SU : changed from get to get_optional in case of reading a sim file 
  bool read_xlm = false;
  try {
    read_xlm = settings_pt.get<bool>("daq_settings.sis3301.global.read_xlm");
  } catch (std::exception const& ex) { read_xlm = false; }
  float version = 6;
  try {
    version = settings_pt.get<float>("daq_settings.global.daq_version");
  } catch(std::exception const& ex) { version = 6; }
  
  if (read_xlm) {
		input.readval(ddc_trigger_timestamp);
		if(version > 7.0){
		  input.readval(ddc_sequence_number);
		}
		input.readval(ddc_max_filter_resp);
		input.readval(ddc_max_chan_ID);
		input.readval(num_ddc_boards);
		ddc_s1.resize(num_ddc_boards);
		ddc_s2.resize(num_ddc_boards);
		for(unsigned short bb=0;bb<num_ddc_boards;bb++){
		  input.readval(ddc_s1[bb]);
		  input.readval(ddc_s2[bb]);
		}
		if(version > 7.0){
		  input.readval(ddc_check_byte);
		}
  }
  input.readval(record_format);
  input.readval(record_size);
  
  input.readval(event_trigger_timestamp);
  
  for(unsigned int cc=0;cc</*1*/gid_num_chans; cc++){
    channels.push_back(new LEvtChannel(input, byteswap));
  }
}

int LEvtEvent::Read(std::ifstream &input, bool byteswap,boost::property_tree::ptree &settings_pt, size_t file_end_position){
  input.readval(gid_date_time);
  input.readval(gid_location);
  input.readval(gid_event_num);
  input.readval(gid_num_chans);
  input.readval(gid_event_size);
  
  // 20130528 - SU : changed from get to get_optional in case of reading a sim file 
  bool read_xlm = false;
  try {
    read_xlm = settings_pt.get<bool>("daq_settings.sis3301.global.read_xlm");
  } catch (std::exception const& ex) { read_xlm = false; }
  float version = 6;
  try {
    version = settings_pt.get<float>("daq_settings.global.daq_version");
  } catch(std::exception const& ex) { version = 6; }
  
  if (read_xlm) {
		input.readval(ddc_trigger_timestamp);
		if(version > 7.0){
		  input.readval(ddc_sequence_number);
		}
		input.readval(ddc_max_filter_resp);
		input.readval(ddc_max_chan_ID);
		input.readval(num_ddc_boards);
		ddc_s1.resize(num_ddc_boards);
		ddc_s2.resize(num_ddc_boards);
		for(unsigned short bb=0;bb<num_ddc_boards;bb++){
		  input.readval(ddc_s1[bb]);
		  input.readval(ddc_s2[bb]);
		}
		if(version > 7.0){
		  input.readval(ddc_check_byte);
		}
  }
  input.readval(record_format);
  input.readval(record_size);
  
  input.readval(event_trigger_timestamp);
  
  LEvtChannel *ch = 0;
  int status = 0;
  size_t current_position;
  for(unsigned int cc=0;cc</*1*/gid_num_chans; cc++){
  	current_position = input.tellg();
  	if (file_end_position < current_position) return -1;
    ch = new LEvtChannel();
    status = ch->Read(input, byteswap, file_end_position);
    if (status == -1) {
    	delete ch;
    	return -1;
    }
    channels.push_back(ch);
  }
  return 0;
}


LEvtEvent::~LEvtEvent(){
  for (std::vector<LEvtChannel*>::iterator iter = channels.begin();
    iter != channels.end(); ++iter) {
    delete *iter;
  }  
}

