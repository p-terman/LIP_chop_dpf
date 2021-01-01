#include <iostream>
#include <exception>
#include "LEvtFile.h"
#include "LEvtChannel.h"
using namespace std;
LEvtFile::LEvtFile(){}

LEvtFile::LEvtFile(std::string path, std::string evt_file){
  byteswap = false;
  std::string fullpath=path+'/'+evt_file;
  input.open(fullpath.data(), std::ios_base::binary);
  if(!input.is_open()){
    std::cerr<<"Failed to open file: "<<fullpath.data()<<std::endl;
  }
  else{
    std::cout<<fullpath.data()<<std::endl;
    this->GetSettings();
    this->GetFileHeader();
    this->GetLivetime();
    this->GetEvents();
    input.close();
  }
}

LEvtFile::~LEvtFile(){
  for (std::vector<LEvtEvent*>::iterator iter = evt_events.begin();
    iter != evt_events.end(); ++iter) {
    delete *iter;
  }
  delete [] settings_string;
}

//Assuming input is an open fstream, check endianness then read settings_string
void LEvtFile::GetSettings() {
  input.readval(endianness);
  byteswap = false;
  if(endianness == 0x04030201)
    byteswap = true;
	
  input.readval(settings_string_len);
  settings_string = new char[settings_string_len+1];
  input.readstr(settings_string, settings_string_len);
  settings_string[settings_string_len] = 0;
  
  //make an sstream, which is an istream, to use the boost xml parser
  std::stringstream settings_sstream(std::stringstream::in | std::stringstream::out);
  settings_sstream<<settings_string;
  read_xml((std::istream&) settings_sstream, settings_pt);
}

//Assuming you've just read in the settings_string, read file header
void LEvtFile::GetFileHeader(){
  input.readval(fh_endianness);
  //if(fh_endianness != 0x01020304)
  //  byteswap = true;
  input.readval(fh_date_time);
  input.readval(fh_location);
  input.readval(fh_number_events);
  
  //Now we know the sizes, resize the event num,loc vectors
  fh_gid_event_num.resize(fh_number_events);
  fh_gid_byte_loc.resize(fh_number_events);
  for(unsigned int ee=0; ee<fh_number_events;ee++){
    input.readval(fh_gid_event_num[ee]);
    input.readval(fh_gid_byte_loc[ee]);
  }
}

//Assuming you've got event locations, read the sequences' livetime
void LEvtFile::GetLivetime(){
  input.readval(lt_num_sequences);
  lh_timestamp_latch.resize(lt_num_sequences);
  lh_timestamp_end.resize(lt_num_sequences);
  
  for(unsigned short ss=0; ss<lt_num_sequences;ss++){
    input.readval(lh_timestamp_latch[ss]);
    input.readval(lh_timestamp_end[ss]);
  }
}


//Assuming you've got livetime, read in the events to fill the vector<LEvent*> evt_events
void LEvtFile::GetEvents(){
  LEvtEvent *ev;
  size_t file_end_position = 0, last;
  input.seekg(0, std::ifstream::end);
  last = input.tellg();
  int status = 0;
  for(unsigned int ee=0;ee<fh_number_events;ee++){
    input.seekg(fh_gid_byte_loc[ee]);
    if (ee+1<fh_number_events) file_end_position = fh_gid_byte_loc[ee+1];
    else file_end_position = last;
    
    ev = new LEvtEvent();
    status = ev->Read(input, byteswap, settings_pt, file_end_position);
    if (status == -1) {	
    	delete ev;
    	ev = new LEvtEvent(122);
    	ev->gid_event_num = fh_gid_event_num[ee];
    }
    evt_events.push_back(ev);
  }
}

void LEvtFile::CalibrateChannels(std::vector<double> mVns_per_phe){
  double amps_gain = 7.5;
  try {
    amps_gain = settings_pt.get<float>("daq_settings.global.preamp")
      * settings_pt.get<float>("daq_settings.global.postamp");
  } catch (std::exception const& ex) { amps_gain = 7.5; }
  
  std::vector<LEvtChannel*> chs;
  for (std::vector<LEvtEvent*>::iterator ev = evt_events.begin();
    ev != evt_events.end(); ++ev) {
    //auto chs = (*ev)->channels;
    for(size_t i=0;(i<chs.size())&&(i<mVns_per_phe.size());++i){
      chs[i]->CalibratePods(mVns_per_phe.at(i),amps_gain);
    }
    // added to converte the water pmts into mV 
    for (size_t i=128; i<136 && i<chs.size(); i++) {
      // counter the 10ns/sample, 25x water pmts amplification
      chs[i]->CalibratePods(10.0, 25.0); 
    }
  }
}
