#include <iostream>
#include <cstring>
#include <exception>
#include <vector>
#include "LCvtFile.h"

#define VERBOSE 0

#define writeval(a) write((char*)&a, sizeof(a))
#define writestr(a,b) write((char*)a, b)
#define writevec(a,b) write((char*)&a[0], sizeof(a[0])*b)

using namespace std;
LCvtFile::LCvtFile(){}

LCvtFile::LCvtFile(std::string path, std::string cvt_file){
  byteswap = false;
  std::string fullpath=path+'/'+cvt_file;
  input.open(fullpath.data(), std::ios_base::binary);
  if(!input.is_open()){
    std::cout<<"Failed to open file: "<<fullpath.data()<<std::endl;
  }
  else{
    std::cout<<fullpath.data()<<std::endl;
    this->GetSettings();
#if VERBOSE>0
    cout<<"gotSettings"<<endl;
#endif
    this->GetFileHeader();
#if VERBOSE>0
    cout<<"gotFileHeader"<<endl;
#endif
    this->GetLivetime();
#if VERBOSE>0
    cout<<"gotLivetime"<<endl;
#endif
    this->GetEvents();
#if VERBOSE>0
    cout<<"gotEvents"<<endl;
#endif
    input.close();
  }
}

LCvtFile::LCvtFile(LEvtFile *evt) :   
  endianness(evt->endianness),
  byteswap(evt->byteswap),
  settings_string_len(evt->settings_string_len),
  settings_pt(evt->settings_pt),
  fh_endianness(evt->fh_endianness),
  fh_date_time(evt->fh_date_time),
  fh_location(evt->fh_location),
  fh_number_events(evt->fh_number_events),
  fh_gid_event_num(evt->fh_gid_event_num),
  fh_gid_byte_loc(evt->fh_gid_byte_loc),
  lt_num_sequences(evt->lt_num_sequences),
  lh_timestamp_latch(evt->lh_timestamp_latch),
  lh_timestamp_end(evt->lh_timestamp_end),
  input_filename(evt->input_filename)
{
  // copy settings_string
  settings_string = new char[settings_string_len+1];
  strcpy(settings_string, evt->settings_string);
  // copy events
  for (unsigned int ee=0; ee<fh_number_events; ee++) {
    cvt_events.push_back(new LCvtEvent(evt->evt_events[ee]));
  }
}


LCvtFile::~LCvtFile(){
  for (std::vector<LCvtEvent*>::iterator iter = cvt_events.begin();
    iter != cvt_events.end(); ++iter) {
    delete *iter;
  }
  delete [] settings_string;
}

//Assuming input is an open fstream, check endianness then read settings_string
void LCvtFile::GetSettings() {
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
void LCvtFile::GetFileHeader(){
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
void LCvtFile::GetLivetime(){
  input.readval(lt_num_sequences);
  lh_timestamp_latch.resize(lt_num_sequences);
  lh_timestamp_end.resize(lt_num_sequences);
  
  for(unsigned short ss=0; ss<lt_num_sequences;ss++){
    input.readval(lh_timestamp_latch[ss]);
    input.readval(lh_timestamp_end[ss]);
  }
}


//Assuming you've got livetime, read in the events to fill the vector<LEvent*> cvt_events
void LCvtFile::GetEvents(){
  for(unsigned int ee=0;ee<fh_number_events;ee++){
    input.seekg(fh_gid_byte_loc[ee]);
    cvt_events.push_back(new LCvtEvent(input, byteswap, settings_pt));
#if VERBOSE>0
    cout<<ee<<endl;
#endif
  }
}


void LCvtFile::WriteCvtFile(std::string path, std::string cvt_file) {
  this->UpdateFileLocations();
  
  string fullpath=path+'/'+cvt_file;
  ofstream output(fullpath.data(), std::ios_base::binary);
  // write the xml part 
  endianness = 0x01020304;
  output.writeval(endianness);
  output.writeval(settings_string_len);
  output.writestr(settings_string, settings_string_len);
  // write the file header
  output.writeval(fh_endianness);  
  output.writeval(fh_date_time);
  output.writeval(fh_location);
  output.writeval(fh_number_events);
  for(unsigned int ee=0; ee<fh_number_events;ee++){
    output.writeval(fh_gid_event_num[ee]);
    output.writeval(fh_gid_byte_loc[ee]);
  }
  // write the livetime 
  output.writeval(lt_num_sequences);
  for(unsigned short ss=0; ss<lt_num_sequences;ss++){
    output.writeval(lh_timestamp_latch[ss]);
    output.writeval(lh_timestamp_end[ss]);
  }
  // write the events 
  LCvtEvent *ev = NULL;
  LCvtChannel *ch = NULL;
  vector<float> zeros;
  for(unsigned int ee=0;ee<fh_number_events;ee++){
    ev = cvt_events[ee];
    
    output.writeval(ev->gid_date_time);
    output.writeval(ev->gid_location);
    output.writeval(ev->gid_event_num);
    output.writeval(ev->gid_num_chans);
    output.writeval(ev->gid_event_size);
    output.writeval(ev->record_format);
    output.writeval(ev->record_size);    
    output.writeval(ev->event_trigger_timestamp);
    // write the channels
    for(unsigned int cc=0;cc<ev->gid_num_chans; cc++){
      ch = ev->channels[cc];
      
      output.writeval(ch->ch_binary_datatype);
      output.writeval(ch->voltage_resolution);
      output.writeval(ch->voltage_offset);
      output.writeval(ch->time_resolution);
      output.writeval(ch->pretrigger);
      output.writeval(ch->event_size);
      output.writeval(ch->pulse_detect_pretrigger);
      output.writeval(ch->pulse_detect_posttrigger);
      output.writeval(ch->number_of_pods);      
      // read pods
      if (ch->number_of_pods==0) continue;
      output.writevec(ch->pod_starts, ch->number_of_pods);
      output.writevec(ch->pod_lengths, ch->number_of_pods);
      zeros.resize(ch->number_of_pods, 0);
      output.writevec(zeros, ch->number_of_pods);
      for(unsigned short pp=0;pp<ch->number_of_pods;pp++){
        output.writevec(ch->pod_data_phe[pp], ch->pod_lengths[pp]);
      }
    } // end of channel loop
  } // end of event loop
  
  output.close();
}

void LCvtFile::UpdateFileLocations() {
  size_t current = 0;
  current += sizeof(endianness);
  current += sizeof(unsigned int) + settings_string_len;
  current += sizeof(fh_endianness);
  current += sizeof(fh_date_time);
  current += sizeof(fh_location);
  current += sizeof(fh_number_events);
  current += sizeof(unsigned int)*fh_number_events;//sizeof(fh_gid_event_num[0])*fh_number_events;
  current += sizeof(unsigned int)*fh_number_events;//sizeof(fh_gid_byte_loc[0])*fh_number_events;
  current += sizeof(lt_num_sequences);
  current += sizeof(unsigned long long)*lt_num_sequences;//sizeof(lh_timestamp_latch[0])*lt_num_sequences;
  current += sizeof(unsigned long long)*lt_num_sequences;////sizeof(lh_timestamp_end[0])*lt_num_sequences;
  for(unsigned int ee=0;ee<fh_number_events;ee++){
    fh_gid_byte_loc[ee] = current;
    current += cvt_events[ee]->GetSize();
  }
}

void LCvtFile::Convert(vector<double> pmt_gains) {
  float preamp = 5;
  try {
    preamp = settings_pt.get<float>("daq_settings.global.preamp");
  } catch(exception const& ex) { preamp = 5; }
  float postamp = 1.5;
  try { 
    postamp = settings_pt.get<float>("daq_settings.global.postamp");
  } catch(exception const& ex) { postamp = 1.5; }
  float w_amp = 25.0;
  for (size_t ee=0; ee<cvt_events.size(); ee++) {
    for (size_t cc=0; cc<cvt_events[ee]->channels.size(); cc++) {
      if (cc < 122) {
        cvt_events[ee]->channels[cc]->Convert(pmt_gains[cc]*preamp*postamp);
      }else if (cc >= 128 and cc < 136) {
        if (cc < pmt_gains.size()) 
        	cvt_events[ee]->channels[cc]->Convert(w_amp*pmt_gains[cc]);
        else cvt_events[ee]->channels[cc]->Convert(w_amp*10);
      }else {
        cvt_events[ee]->channels[cc]->Convert(10);
      }
    }
  }  
}

