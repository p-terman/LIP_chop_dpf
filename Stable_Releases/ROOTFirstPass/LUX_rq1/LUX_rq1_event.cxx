#include <iostream>
#include "LUX_rq1_pulse.h"
#include "LUX_rq1_header.h"
#include "LUX_rq1_event.h"
#include "LUX_rq1_channel.h"
using namespace std;

ClassImp(LUX_rq1_event)

LUX_rq1_event::LUX_rq1_event(){
  num_pulses = 0;
  pulse = new TClonesArray("LUX_rq1_pulse");
  pulse->BypassStreamer();
  num_channels = 0;
  channel = new TClonesArray("LUX_rq1_channel");
  channel->BypassStreamer();
}

void LUX_rq1_event::Clear(const Option_t* option){
  num_pulses=0;
  num_channels=0;
}

LUX_rq1_event::LUX_rq1_event(UInt_t a, UInt_t b){
  num_pulses=a;
  pulse = new TClonesArray("LUX_rq1_pulse");
  pulse->BypassStreamer();
  num_channels=b;
  channel = new TClonesArray("LUX_rq1_channel;");
  channel->BypassStreamer();
}

LUX_rq1_event::LUX_rq1_event(LUX_rq1_event *lrq1e){
  num_pulses = lrq1e->Get_num_pulses();
  num_channels = lrq1e->Get_num_channels();
}

LUX_rq1_event::~LUX_rq1_event(){delete pulse;delete channel;}

void LUX_rq1_event::AddPulse(LUX_rq1_pulse *lrq1p, Int_t index){
  TClonesArray &pls = *pulse;
  pls[index] = new(pls[index]) LUX_rq1_pulse(lrq1p);
}

LUX_rq1_pulse* LUX_rq1_event::Get_pulse(UInt_t num){
  if (num<num_pulses) return (LUX_rq1_pulse*)pulse->At(num);
  else {
    cerr << "Pulse index out of range, returning 0"<<endl;
    return 0;
  }
}

void LUX_rq1_event::AddChannel(LUX_rq1_channel *lrq1c, Int_t index){
  TClonesArray &chan = *channel;
  chan[index] = new(chan[index]) LUX_rq1_channel(lrq1c);
}

LUX_rq1_channel* LUX_rq1_event::Get_channel(UInt_t n_pulse, UInt_t n_channel){
  UInt_t chan_number = n_pulse*num_channels+n_channel;
  UInt_t total_chans = num_channels*num_pulses;
  if (chan_number<total_chans) return (LUX_rq1_channel*)channel->At(chan_number);
  else {
    cerr << "Channel index out of range, returning 0"<<endl;
    return 0;
  }
}
