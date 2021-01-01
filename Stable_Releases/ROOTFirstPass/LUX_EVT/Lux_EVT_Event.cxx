#include <iostream>
#include "Lux_EVT_Event.h"
#include "Lux_EVT_Channel.h"
#include "Lux_EVT_Pulse.h"
#include "TClonesArray.h"

using namespace std;

ClassImp(Lux_EVT_Event)

Lux_EVT_Event::Lux_EVT_Event(){
  eDateTime = 0; 
  eLocation = 0;
  eventNum = 0;
  numChannels = 0;
  numPulses = 0;
  eventSize = 0;
  recordForm = 0;
  recordSize = 0;
  timeStamp = 0;
  channel = new TClonesArray("Lux_EVT_Channel"); 
  channel->BypassStreamer();
  pulse = new TClonesArray("Lux_EVT_Pulse"); 
  pulse->BypassStreamer();
}

Lux_EVT_Event::Lux_EVT_Event(ULong_t edt, ULong_t el, ULong_t en, ULong_t nc, ULong_t np, ULong_t es, ULong_t rf, ULong_t rs, ULong64_t ts){
  eDateTime = edt;
  eLocation = el;
  eventNum = en;
  numChannels = nc;
  numPulses = np;
  eventSize = es;
  recordForm = rf;
  recordSize = rs;
  timeStamp = ts;
  channel = new TClonesArray("Lux_EVT_Channel");
  channel->BypassStreamer();
  pulse = new TClonesArray("Lux_EVT_Pulse"); 
  pulse->BypassStreamer();
}

Lux_EVT_Event::Lux_EVT_Event(Lux_EVT_Event *e){
  eDateTime = e->Get_eDateTime();
  eLocation = e->Get_eLocation();
  eventNum = e->Get_eventNum();
  numChannels = e->Get_numChannels();
  numPulses = e->Get_numPulses();
  eventSize = e->Get_eventSize();
  recordForm = e->Get_recordForm();
  recordSize = e->Get_recordSize();
  timeStamp = e->Get_timeStamp();
  channel = e->Get_channel();
}

Lux_EVT_Event::~Lux_EVT_Event(){delete channel; delete pulse;}

void Lux_EVT_Event::Fill(ULong_t edt, ULong_t el, ULong_t en, ULong_t nc, ULong_t np, ULong_t es, ULong_t rf, ULong_t rs, ULong64_t ts){
  eDateTime = edt; 
  eLocation = el;
  eventNum = en;
  numChannels = nc;
  numPulses = np;
  eventSize = es;
  recordForm = rf;
  recordSize = rs;
  timeStamp = ts;
}

void Lux_EVT_Event::AddChannel(Lux_EVT_Channel *c, Int_t index){
  TClonesArray &chan = *channel;
  chan[index] = new(chan[index]) Lux_EVT_Channel(c);
}

void Lux_EVT_Event::Clear(Option_t* /*option*/){
  eDateTime = 0; 
  eLocation = 0;
  eventNum = 0;
  numChannels = 0;
  numPulses = 0;
  eventSize = 0;
  recordForm = 0;
  recordSize = 0;
  timeStamp = 0;
  channel->Clear();
  pulse->Clear();
}

void Lux_EVT_Event::AddPulse(Lux_EVT_Pulse *p, Int_t index){
  TClonesArray &pls = *pulse;
  pls[index] = new(pls[index]) Lux_EVT_Pulse(p);
}

Lux_EVT_Pulse* Lux_EVT_Event::Get_pulse(ULong_t index){
  if (index<numPulses) return (Lux_EVT_Pulse*)pulse->At(index);
  else {
    cerr<<"Pulse index out of range: returning 0"<<endl;
    return 0;
  }
}

Lux_EVT_Channel* Lux_EVT_Event::Get_channel(ULong_t index){
  if (index<numChannels) return (Lux_EVT_Channel*)channel->At(index);
  else {
    cerr<<"Channel index out of range: returning 0"<<endl;
    return 0;
  }
}


