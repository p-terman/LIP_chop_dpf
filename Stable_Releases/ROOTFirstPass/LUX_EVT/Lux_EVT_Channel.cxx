#include <iostream>

#include "Lux_EVT_Channel.h"
#include "Lux_EVT_Event.h"
#include "Lux_EVT_Pulse.h"
#include "TClonesArray.h"
using namespace std;

ClassImp(Lux_EVT_Channel)

Lux_EVT_Channel::Lux_EVT_Channel(){
  binDataType = '0';
  channelNumber = 0;
  voltageRes = 0;
  voltageOff = 0;
  timeRes = 0;
  preTrigger = 0;
  eventSize = 0;
  pulseDetect = 0;
  pulseEnd = 0;
  numPulses = 0;
}
Lux_EVT_Channel::Lux_EVT_Channel(Char_t bdt, ULong_t cn, Double_t vr, Double_t vo, Double_t tr, Long_t pt, ULong_t es, ULong_t pd, ULong_t pe, ULong_t np){
  binDataType = bdt;
  channelNumber = cn;
  voltageRes = vr;
  voltageOff = vo;
  timeRes = tr;
  preTrigger = pt;
  eventSize = es;
  pulseDetect = pd;
  pulseEnd = pe;
  numPulses = np;
}

Lux_EVT_Channel::Lux_EVT_Channel(Lux_EVT_Channel *c){
  binDataType = c->Get_binDataType();
  channelNumber = c->Get_channelNumber(); 
  voltageRes = c->Get_voltageRes();
  voltageOff = c->Get_voltageOff();
  timeRes = c->Get_timeRes();
  preTrigger = c->Get_preTrigger();
  eventSize = c->Get_eventSize();
  pulseDetect = c->Get_pulseDetect();
  pulseEnd = c->Get_pulseEnd();
  numPulses = c->Get_numPulses();
}

void Lux_EVT_Channel::Fill(Char_t bdt, ULong_t cn, Double_t vr, Double_t vo, Double_t tr, Long_t pt, ULong_t es, ULong_t pd, ULong_t pe, ULong_t np){
  binDataType = bdt;
  channelNumber = cn;
  voltageRes = vr;
  voltageOff = vo;
  timeRes = tr;
  preTrigger = pt;
  eventSize = es;
  pulseDetect = pd;
  pulseEnd = pe;
  numPulses = np;
}

void Lux_EVT_Channel::Clear(Option_t* /*option*/){

  binDataType = '0';
  channelNumber = 0; 
  voltageRes = 0; 
  voltageOff = 0; 
  timeRes = 0; 
  preTrigger = 0; 
  eventSize = 0; 
  pulseDetect = 0; 
  pulseEnd = 0; 
  numPulses = 0; 
}


