#include <iostream>
#include "Lux_EVT_Event.h"
#include "Lux_EVT_Channel.h"
#include "Lux_EVT_Pulse.h"
#include "TClonesArray.h"
using namespace std;

ClassImp(Lux_EVT_Pulse);

Lux_EVT_Pulse::Lux_EVT_Pulse(){
  pulseNumber = 0;
  pulseStarts = 0;
  pulseLength = 0;
  pulseBaseline = 0;
  channelNumber = 0;
  pulseData = new TH1I;
}

Lux_EVT_Pulse::Lux_EVT_Pulse(ULong_t pn, Long_t ps, ULong_t pl, ULong_t pb, ULong_t cn, TH1I* pd){
  pulseNumber = pn;
  pulseStarts = ps;
  pulseLength = pl;
  pulseBaseline = pb;
  channelNumber = cn;
  pulseData = pd;
}

Lux_EVT_Pulse::Lux_EVT_Pulse(Lux_EVT_Pulse *p){
  pulseNumber = p->Get_pulseNumber();
  pulseStarts = p->Get_pulseStarts();
  pulseLength = p->Get_pulseLength();
  pulseBaseline = p->Get_pulseBaseline();
  channelNumber = p->Get_channelNumber();
  pulseData = p->Get_pulseData(); 
}

Lux_EVT_Pulse::~Lux_EVT_Pulse(){}

void Lux_EVT_Pulse::Fill(ULong_t pn, Long_t ps, ULong_t pl, ULong_t pb, ULong_t cn, TH1I *pd ){
  pulseNumber = pn;
  pulseStarts = ps;
  pulseLength = pl;
  pulseBaseline = pb;
  channelNumber = cn;
  pulseData = pd;
}

void Lux_EVT_Pulse::Clear(Option_t* /*option*/){
  pulseNumber= 0; 
  pulseStarts = 0; 
  pulseLength = 0; 
  pulseBaseline = 0; 
  channelNumber = 0;
  /*pulseData->Clear();*/
}
