#include <iostream>
#include "Lux_EVT_Header.h"
#include "Lux_EVT_Event.h"
#include "Lux_EVT_Channel.h"
#include "Lux_EVT_Pulse.h"
using namespace std;

ClassImp(Lux_EVT_Header);


Lux_EVT_Header::Lux_EVT_Header(){
  dateTime = 0;
  location = 0;
  firstEvent = 0;
  numEvents = 0;
  numSeq = 0;
  timeStampLatch = new TArrayL();
  timeStampEnd = new TArrayL();
}

Lux_EVT_Header::Lux_EVT_Header(ULong_t dt, ULong_t l, ULong_t fe, ULong_t ne, ULong_t ns){
  dateTime = dt;
  location = l;
  firstEvent = fe;
  numEvents = ne;
  numSeq = ns;
  timeStampLatch = new TArrayL();
  timeStampEnd = new TArrayL();
}

Lux_EVT_Header::Lux_EVT_Header(Lux_EVT_Header *h){
  dateTime = h->Get_dateTime();
  location = h->Get_location();
  firstEvent = h->Get_firstEvent();
  numEvents = h->Get_numEvents();
  numSeq = h->Get_numSeq();
  timeStampLatch = h->Get_timeStampLatch();
  timeStampEnd = h->Get_timeStampEnd();
}

Lux_EVT_Header::~Lux_EVT_Header(){}

void Lux_EVT_Header::Fill(ULong_t dt, ULong_t l, ULong_t fe, ULong_t ne, ULong_t ns){
  dateTime = dt;
  location = l;
  firstEvent = fe;
  numEvents = ne;
  numSeq = ns;
}

void Lux_EVT_Header::Clear(Option_t* /*option*/){
	dateTime = 0; 
	location = 0; 
	firstEvent = 0; 
	numEvents = 0;
	numSeq = 0;
} 
