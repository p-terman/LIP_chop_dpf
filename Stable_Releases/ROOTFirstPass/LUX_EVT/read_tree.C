#include <stdlib.h>

#include "Riostream.h"
#include "TROOT.h"
#include "TFile.h"
#include "TNetFile.h"
#include "TRandom.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "TStopwatch.h"
#include "TTreeCacheUnzip.h"
#include "TRandom3.h"
#include "TH1I.h"
#include "TCanvas.h"

#include "Lux_EVT_Header.h"
#include "Lux_EVT_Channel.h"
#include "Lux_EVT_Event.h"
#include "Lux_EVT_Pulse.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
using namespace std;

int main(int argc, char* argv[]){    // arguments are the input root file name
  char *inFileName;
  inFileName = argv[1];
  TFile f(inFileName);				// open the input root file
  TTree *T = (TTree*)f.Get("Lux_EVT_Tree");	// get the data tree from the root file
  TCanvas *c1 = new TCanvas("c1","c1");

  Lux_EVT_Event *event = new Lux_EVT_Event();		// create the necessayr Lux_EVT class objects
  Lux_EVT_Header *header = new Lux_EVT_Header();
  Lux_EVT_Channel *channel = new Lux_EVT_Channel();
  Lux_EVT_Pulse *pulse = new Lux_EVT_Pulse();
  char filename[100]; 

  T->SetBranchAddress("event", &event);			// set the branch addresses of the tree so the data can be located by the Lux_EVT classes
  T->SetBranchAddress("header", &header);
  Int_t nentries = (Int_t)T->GetEntries();		// get the number of entries in the tree
  //  HEADER
  //nentries = 1;
  T->GetEntry(0);
  cerr<<"Header:        "<<endl;
  cerr<<"          dateTime: "<<header->Get_dateTime()<<endl;
  cerr<<"          location: "<<header->Get_location()<<endl;
  cerr<<"        firstEvent: "<<header->Get_firstEvent()<<endl;
  cerr<<"         numEvents: "<<header->Get_numEvents()<<endl;
  cerr<<"            numSeq: "<<header->Get_numSeq()<<endl;
  cerr<<"    timeStampLatch: ";
  if (header->Get_numSeq()>0) {
    for (Int_t seq = 0;seq<header->Get_numSeq();seq++){
      cerr<<header->Get_timeStampLatch(seq)<<", ";
    }
  }
  cerr<<endl;
  cerr<<"    timeStampEnd: ";
  if (header->Get_numSeq()>0){
    for (Int_t seq = 0;seq<header->Get_numSeq();seq++){
      cerr<<header->Get_timeStampEnd(seq)<<", ";
    }
  }
  cerr<<endl;
  for (Int_t ev=0;ev<nentries;ev++) {			// scan through the tree entries
    T->GetEntry(ev);					// get the current entries data
    //  if (ev) continue; //dump first event only
    cerr<<" Event:        "<<ev<<endl;				// display the data
    cerr<<"    eDateTime: "<<event->Get_eDateTime()<<endl;
    cerr<<"    eLocation: "<<event->Get_eLocation()<<endl;
    cerr<<"     eventNum: "<<event->Get_eventNum()<<endl;
    cerr<<"  numChannels: "<<event->Get_numChannels()<<endl;
    cerr<<"    numPulses: "<<event->Get_numPulses()<<endl;
    cerr<<"    eventSize: "<<event->Get_eventSize()<<endl;
    cerr<<"   recordForm: "<<event->Get_recordForm()<<endl;
    cerr<<"    timeStamp: "<<event->Get_timeStamp()<<endl<<endl;
    for (ULong_t ch = 0;ch < event->Get_numChannels(); ch++){		// scan through this events channels
      channel->Clear();
      channel = event->Get_channel(ch);					// fill the Lux_EVT_Channel instance channel with this channels data
      cerr<<" Channel:         "<<ch<<endl;
      cerr<<"   channelNumber: "<<channel->Get_channelNumber()<<endl;	// display the data
      cerr<<"      voltageRes: "<<channel->Get_voltageRes()<<endl;
      cerr<<"      voltageOff: "<<channel->Get_voltageOff()<<endl;
      cerr<<"         timeRes: "<<channel->Get_timeRes()<<endl;
      cerr<<"      preTrigger: "<<channel->Get_preTrigger()<<endl;
      cerr<<"       eventSize: "<<channel->Get_eventSize()<<endl;	
      cerr<<"     pulseDetect: "<<channel->Get_pulseDetect()<<endl;
      cerr<<"        pulseEnd: "<<channel->Get_pulseEnd()<<endl;
      cerr<<"       numPulses: "<<channel->Get_numPulses()<<endl<<endl;
    }
    for (ULong_t p = 0;p < event->Get_numPulses(); p++){		// scan through this events pulses
      pulse->Clear();
      pulse = event->Get_pulse(p);					// fill the Lux_EVT_Pulse instance pulse with this pulses data
      cerr<<" Pulse:           "<<p<<" of "<<event->Get_numPulses()<<endl;  	//display the data
      cerr<<"     pulseNumber: "<<pulse->Get_pulseNumber()<<endl;
      cerr<<"     pulseStarts: "<<pulse->Get_pulseStarts()<<endl;
      cerr<<"     pulseLength: "<<pulse->Get_pulseLength()<<endl;
      cerr<<"   pulseBaseline: "<<pulse->Get_pulseBaseline()<<endl;
      cerr<<"   channelNumber: "<<pulse->Get_channelNumber()<<endl<<endl;
      pulse->Get_pulseData()->Draw();						// draw the pulse trace for this pulse
      sprintf(filename,"plots/pulse_hist_%d_%d_%d.gif",(int)event->Get_eventNum(),(int)pulse->Get_channelNumber(),(int)pulse->Get_pulseNumber());
      c1->Print(filename);			// save the pulse trace for this pulse.
    }
  }
}

