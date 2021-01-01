#include <stdlib.h>
#include <algorithm>

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

#include "Lux_EVT_Event.h"
#include "Lux_EVT_Channel.h"
#include "Lux_EVT_Pulse.h"
#include "Lux_EVT_Header.h"

#include "XmlParser.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <vector>
using namespace std;

#define EndianByteSwap(x) ByteSwap((unsigned char *) &x,sizeof(x))

void ByteSwap(unsigned char * b, int n){ // taken from http://www.codeproject.com/KB/cpp/endianness.aspx
   register int i = 0;
   register int j = n-1;
   while (i<j){
      std::swap(b[i], b[j]);
      i++, j--;
   }
}

// void EndianByteSwapUlong(unsigned int *var){  //endian-ness fixing routine
// 	*var = (*var>>24) |
// 		((*var<<8) & 0x00FF0000) |
// 		((*var>>8) & 0x0000FF00) |
// 		(*var<<24);
// }

int main(int argc, char* argv[]){  //  inputs are the evt file name and verbosity 
	TStopwatch timer;
	timer.Start();
	Int_t verbosity = 0;	   // verbosity defaults to 0, which means no printed output.  A value greater than 0 results in printed output
	char *inFileName;
	Int_t version = 2;    // 2 for newer files

	if (argc>2) verbosity = atoi(argv[2]); 		// handle the runtime arguments
	if (argc>1) inFileName = argv[1];	
	else{
		cerr<<"ERROR: NO FILENAME ENTERED"<<endl;
		return 0;
	}
	if (verbosity>0) cerr<<"inFileName is "<<inFileName<<endl;	
	Int_t nptot=0;					// counts the number of pulses in an event
    
    //If your int size or long long is different you may need to do a find/replace
    unsigned int uintTest;
    unsigned long long int ullintTest;
    if( (sizeof(uintTest) != 4) || (sizeof(ullintTest) != 8) ) {
        cout << "* Your compiler's UInt and long long is not compatible *" << endl;
        cout << "* If you know uints are 4 and 8 bytes, you can do a    *" << endl;
        cout << "* find-and-replace search in ./Make_EVT_Tree.cxx, or   *" << endl;
        cout << "* email Nick at niwalsh@ucdavis.edu                 *" << endl;
        exit(0);
    }

	Lux_EVT_Event *event = new Lux_EVT_Event();	// intialize a Lux_EVT_Event object
	Lux_EVT_Channel *cnl = new Lux_EVT_Channel();	// intialize a Lux_EVT_Channel object
	Lux_EVT_Pulse *puls = new Lux_EVT_Pulse();	// intialize a Lux_EVT_Pulse object
	Lux_EVT_Header *header = new Lux_EVT_Header();	// intialize a Lux_EVT_Header object

	Int_t split = 2;				// the split level must be at least two so the tree will create all branches and subbranches
	Int_t bsize = 32000;				// buffer size for each branch

	char outFileName[200];
	strcpy(outFileName,inFileName);			// copy the input file name to the output file name
	strtok(outFileName,".");			// remove the file type identifier (".evt") from the output file name
	strcat(outFileName,".root");			// add a new file type identifier (".root") to the output file name
	if (verbosity>0) cerr<<"outFileName is "<<outFileName<<endl;
	TFile *f = new TFile(outFileName,"recreate");	// creation of the tree
	TTree *Lux_EVT_Tree = new TTree("Lux_EVT_Tree","containing Lux EVT data");  //create the tree
	Lux_EVT_Tree->Branch("header", "Lux_EVT_Header", &header, bsize,1);	// adds the header branch to the tree
	Lux_EVT_Tree->Branch("event", "Lux_EVT_Event", &event, bsize,split);	// adds the event branch (automatically contains channels and pulses) to the tree


	ifstream inFile;				
	inFile.open (inFileName, ios::binary);			// open the input file
	unsigned int endianTest;
	inFile.read( (char*)&endianTest, sizeof(endianTest) );  // read the endian-ness from the input file
	bool BYTESWAP = true;					
	if (endianTest == 0x01020304) BYTESWAP = false;		// determines whether the byte order of the input variables must be reversed

	if(verbosity==-1){
		if(BYTESWAP) cerr<<"Endianness must be changed"<<endl<<endl;
		else cerr<<"Endianness matched"<<endl<<endl;
	}

    //XML Header
    unsigned int xmlLength;
    inFile.read( (char*)&xmlLength, sizeof(xmlLength) );
    if(verbosity>0) cout << "XML Length: " << xmlLength << endl;
    char * xmlBuffer= new char [xmlLength];
    inFile.read( xmlBuffer, xmlLength );
    if(verbosity>1) cout << *xmlBuffer << endl;
    XmlParser *xmlParse = new XmlParser(xmlBuffer);
    double evt_bldr_version = xmlParse->GetEvtSettingsEventBuilderVersion() ;
    if(verbosity>0) 
        cout << "XML:evt_settings:event_builder_version " 
             <<  evt_bldr_version << endl;
    evt_bldr_version = xmlParse->GetGlobalEventBuilderVersion() ;
    if(verbosity>0)
        cout <<"XML:global:event_builder_version " <<  evt_bldr_version << endl;
    double daq_version = xmlParse->GetGlobalDAQVersion() ;
    if(verbosity>0) cout << "XML:global:daq_version " <<  daq_version << endl;
    delete [] xmlBuffer;
    delete xmlParse;

    //Re-do endian test--format current as of April 2010
    //ifFile.read( (char*)&endiantTest, sizeof(endianTest));
    inFile.ignore( sizeof(endianTest) );

//	file header info
								// These 4 lines of code will be repeated throughout the program for different
								// variables, but will only be commented here.
	unsigned int dateTime;				// defines the variable dateTime
	inFile.read( (char*)&dateTime, sizeof(dateTime) );      // read the variable dateTime from the binary file
	if(BYTESWAP) EndianByteSwap(dateTime);		// changes the byte order of the variable dateTime if necessary 
	if (verbosity>0) cout << "Date and Time = " << dateTime << endl;  // displays the value of dateTime if requested

	unsigned int location;
	inFile.read( (char*)&location, sizeof(location) );
	if(BYTESWAP) EndianByteSwap(location);
	if (verbosity>0) cout << "Location = " << location << endl;

	unsigned int firstEvent;
	inFile.read( (char*)&firstEvent, sizeof(firstEvent) );	
	if(BYTESWAP) EndianByteSwap(firstEvent);
	if (verbosity>0) cout << "First Event = " << firstEvent << endl;
	
	unsigned int numEvents;
	inFile.read( (char*)&numEvents, sizeof(numEvents) );
	if(BYTESWAP) EndianByteSwap(numEvents);
	if (verbosity>0) cout << "Number of Events = " << numEvents << endl;

	unsigned short int numSeq = 0;
	if(version == 2){
	  //	live time header
	  inFile.read( (char*)&numSeq,sizeof(numSeq) );
	  if(BYTESWAP) EndianByteSwap(numSeq);
	  if(verbosity>0) cout << "Number of Sequences "<<numSeq << endl;;
	
	  unsigned long long int timeStampLatch[numSeq];
	  unsigned long long int timeStampEnd[numSeq];
	  for(int s=0; s<numSeq; s++){	
	    inFile.read( (char*)&timeStampLatch[s], sizeof(timeStampLatch[s]) );
	    if(BYTESWAP) EndianByteSwap(timeStampLatch[s]);
	    header->Set_timeStampLatch((ULong_t)timeStampLatch[s],(Int_t)s);
	    if(verbosity>0) cout << timeStampLatch[s] << "\t";
	 
	    inFile.read( (char*)&timeStampEnd[s], sizeof(timeStampEnd[s]) );
	    if(BYTESWAP) EndianByteSwap(timeStampEnd[s]);
	    header->Set_timeStampEnd((ULong_t)timeStampEnd[s],(Int_t)s);
	    if(verbosity>0) cout << timeStampEnd[s] << "\t";
	 	  }
	  if(verbosity>0) cout << endl;
	}
	
	header->Fill(ULong_t(dateTime),ULong_t(location),ULong_t(firstEvent),ULong_t(numEvents),ULong_t(numSeq));  // fill the Lux_EVT_Header instance header with the appropriatge information
	
	//	event GID

	for(unsigned int e=0; e<numEvents; e++){			// scan through the events in the binary file
	  event->Clear();						// make sure the Lux_EVT_Event instance event is empty
	 
	  unsigned int eDateTime;
	  inFile.read( (char*)&eDateTime, sizeof(eDateTime) );          // read event level information from the binary file
	  if(BYTESWAP) EndianByteSwap(eDateTime);
	  if (verbosity>0) cout << "Event Date and Time = " << eDateTime << endl;
	  
	  unsigned int eLocation;
	  inFile.read( (char*)&eLocation, sizeof(eLocation) );
	  if(BYTESWAP) EndianByteSwap(eLocation);
	  if (verbosity>0) cout << "event Location = " << eLocation << endl;
	  
	  unsigned int eventNum;
	  inFile.read( (char*)&eventNum, sizeof(eventNum) );
	  if(BYTESWAP) EndianByteSwap(eventNum);
	  if (verbosity>0) cout << "Event number = " << eventNum << " of "<<numEvents<<endl;
	  
	  unsigned int numChannels;
	  inFile.read( (char*)&numChannels, sizeof(numChannels) );
	  if(BYTESWAP) EndianByteSwap(numChannels);
	  if (verbosity>0) cout << "Number of channels = " << numChannels << endl;
	  
	  //	record info
	  
	  unsigned int eventSize;	//not in documentation (event_builder_v3.c line 516)
	  inFile.read( (char*)&eventSize, sizeof(eventSize) );
	  if(BYTESWAP) EndianByteSwap(eventSize);
	  if (verbosity>0) cout << "Event size = " << eventSize << endl;
	  
	  unsigned int recordForm;
	  inFile.read( (char*)&recordForm, sizeof(recordForm) );
	  if(BYTESWAP) EndianByteSwap(recordForm);
	  if (verbosity>0) cout << "Record Format = " << recordForm << endl;
	  
	  unsigned int recordSize;
	  inFile.read( (char*)&recordSize, sizeof(recordSize) );
	  if(BYTESWAP) EndianByteSwap(recordSize);
	  if (verbosity>0) cout << "Record Size = " << recordSize << endl;
	  
	  //	timestamp
	  
	  unsigned long long int timeStamp;
	  inFile.read( (char*)&timeStamp, sizeof(timeStamp) );
	  if(BYTESWAP) EndianByteSwap(timeStamp);
	  if (verbosity>0) cout << "timeStamp = " << timeStamp << endl;
	  
	  //	event channel header
	  
	  for(unsigned int c=0; c<numChannels; c++){		// scan through the channels in this event
	    cnl->Clear();						// make sure the Lux_EVT_Channel instance cnl is empty
									// read the channel information from the binary file
	    char binDataType;	//should be 14, they are not using single quotes (event_builder_v3.c line 132)
	    inFile.read( (char*)&binDataType, sizeof(binDataType));
	    if(BYTESWAP) EndianByteSwap(binDataType);
	    if (verbosity>0) cout << "binary data type = " << binDataType << endl;
	    
	    double voltageRes;
	    inFile.read( (char*)&voltageRes, sizeof(voltageRes));
	    if(BYTESWAP) EndianByteSwap(voltageRes);
	    if (verbosity>0) cout << "voltage resolution = " << voltageRes << endl;
	    
	    double voltageOff;
	    inFile.read( (char*)&voltageOff, sizeof(voltageOff));
	    if(BYTESWAP) EndianByteSwap(voltageOff);
	    if (verbosity>0) cout << "voltage offset = " << voltageOff << endl;
	    
	    double timeRes;
	    inFile.read( (char*)&timeRes, sizeof(timeRes));
	    if(BYTESWAP) EndianByteSwap(timeRes);
	    if (verbosity>0) cout << "time resolution = " << timeRes << endl;
	    
	    signed int preTrigger;
	    inFile.read( (char*)&preTrigger, sizeof(preTrigger) );
	    if(BYTESWAP) EndianByteSwap(preTrigger);
	    if (verbosity>0) cout << "pretrigger = " << preTrigger << endl;
	    
	    unsigned int eventSize;
	    inFile.read( (char*)&eventSize, sizeof(eventSize) );
	    if(BYTESWAP) EndianByteSwap(eventSize);
	    if (verbosity>0) cout << "Event size = " << eventSize << endl;
	    
	    unsigned int pulseDetect;
	    inFile.read( (char*)&pulseDetect, sizeof(pulseDetect) );
	    if(BYTESWAP) EndianByteSwap(pulseDetect);
	    if (verbosity>0) cout << "pulse detect pretrigger = " << pulseDetect << endl;
	    
	    unsigned int pulseEnd;
	    inFile.read( (char*)&pulseEnd, sizeof(pulseEnd) );
	    if(BYTESWAP) EndianByteSwap(pulseEnd);
	    if (verbosity>0) cout << "pulse end posttrigger = " << pulseEnd << endl;
	    
	    unsigned int numPulses;
	    inFile.read( (char*)&numPulses, sizeof(numPulses) );
	    if(BYTESWAP) EndianByteSwap(numPulses);
	    if (verbosity>0) cout << "number of pulses = " << numPulses << endl;
	    
	    cnl->Fill(Char_t(binDataType), ULong_t(c), Double_t(voltageRes), Double_t(voltageOff), Double_t(timeRes), Long_t(preTrigger), ULong_t(eventSize), ULong_t(pulseDetect), ULong_t(pulseEnd), ULong_t(numPulses));     // fill the Lux_EVT_Channel instance cnl with the approproate information
	 
	    const unsigned int numPulsesconst = numPulses; 
	    signed int pulseStarts[numPulsesconst];			// scan through the pulses in the current channel
	    for(unsigned int np=0; np<numPulses; np++){
	      inFile.read( (char*)&pulseStarts[np],sizeof(pulseStarts[np]) );
	      if (verbosity>0) cout << "pulse start for pulse " << np << " = " << pulseStarts[np] <<endl;
	      
	    }
	    unsigned int pulseLengths[numPulsesconst];
	    for(unsigned int np=0; np<numPulses; np++){
	      inFile.read( (char*)&pulseLengths[np],sizeof(pulseLengths[np]) );
	      if(BYTESWAP) EndianByteSwap(pulseLengths[np]);
	      if (verbosity>0) cout << "pulse length for pulse " << np << " = " << pulseLengths[np] <<endl;
	      
	    }
	    unsigned int pulseBaselines[numPulsesconst];
	    for(unsigned int np=0; np<numPulses; np++){
	      inFile.read( (char*)&pulseBaselines[np],sizeof(pulseBaselines[np]) );
	      if(BYTESWAP) EndianByteSwap(pulseBaselines[np]);
	      if (verbosity>0) cout << "pulse baseline for pulse " << np << " = " << pulseBaselines[np] <<endl;
	    }					
												
	    for(unsigned int np=0; np<numPulses; np++){		// pulse data is stored in TH1I root histograms
	      puls->Clear();						// make sure the Lux_EVT_Pulse instance puls is empty
	      TH1I *pulse_data_hist = new TH1I("pulse_data_hist","pulse_data_hist",(Double_t)pulseLengths[np],(Double_t)pulseStarts[np],(Double_t)pulseStarts[np]+(Double_t)pulseLengths[np]);				// create the histogram for the current pulse
	      for(int pl=0; pl<(int)pulseLengths[np]; pl++){			// scan trough the pulse data and fill the histogram
		unsigned short int pulseData;
		inFile.read( (char*)&pulseData,sizeof(pulseData));
	        if(BYTESWAP) EndianByteSwap(pulseData);
		pulse_data_hist->Fill(pl+pulseStarts[np],pulseData-(unsigned short int)pulseBaselines[np]); 
	      }
	      puls->Fill(ULong_t(np), Long_t(pulseStarts[np]), ULong_t(pulseLengths[np]), ULong_t(pulseBaselines[np]), ULong_t(c) ,pulse_data_hist);  // fill the Lux_EVT_Pulse instance puls with the appropriate data
	      event->AddPulse(puls,nptot);	// add the Lux_EVT_Pulse puls to the Lux_EVT_Event event
	      nptot++;				// counts the number of pulses in the event
	      pulse_data_hist->Clear();		// clear the histogram
	    }
	    event->AddChannel(cnl,c);		// add the Lux_EVT_Channel cnl to the Lux_EVT_Event event
	    numPulses = 0;
	  }
	  event->Fill(ULong_t(eDateTime), ULong_t(eLocation), ULong_t(eventNum), ULong_t(numChannels), ULong_t(nptot), ULong_t(eventSize), ULong_t(recordForm), ULong_t(recordSize), ULong64_t(timeStamp));   //fill the Lux_EVT_Event event with the appropriate data
	  nptot = 0;
	  Lux_EVT_Tree->Fill();			    // write the data to the tree
	}
	f->Write();						// write the tree to the root file
	f->Close();						// close the root file
	timer.Stop();
	Double_t rtime = timer.RealTime();
	Double_t ctime = timer.CpuTime();
	if (verbosity>0) printf("RealTime=%f seconds, CpuTime=%f seconds\n",rtime,ctime);
	return 0;
}


