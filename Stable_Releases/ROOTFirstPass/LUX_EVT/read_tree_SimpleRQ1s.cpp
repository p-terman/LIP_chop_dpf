 #include <stdlib.h>
#include <cstdlib>

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

#include <functional>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <ctime>
using namespace std;

Long_t pulseStarts;
Long_t pulseLength;

void reverse(Double_t orig[], Int_t b)
{
    Int_t a=0;
    double swap;
    for(a;a<--b;a++) //increment a and decrement b until they meet eachother
    {
        swap=orig[a];       //put what's in a into swap space
        orig[a]=orig[b];    //put what's in b into a
        orig[b]=swap;       //put what's in the swap (a) into b
    }
    return;    
}

Double_t n;
Double_t w;
Int_t nd;

double peakfind (Double_t arg[],Int_t length, Double_t arg2[]) { //find the maximum value of the pulse and the bin at which it occurs
   n = arg[0];
   w = arg2[0];
   for (Int_t d=3; d<length; d++) { 
     if (arg[d]>n) { n=arg[d]; w=arg2[d]; nd=d; }
   }
   return n;
}

Double_t sigma;
Double_t mean;
void stdev (Double_t arg[], Int_t length) { //compute the standard deviation of the end of the pulse (used in the pulsestart/end)
   Double_t sum = 0;
   Int_t s = 14;
   for (Int_t m=(length-s); m<length; m++) { //take the end of the pulse and find the mean
     sum += arg[m];
   }
   mean = (sum/s);
   Double_t stdevalue[s]; //reassign values from pdata to new array
   for (Int_t m=0; m<(s+1); m++) {
     stdevalue[m] = arg [m+(length-s)];
   }
   for (Int_t m=0; m<s; m++) { //compute deviation 
        stdevalue[m] = (stdevalue[m]-mean);
   }
   for (Int_t m=0; m<s; m++) {  //square deviation
     stdevalue[m] *= stdevalue[m]; 
   }
   Double_t sum2 = 0;
   for (Int_t m=0; m<s; m++) { //add up squared deviations
        sum2 +=stdevalue[m];
   }
   sigma = (sqrt(sum2/(s-1))); //divide by one less than sample size and take the square root
   return;
}


double sigmaev;
double meanev;


double meanevent (Double_t arg[], Int_t length, /*Int_t numberpulses, Int_t currentpulse,*/ Double_t evvalue[], Int_t lev) {
  //find the mean and st. deviation of all non-zero entries in an event
double sumev;

   for (int m=0; m<lev; m++) {
     if (m==0) sumev=0;
     else sumev += evvalue[m];
   }
  meanev = (sumev/lev);
  cerr << "Mean for this event: " << meanev << endl;
  return meanev;
}

double sigmaevent (Double_t arg[], Int_t length, Double_t evvalue[], Int_t lev) {

  for (Int_t m=0; m<lev; m++) { //compute deviation 
        evvalue[m] = (evvalue[m]-meanev);
   }
   for (Int_t m=0; m<lev; m++) {  //square deviation
     evvalue[m] = pow(evvalue[m],2); 
   }
   Double_t sum2;
   for (Int_t m=0; m<lev; m++) { //add up squared deviations
     if (m==0) sum2=0;
     else sum2 += evvalue[m];
   }
   sigmaev = (sqrt(sum2/(lev-1))); //divide by one less than sample size and take the square root

   cerr << "Std Deviation for this event: " << sigmaev << endl;
   cerr << "Test. This is the length of the event. " << lev << endl;

   return sigmaev;
   }

Double_t pst,pet;
Double_t ps,pe;
Int_t psd, ped;



int pstartend (Double_t arg[], Int_t length, Double_t arg2[], double &pstarttime, double &pendtime, double &pstartvalue, double &pendvalue, int &pendbin) { //start/end bin of pulse
   stdev (arg, length);
   Int_t save=0;
   Double_t t=0;
   if (n>=10) t=2;
   else t=.2;
   pstarttime=0; pendtime=0; pstartvalue=0; pendvalue=0; psd=0; pendbin=0;
   for (Int_t d=4; d<length; d++) { //find a rough beginning of pulse
     if (/*(arg[d+4]>t)and(arg[d+3]>t)and(arg[d+2]>t)and*/(arg[(d+1)]>=t)and(arg[d]>=t)and(arg[d-1]>=t)) 
       {  save=d; pstarttime=arg2[d]; pstartvalue=arg[d]; psd=d; break; }
   }
   for (Int_t d=(save-4); d<save; d++) { //find a more accurate beginning of pulse before the rough one.. first pt that exceeds mean+5sigma
      if (arg[d]>(mean+(5*sigma))and(d!=0)and(d!=1))  
 	    { pstarttime=arg2[d]; pstartvalue=arg[d]; psd=d; break; } 
   }
  reverse (arg,length);
  reverse (arg2,length);

  Int_t firstbin = 0;
  Int_t lastbin = 0;
  if (n<5) { firstbin=length-psd-15; lastbin=(length-nd);}
  else { firstbin=14; lastbin=length; }

  for (Int_t d=firstbin; d<lastbin; d++)
    if (/*(arg[d+4]>t)and(arg[d+3]>t)and(arg[d+2]>t)and*/(arg[(d+1)]>=t)and(arg[d]>=t)and(arg[d-1]>=t)) 
       {  save=d; pendtime=arg2[d]; pendvalue=arg[d]; pendbin=(length-d-1); break; }
       
    for (Int_t d=(save-4); d<save; d++) { //find a more accurate end of pulse before the rough one.. first pt that exceeds mean+5sigma
     if (arg[d]>(mean+(5*sigma)))
             { pendtime=arg2[d]; pendvalue=arg[d]; pendbin=(length-d-1); break;}
    } 

  reverse (arg,length);
  reverse (arg2,length);
return psd;
}

Double_t y5,y6,y2,y4=0;
Int_t t10l,t10r, t50l, t50r=0;


void ppercent (Double_t arg[], Int_t length, Double_t arg2[], int &startfallback, int &endfallback, int &fiftyl, int &fiftyr) { //find 10%/50% of max value on rising/falling edges
  Double_t x1,x3,y1,y3 = 0;
  double pten=.10*n;
  for (int i=psd; i<=nd; i++) { //find 10% of max peak on rising edge
    if ((arg[i]>=pten)and(arg[i-1]<=pten)) { x3=arg[i]; y3=arg2[i]; x1=arg[(i-1)]; y1=arg2[(i-1)]; startfallback=(i-1); break; }
  }
  if (x3==0) {  for (int i=(psd-4); i<=nd; i++) { //if t10l doesn't exist in found pulse, look before it by 4 bins
    if ((arg[i]>=pten)and(arg[i-1]<=pten)) { x3=arg[i]; y3=arg2[i]; x1=arg[(i-1)]; y1=arg2[(i-1)]; startfallback=(i-1); break; }
    } }
  if (x3==0) { for (int i=2; i<=nd; i++) { //if t10l doesn't exist in found pulse or 4 bins before, look to beginning of data
    if ((arg[i]>=pten)and(arg[i-1]<=pten)) { x3=arg[i]; y3=arg2[i]; x1=arg[(i-1)]; y1=arg2[(i-1)]; startfallback=(i-1); break; }
    } }
  y5 = (((pten-x1)*(y3-y1))/(x3-x1))+y1; //linearly interpolate
  
  x3=0;
  for (int i=nd; i<=ped; i++) { //find 10% of max peak on falling edge
    if ((arg[i]<=pten)and(arg[i-1]>=pten)) { x3=arg[i]; y3=arg2[i]; x1=arg[(i-1)]; y1=arg2[(i-1)]; endfallback=i; break; }
  }
  if (x3==0) { for (int i=nd; i<=(ped+4); i++) { 
    if ((arg[i]<=pten)and(arg[i-1]>=pten)) { x3=arg[i]; y3=arg2[i]; x1=arg[(i-1)]; y1=arg2[(i-1)]; endfallback=i; break; }
    } }
  if (x3==0) { for (int i=nd; i<=length; i++) { 
    if ((arg[i]<=pten)and(arg[i-1]>=pten)) { x3=arg[i]; y3=arg2[i]; x1=arg[(i-1)]; y1=arg2[(i-1)]; endfallback=i; break; }
    } }
  y6 = (((pten-x1)*(y3-y1))/(x3-x1))+y1;

  x3=0;
  double pfifty=.5*n;  
  for (int i=psd; i<=nd; i++) { //find 50% of max peak on rising edge
    if ((arg[i]>=pfifty)and(arg[i-1]<=pfifty)) { x3=arg[i]; y3=arg2[i]; x1=arg[(i-1)]; y1=arg2[(i-1)]; fiftyl=i; break; }
  } 
  if (x3==0) {  for (int i=(psd-4); i<=nd; i++) { 
    if ((arg[i]>=pfifty)and(arg[i-1]<=pfifty)) { x3=arg[i]; y3=arg2[i]; x1=arg[(i-1)]; y1=arg2[(i-1)]; fiftyl=i; break; }
    } }
  if (x3==0) { for (int i=2; i<=nd; i++) { 
    if ((arg[i]>=pfifty)and(arg[i-1]<=pfifty)) { x3=arg[i]; y3=arg2[i]; x1=arg[(i-1)]; y1=arg2[(i-1)]; fiftyl=i; break; }
    } }
  y2 = (((pfifty-x1)*(y3-y1))/(x3-x1))+y1;

  x3=0;
  for (int i=nd; i<=ped; i++) {  // find 50% of max peak on falling edge
    if ((arg[i]<=pfifty)and(arg[i-1]>=pfifty)) { x3=arg[i]; y3=arg2[i]; x1=arg[(i-1)]; y1=arg2[(i-1)]; fiftyr=i; break; }
  } 
  if (x3==0) { for (int i=nd; i<=(ped+4); i++) {
    if ((arg[i]<=pfifty)and(arg[i-1]>=pfifty)) { x3=arg[i]; y3=arg2[i]; x1=arg[(i-1)]; y1=arg2[(i-1)]; fiftyr=i; break; }
    } }
  if (x3==0) { for (int i=nd; i<=length; i++) { 
    if ((arg[i]<=pfifty)and(arg[i-1]>=pfifty)) { x3=arg[i]; y3=arg2[i]; x1=arg[(i-1)]; y1=arg2[(i-1)]; endfallback=i; break; }
    } }
  y4 = (((pfifty-x1)*(y3-y1))/(x3-x1))+y1;

  return;
}

double riemann (Double_t arg[], Double_t arg2[]) { //compute the area underneathe a pulse using a trapezoidal riemann sum
  Double_t sum = 0;
  for (Int_t d=psd; d<(ped+1); d++) {
    sum += 2*arg[d];
  }
  sum = (/*((arg2[1]-arg2[0])/2)*/.5*(sum-arg[psd]-arg[ped]));
  cerr << " traparea approx: " << sum << endl;
  return sum;
}

double riemann2 (Double_t arg[], Double_t arg2[]) { //cpmpute the area underneath a pulse using the midpoint sum method
  Double_t tot = 0;
  Double_t y1,y3 = 0;
  for (Int_t d=psd; d<(ped+1); d++) {
    y1=arg[d];
    y3=arg[(d+1)];
    Double_t y2 = (/*((arg2[1]-arg2[0])/2)*/.5*(y3-y1))+y1;
    tot += y2;
  }
  // tot *= (arg2[1]-arg2[0]);
  cerr << "midptarea approx: " << tot << endl;
  return tot;
}

Double_t simp=0;
Double_t area;
double simpsons (Double_t arg[], Double_t arg2[]) { //compute the area underneath a pulse using simpson's rule
  simp=0;
  Int_t end=0;
  if (((ped-psd)%2)==0) { end = (ped-psd)+1; }
  else if (((ped-psd)%2)==1) { end = (ped-psd); }
  
  Double_t add = 0;  
 for (int d=0; d<end; d++) {
  if ((d%2)==0) { add =( 2*arg[(psd+d)]); }
  else if ((d%2)==1) { add = (4*arg[(psd+d)]); }
    simp += add;
  }
 simp -= arg[psd];
 simp -= arg[(psd+(end-1))];
 //simp *= /*(arg2[1]-arg2[0])*/1/3;
 simp *=(9.52/3);
 area=simp;
 cerr << " simparea approx: " << simp << endl <<endl;
 return simp;
}

double total=0; 
double chansum (Double_t arg[]) { //compute the total area underneath the pulses for a channel for any given event
  total +=simp;
  cerr << "Running channel total: " << total << endl;
  return total;
}

//Samples from Melinda's code

Double_t aream;
double integral (Double_t arg[]) {
  aream=0;
  for (int i=psd; i<(ped+1); i++) 
    { aream += arg[i]; }

  cerr << "Melinda's integral: " << aream << endl;
return aream;
}

void FindEventStartAndEnd(Double_t pulseData[], 
				   Int_t numBins, 
				   Double_t threshold, 
				   Int_t& pulseStartPoint, 
				   Int_t& pulseEndPoint) {
	for( int b=pulseStartPoint; b<numBins; b++ ) {
		if( pulseData[b] > threshold ) {
			pulseStartPoint = b;
			while( pulseData[b] >= threshold/10. && b < numBins )
				b++;
			pulseEndPoint = b;
			break;
		}
	}
	cerr << "Melinda's beginning: " << pulseStartPoint << endl;
	cerr << "Melinda's end: " << pulseEndPoint << endl << endl;
}

void FindMaxBinAndPeakValue(Double_t pulseData[],
						   Int_t binMin,
						   Int_t binMax,
						   Int_t& maxBin,
						   Double_t& peakValue) {
	maxBin = 0;
	vector <Int_t> maxBins;
	int lastElement=0;
	for(int b=binMin; b<binMax; b++){
		if(pulseData[b] > pulseData[maxBin])
			maxBin = b;
	}
	maxBins.push_back(maxBin);
	for(int b=binMin; b<binMax; b++){
		if(pulseData[b] == pulseData[maxBin] && b != maxBin) {
			maxBins.push_back(b);
			lastElement = maxBins.size()-1;
		}
	}
	bool SEQ=false;
	for(unsigned int i=0; i<maxBins.size()-1; i++) 
		if( (float)(maxBins[lastElement-i]-maxBins[lastElement-(i+1)]) == 1.) 
			SEQ=true;	
	if(SEQ)	maxBin += maxBins.size()/3;
	peakValue = pulseData[maxBin];

	cerr << "Melinda's peak bin: " << maxBin << endl;
	cerr << "Melinda's peak value: " << peakValue << endl;
}


void FindFractionRiseFallBin(Double_t pulseData[],
							Int_t edgeBin,
							Int_t maxBin,
							Double_t fraction,
							Bool_t RISE, 
							Double_t fractionBin){
	fractionBin=0;
	Float_t fractionValue = fraction*pulseData[maxBin];
	vector<int> bins;
	if(RISE){
		for(int b=edgeBin; b<maxBin; b++){
			if(pulseData[b] >= (fractionValue-75.) &&
			   pulseData[b] <= (fractionValue+75.)){
				bins.push_back(b);
			}
		}
	}
	else {
		for(int b=maxBin; b<edgeBin; b++){
			if(pulseData[b] >= (fractionValue-75.) &&
			   pulseData[b] <= (fractionValue+75.)){
				bins.push_back(b);
			}
		}
	}
	for(unsigned int i=0; i<bins.size(); i++) fractionBin+=bins[i];
	if(bins.size()) fractionBin/=bins.size();

	cerr << "Melinda's fraction bin: " << fractionBin << endl;
	
}


//End of Melinda's code


double meanofevent;
double sigmaofevent;
//Some of Justin's code that I added to
void mean_sigma_of_event(Lux_EVT_Event* event, double &meanofev, double&sigmaofev) {  // Sum of all pulses from all channels

  Lux_EVT_Pulse *pulse;
  Int_t start, end;
  start = 0;
  end = 0;
  TH1I* pdata;

  //Get the total time of event 
  for (ULong_t p = 0;p < event->Get_numPulses(); p++){		// scan through this event's pulses
    pulse = event->Get_pulse(p);
    Int_t pstart = pulse->Get_pulseStarts();
    //cerr << "pulse " << p << " starts: " << pstart << endl;
    Int_t pend = pstart + pulse->Get_pulseLength();
    //cerr << "1st length: " << pulse->Get_pulseLength() << endl;
    if (start > pstart) {
      start = pstart;
    }
    if (end < pend) {
      end = pend;
    }
  }

  Lux_EVT_Channel* channel;
  channel = event->Get_channel(0);
  Double_t voltRes = channel->Get_voltageRes();
  Double_t timeRes = channel->Get_timeRes();
  Int_t totalTime = end - start;
  Int_t chan_sum_length = totalTime + 5;
  TH1D* chan_sums = new TH1D("chan_sums","All Channels", chan_sum_length, start*timeRes, end*timeRes);
  for (ULong_t p = 0;p < event->Get_numPulses(); p++){    // scan through this event's pulses again
    pulse = event->Get_pulse(p);
    Int_t offset = pulse->Get_pulseStarts() - start;
    Int_t length = pulse->Get_pulseLength();
    pdata = pulse->Get_pulseData();
    
    Int_t baseline = pulse->Get_pulseBaseline();
    baseline *= -1;
    for (Int_t i = 0; i < length; i++) { // Add pulse data to sum
      if (pdata->GetBinContent(i) != baseline)
        chan_sums->AddBinContent(i + offset, -(pdata->GetBinContent(i))*voltRes);
    }
  }
double sumev;
 int numentries;

  for (int m=0; m<end; m++) {
    if( (chan_sums->GetBinContent(m))!=0 ) {
      sumev += chan_sums->GetBinContent(m); numentries++; }
   }

  meanofev = (sumev/numentries);
  cerr << "Mean for this event: " << meanofev << endl;

  for (Int_t m=0; m<end; m++) { //compute deviation 
        if( chan_sums->GetBinContent(m)!=0 )
    chan_sums->SetBinContent(m, chan_sums->GetBinContent(m)-meanofev);
   }
   for (Int_t m=0; m<end; m++) {  //square deviation
        if( chan_sums->GetBinContent(m)!=0 )
     chan_sums->SetBinContent(m, pow(chan_sums->GetBinContent(m),2));
   }
   Double_t sum2;
   for (Int_t m=0; m<end; m++) { //add up squared deviations
        if( chan_sums->GetBinContent(m)!=0 )
	  sum2 += chan_sums->GetBinContent(m);
   }
  sigmaofev = (sqrt(sum2/(numentries-1))); //divide by one less than sample size and take the square root

   cerr << "Std Deviation for this event: " << sigmaofev << endl;
   cerr << "Test. This is the length of the event. " << numentries << endl;


  
  return;
}


int main(int argc, char* argv[]){    // arguments are the input root file name
  char *inFileName;
  inFileName = argv[1];
  TFile f(inFileName);				// open the input root file
  TTree *T = (TTree*)f.Get("Lux_EVT_Tree");	// get the data tree from the root file
  TCanvas *c2 = new TCanvas("c2","c2");

  Lux_EVT_Event *event = new Lux_EVT_Event();		// create the necessayr Lux_EVT class objects
  Lux_EVT_Header *header = new Lux_EVT_Header();
  Lux_EVT_Channel *channel = new Lux_EVT_Channel();
  Lux_EVT_Pulse *pulse = new Lux_EVT_Pulse();
  char filename[100]; 

  Double_t eventvalues[15000];
  int nump0;

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

  Int_t numberofPulses; 
  Double_t histarea[15000]; //initialize array that will be used to create a histogram of the areas for this run
  Double_t histareasimp[15000];
  Int_t histt0[15000]; //used to create histogram of pulse t0 bins
  Int_t histt2[15000]; //used to create histogram of pulse t2 bins
  Int_t histt1[15000];
  Int_t histt10l[15000];
  Int_t histt10r[15000];
  Int_t histt50l[15000];
  Int_t histt50r[15000];

  Int_t numchan;
  Double_t evtchansum[15000];

  Int_t nument;
  Double_t histmean[15000]; //used to create histogram of event means
  Double_t histsigma[15000]; //used to create histogram of event std deviations
  Double_t histpeakheight[15000];
  

  for (Int_t ev=0;ev<nentries;ev++) {			// scan through the tree entries
    T->GetEntry(ev);					// get the current entries data
    //  if (ev) continue; //dump first event only





    // mean_sigma_of_event(event,meanofevent,sigmaofevent);





    for (ULong_t ch = 0;ch < event->Get_numChannels(); ch++){		// scan through this events channels
      //channel->Clear();
      channel = event->Get_channel(ch);					// fill the Lux_EVT_Channel instance channel with this channels data
    }

    for (ULong_t p = 0;p < event->Get_numPulses(); p++){		// scan through this events pulses
      // pulse->Clear();
      pulse = event->Get_pulse(p);					// fill the Lux_EVT_Pulse instance pulse with this pulses data
      Double_t voltRes = channel->Get_voltageRes();                     //get the voltage resolution and time resolution for this pulse's channel
      Double_t timeRes = channel->Get_timeRes();
      pulseStarts = pulse->Get_pulseStarts();
      pulseLength = pulse->Get_pulseLength();

      if (pulseLength==0) { cerr << " Pulse:           " <<p<< " of "<<event->Get_numPulses()<< " contains no data. " << endl; nump0++; continue; }

      Double_t ptime[pulseLength]; //initialize arrays that will hold the pulse data
      Double_t pdata[pulseLength];
      Int_t bin[pulseLength];  //fill an array with the bin numbers
      for (Int_t i=0; i<pulseLength; i++) {
	bin[i]=i; }

      Double_t channelgain[4]; //channel's gain for the data set I am running this on, how to generalize??
      channelgain[0]=4.04; channelgain[1]=4.13; channelgain[2]=4.07; channelgain[3]=4.09;

      for (Int_t d=0; d<pulseLength; d++){ //fill arrays with pulse data
	if ((pulse->Get_pulseData()->GetBinContent(d))!=0)
	  {pdata[d] = -((2000/pow(2,14))*(pulse->Get_pulseData()->GetBinContent(d)));}//-(1000*(voltRes)*(pulse->Get_pulseData()->GetBinContent(d)));}
        else {pdata[d] = 0;}
       }
      for (Int_t t = 0; t<pulseLength; t++){
        ptime[t]=((timeRes)*(pulseStarts+t));
      }




       cerr<<" Pulse:           "<<p<<" of "<<event->Get_numPulses()<<endl;  	//display the data
      //pulse->Get_pulseData()->Draw();		


      /* TH1F *pulsetrace = new TH1F("pulse","pulse trace", pulseLength, timeRes*pulseStarts, (timeRes*pulseStarts+timeRes*pulseLength));

  for (Int_t i = 0; i < pulseLength; i++) { 
     if ((i!=0)and(i!=1))
      pulsetrace->SetBinContent(i, pdata[i]); //***
  }
 
  pulsetrace->Draw();   // draw the pulse trace for this pulse*/
      sprintf(filename,"plots/pulse_hist_%d_%d_%d.gif",(int)event->Get_eventNum(),(int)pulse->Get_channelNumber(),(int)pulse->Get_pulseNumber());
      /*c2->Print(filename);			// save the pulse trace for this pulse.
      cerr << "File " << filename << " created."<<endl;
      delete pulsetrace;

      pulse->Get_pulseData()->SetBinContent(1,0);
      pulse->Get_pulseData()->Draw();
      sprintf(filename,"plots2/pulse_hist_%d_%d_%d.gif",(int)event->Get_eventNum(),(int)pulse->Get_channelNumber(),(int)pulse->Get_pulseNumber());
      c2->Print(filename);			// save the pulse trace for this pulse.*/
      cerr << "File " << filename << " created."<<endl; 
	
      cerr<<"     pulseNumber: "<<pulse->Get_pulseNumber()<<endl;
      cerr<<"     pulseStarts: "<<pulse->Get_pulseStarts()<<endl;
      cerr<<"     pulseLength: "<<pulse->Get_pulseLength()<<endl;
      cerr<<"   pulseBaseline: "<<pulse->Get_pulseBaseline()<<endl;
      cerr<<"   channelNumber: "<<pulse->Get_channelNumber()<<endl;
      // if ((sizeof(pdata)/sizeof(pulseLength))==0) { cerr << "pulse contains no data" << endl; /*nump0++;*/ continue; }//if the pulse contains no data, skip the rest of the steps and count up the number of these empty pulses


      

  peakfind (pdata,pulseLength,ptime);
  cerr << "        pulsePeak: " << n <<endl;
  cerr << "     pulsePeak'Bin': " << nd <<endl<<endl;
  pstartend (pdata,pulseLength,ptime,pst,pet,ps,pe,ped);
  //cerr << " pulse start bin:  " << pulse->Get_pulseData()->FindBin(ps/(voltRes*1000)) << endl;
  ppercent (pdata,pulseLength,ptime,t10l,t10r,t50l,t50r);
  if ((t10l<psd)or(psd<=0)) { psd=t10l; pst=ptime[t10l]; }
  if ((ped<t10r)or(ped==0)) { ped=t10r; pet=ptime[ped]; }
  if (ped>(pulseLength-1)) { ped=t10r; pet=ptime[ped]; }
   cerr << "       apstartbin: " << psd << " time " << pst << endl;
   cerr << "         apendbin: " << ped << " time " << pet << endl << endl;
   cerr << "   10percentTime1: " << y5 << endl;
   // cerr << "    1st bin above: " << pulse->Get_pulseData()->FindFirstBinAbove(.1*(n/(voltRes*1000)),1) << endl;
   // cerr << "interpolate:       " << pulse->Get_pulseData()->Interpolate(.1*(n/(voltRes*1000))) << endl;
   cerr << "   10percentTime2: " << y6 << endl;
   if ((y6==y5)or(sizeof(y5)==0)) { cerr << "10% does not exist on falling edge" << endl; }
   cerr << "   50percentTime1: "<< y2 << endl;
   cerr << "   50percentTime2: "<< y4 << endl;
  if (y4==y2) { cerr << "50% does not exist on either rising or falling edge" << endl; }
  riemann (pdata,ptime);
  riemann2 (pdata,ptime);
  simpsons (pdata,ptime);  

  unsigned int chan;
  if ((pulse->Get_channelNumber())==chan) chansum(pdata);
  else if ((pulse->Get_channelNumber())!=chan) { total=0; numchan++; chansum(pdata); evtchansum[numchan]=total; }

  chan=pulse->Get_channelNumber();


  /* if ((pulse->Get_channelNumber())==0) {
     if (isnan(simp)) { histarea[numberofPulses]=0; }
     else //histarea[numberofPulses]=simp;*/


        //}
  /*histt0[numberofPulses]=psd+pulseStarts+4970;
  histt2[numberofPulses]=ped+pulseStarts+4970;  
  histt10l[numberofPulses]=t10l+pulseStarts+4970;
  histt10r[numberofPulses]=t10r+pulseStarts+4970;
  histt50l[numberofPulses]=t50l+pulseStarts+4970;
  histt50r[numberofPulses]=t50r+pulseStarts+4970;*/

  numberofPulses++;

  cerr << " editt0: " << psd+pulseStarts+4970 << endl;
  cerr << " editt2: " << ped+pulseStarts+4970 << endl;
  cerr << " editt10l: " << t10l+pulseStarts+4970 << endl;
  cerr << " editt10r: " << t10r+pulseStarts+4970 << endl;
  cerr << " editt50l: " << t50l+pulseStarts+4970 << endl;
  cerr << " editt50r: " << t50r+pulseStarts+4970 << endl;

  //Run Melinda's code
  integral(pdata);
  FindMaxBinAndPeakValue (pdata,psd,ped,nd,n);
  FindEventStartAndEnd (pdata,pulseLength,.00035,psd,ped);

  //convert area in mV*ns -> phe
  cerr << " Area*9.52: " << (aream*9.52) << endl;
  cerr << " Area*9.52/gain: " << ((aream*9.52)/(channelgain[chan])) << endl;
  cerr << " (Area*9.52/gain)/.333 " << (((aream*9.52)/(channelgain[chan]))/.333) << endl;
  cerr << " (Area*9.52/gain)*.333 "  << (((aream*9.52)/(channelgain[chan]))*.333) << endl;
  histareasimp[numberofPulses] = (area/(channelgain[chan]));
  histarea[numberofPulses]= (((aream*9.52)/(channelgain[chan]))*.333);
  cerr << "What's in the array: " << histarea[numberofPulses] << endl << endl;

  
      int lengthev;
      if (p==0) { lengthev=0;
         for (int m=0; m<pulseLength; m++) {
            lengthev++;
            if ((pdata[m]!=0)and(m!=1)) {
	      eventvalues[lengthev]=0;
              eventvalues[lengthev] = pdata[m]; }
            else if ((pdata[m]==0)) eventvalues[lengthev]=0; //lengthev--; 
            else if ((m==1)) eventvalues[lengthev]=0;
         }
      }
      else if (p!=0) {
         for (int m=0; m<pulseLength; m++) {
            lengthev++;
            if ((pdata[m]!=0)and(m!=1)) {
              eventvalues[lengthev]=0;
	      eventvalues[lengthev] = pdata[m]; }
            else if ((pdata[m]==0)) eventvalues[lengthev]=0;//lengthev--; 
            else if ((m==1)) eventvalues[lengthev]=0;
         }
	 }

  //find minimum t0 of the event, save t10s, t2s, t50s of each event
  int t0ev[event->Get_numPulses()];
  t0ev[p]=psd+pulseStarts;
  int t10lev[event->Get_numPulses()];
  int t10rev[event->Get_numPulses()];
  int t2ev[event->Get_numPulses()];
  int t50lev[event->Get_numPulses()];
  int t50rev[event->Get_numPulses()];
  int t1ev[event->Get_numPulses()];
  
  int peakheightev[event->Get_numPulses()];
  peakheightev[p]=n;

  t10lev[p]=t10l;
  t10rev[p]=t10r;
  t2ev[p]=ped;
  t50lev[p]=t50l;
  t50rev[p]=t50r;
  t1ev[p]=nd;
  int mint0=10000;
  int minp=0;

  if (p==((event->Get_numPulses())-1)) { 
  nument++;
  meanevent(pdata,pulseLength,eventvalues,lengthev);
  sigmaevent(pdata,pulseLength,eventvalues,lengthev);
  histmean[nument] = (meanev/(voltRes*1000)); histsigma[nument] = (sigmaev/(voltRes*1000)); 
  histpeakheight[nument] = (n*((pow(2,14))/2000))+200;
  for (unsigned int i=0; i<event->Get_numPulses(); i++) {
    if (t0ev[i]<mint0) mint0=t0ev[i]; minp=i; }
  //histt0[nument]=mint0+4970;
  histt10l[nument]=t10lev[minp]+pulseStarts+4970;
  histt10r[nument]=t10rev[minp]+pulseStarts+4970;
  histt2[nument]=t2ev[minp]+pulseStarts+4970;
  histt50l[nument]=t50lev[minp]+pulseStarts+4970;
  histt50r[nument]=t50rev[minp]+pulseStarts+4970;
  histt0[nument]=psd+pulseStarts+4970;
  histt1[nument]=t1ev[minp]+pulseStarts+4970;
  //histpeakheight[nument] = (peakheightev[minp]*(pow(2,14)/2000))+200;
  }


  for (Int_t i = 0; i<pulseLength; i++){   //print arrays to screen
    cerr<<  bin[i] << "\t" << ptime[i] << "\t" << pdata[i] << endl;
  }

    }} 

  //end of analyzation

  cerr << "This is the number of p0s: " << nump0 << endl;
  //create histograms of different rq1s
  
  TH1F *harea = new TH1F("harea","Histogram of areas",10000,-1,99);
  TH1F *hareasimp = new TH1F("hareasimp","Histogram of simp areas",5100,-1,50);

 TH1F *hmean = new TH1F("hmean","Histogram of event means", 100, -5, 1900);
 TH1F *hsigma = new TH1F("hsigma","Histogram of event std deviations", 100, -5, 3500);
 TH1F *ht0 = new TH1F("t0","Histogram of pulse t0s",83,-100,10000);
 TH1F *ht2 = new TH1F("t2","Histogram of pulse t2s",100,-1000,11000);
 TH1F *ht10l = new TH1F("t10l","Histogram of pulse t10ls",83,0,10000);
 TH1F *ht10r = new TH1F("t10r","Histogram of pulse t10rs",83,0,10000);
 TH1F *ht50l = new TH1F("t50l","Histogram of pulse t50ls",83,0,10000);
 TH1F *ht50r = new TH1F("t50r","Histogram of pulse t50rs",83,0,10000);
 TH1F *eventchansum = new TH1F("eventchansum","Histogram of channel totals per event", 100, -40, 7000);
 TH1F *peakheight = new TH1F("peakheight","Histogram of pulse peaks in mV", 100, -1200, 13600);
 TH1F *ht1 = new TH1F("t1","Histogram of pulse t1s",100,-1000,11000);

 for (int d=0; d<numchan; d++) {
   eventchansum->Fill(evtchansum[d]);
 }

 for (Int_t d=0; d<numberofPulses; d++) {
   if (histarea[d]==0) continue;
   else  {
     harea->Fill(histarea[d]);
     hareasimp->Fill(histareasimp[d]);
 }
 }
  
 cerr << "These are the areas in phe" << endl;
 for (Int_t d=0; d<numberofPulses; d++) {
   cerr << histarea[d] << endl; }

 for (int d=0; d<nument; d++) { 
   hmean->Fill(histmean[d]);
   hsigma->Fill(histsigma[d]);
 }

 cout << "histt0\thistt2\thst10l\thst10r\thst50l\thst50r\thistt1" << endl;
 for (int d=0; d</*numberofPulses*/nument; d++) {
   ht0->Fill(histt0[d]);
   ht2->Fill(histt2[d]);
   ht10l->Fill(histt10l[d]);
   ht10r->Fill(histt10r[d]); 
   ht50l->Fill(histt50l[d]);
   ht50r->Fill(histt50r[d]);
   ht1->Fill(histt1[d]);
 cout << histt0[d] <<"\t"<<histt2[d]<<"\t"<<histt10l[d]<<"\t"<<histt10r[d]<<"\t"<<histt50l[d]<<"\t"<<histt50r[d]<<"\t"<<histt1[d] << endl;
 }

 for (int d=0; d<numberofPulses; d++) { 
   peakheight->Fill(histpeakheight[d]);
 }
 
 TFile rootpeakheight("histpeakheight.root","RECREATE");
 peakheight->Write("histpeakheight.root");
 
 TFile roott1bin("histt1bin.root","RECREATE");
 ht1->Write("histt1bin.root");
 
 TFile rootarea("histareaconv.root","RECREATE");
 harea->Write("histareaconv.root");

 TFile rootareasimp("histareasimpconv.root","RECREATE");
 hareasimp->Write("histareasimpconv.root");
 
  TFile rootmean("histmean.root","RECREATE");
 hmean->Write("histmean.root");
 
  TFile rootstdev("histstddev.root","RECREATE");
 hsigma->Write("histstddev.root");

 TFile roott0bin("histt0bin.root","RECREATE");
 ht0->Write("histt0bin.root");

 TFile roott2bin("histt2bin.root","RECREATE");
 ht2->Write("histt2bin.root");

 TFile roott10lbin("histt10lbin.root","RECREATE");
 ht10l->Write("histt10lbin.root");

 TFile roott10rbin("histt10rbin.root","RECREATE");
 ht10r->Write("histt10rbin.root");

 TFile rootchansum("histchansum.root","RECREATE");
  eventchansum->Write("histchansum.root");
 
 TFile roott50lbin("histt50lbin.root","RECREATE");
  ht50l->Write("histt50lbin.root");
 
 TFile roott50rbin("histt50rbin.root","RECREATE");
 ht50r->Write("histt50rbin.root"); 

 /* TCanvas *c12 = new TCanvas("c12","Histogram of pulse t50r bins",200,10,600,400);
 c12->cd();
 ht50r->Draw();
 c12->Print("histogram24.gif");*/



 }
