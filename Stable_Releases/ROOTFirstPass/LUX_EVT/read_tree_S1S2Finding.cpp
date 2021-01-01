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
#include "TH2D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TSpectrum.h"
#include "TPolyMarker.h"
#include "TF1.h"

#include "Lux_EVT_Header.h"
#include "Lux_EVT_Channel.h"
#include "Lux_EVT_Event.h"
#include "Lux_EVT_Pulse.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
using namespace std;

TH1D* sum_of_channels(Lux_EVT_Event* event);
TH1D* first_deriv(TH1D* data, Double_t dx);

Int_t count_peaks(TH1D* data); // return number of peaks in an event
void get_limits(Int_t n_bins, Int_t center, Int_t bin_threshold, Int_t &lower, Int_t &upper); 

Double_t gains[4] = {1, 0.9782, 0.9926, 0.9878}; // Only for 2216

int main(int argc, char* argv[]){    // arguments are the input root file name
  char *inFileName;
  Int_t outFolder;
  inFileName = argv[1];
  outFolder = atoi(argv[2]);
  TFile f(inFileName);				// open the input root file
  TTree *T = (TTree*)f.Get("Lux_EVT_Tree");	// get the data tree from the root file
  TCanvas *c1 = new TCanvas("c1","c1");

  Lux_EVT_Event *event = new Lux_EVT_Event();		// create the necessayr Lux_EVT class objects
  Lux_EVT_Header *header = new Lux_EVT_Header();
  //Lux_EVT_Channel *channel = new Lux_EVT_Channel();
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
  T->GetEntry(0);

  Lux_EVT_Channel* channel;
  channel = event->Get_channel(0);
  Double_t voltRes = channel->Get_voltageRes();
  Double_t timeRes = channel->Get_timeRes();

  //TH2D* peak_shapes = new TH2D("peak_shapes", "Peak Shapes",256, -0.0000001, 0.0000005, 256, 0, 1.5);
  TH2D* dt_s2 = new TH2D("dt_s2", "S2 vs. Delta t", 128, 0, 0.00003, 120, 0, 70);
  //ifstream list;
  //list.open("list.txt");
  //if (!list) {
    //cerr << "Error opening list" << endl;
    //exit(1);
  //}

  //Int_t ev = 0;
  for (Int_t ev=0;ev<nentries;ev++) {			// scan through the tree entries
  //while (!list.eof()) {
    //list >> ev;
    T->GetEntry(ev);					// get the current entries data
    //  if (ev) continue; //dump first event only
    cerr<<" Event:        "<<ev<<endl;				// display the data
    //cerr<<"    eDateTime: "<<event->Get_eDateTime()<<endl;
    //cerr<<"    eLocation: "<<event->Get_eLocation()<<endl;
    //cerr<<"     eventNum: "<<event->Get_eventNum()<<endl;
    //cerr<<"  numChannels: "<<event->Get_numChannels()<<endl;
    //cerr<<"    numPulses: "<<event->Get_numPulses()<<endl;
    //cerr<<"    eventSize: "<<event->Get_eventSize()<<endl;
    //cerr<<"   recordForm: "<<event->Get_recordForm()<<endl;
    //cerr<<"    timeStamp: "<<event->Get_timeStamp()<<endl<<endl;
    
    if (event->Get_numPulses() > 0) {
      Int_t times_counter = 0;
      Int_t not_new;
      Double_t times[event->Get_numPulses()][2];
      for (UInt_t p = 0; p < event->Get_numPulses(); p++) {
        not_new = 0;
	pulse = event->Get_pulse(p);
	Double_t pstart = pulse->Get_pulseStarts();
	Double_t pend = pstart + pulse->Get_pulseLength();
	pstart *= timeRes;
	pend *= timeRes;
	if (!times_counter) {
	  times[0][0] = pstart;
	  times[0][1] = pend;
          times_counter++;
	} else {
	  for (Int_t i = 0; i < times_counter; i++) {
	    if (pstart >= times[i][0] && pstart <= times[i][1]){
              not_new = 1;
	      if ( pend > times[i][1]) {
	        times[i][1] = pend;
		}
	      continue;
	    } else if (pend >= times[i][0] && pend <= times[i][1]) {
	      times[i][0] = pstart;
              not_new = 1;
	      continue;
	    }
	  }
	  if (!not_new) {
	    times[times_counter][0] = pstart;
	    times[times_counter][1] = pend;
	    times_counter++;
	  }
	}
      }

      TH1D* sum = sum_of_channels(event);
      sum->SetStats(kFALSE);
      for (Int_t time = 0; time < times_counter; time++) {
        //Fits
	TF1* g = new TF1("g", "gaus", times[time][0], times[time][1]);
	g->SetLineColor(kBlue);
	sum->Fit(g, "R+");
	//if (times_counter == 2) peak_shapes->Fill(g->GetParameter(2), g->GetParameter(0) );
	delete g;
      }
      Int_t num_peaks = count_peaks(sum);
      TList* functions = sum->GetListOfFunctions();
      TPolyMarker* pm = (TPolyMarker*)functions->FindObject("TPolyMarker");
      //cerr << "peak locations: " << endl;
      if (pm && pm->GetN() == 2) {
        Double_t delta_t = pm->GetX()[0] - pm->GetX()[1];
	Double_t s2 = sum->Integral((sum->GetNbinsX()/3)*2, sum->GetNbinsX());
	dt_s2->Fill(delta_t, s2);
	//cerr << "****************" << endl;
	//cerr << "Delta t: " << delta_t << ", s2: " << s2 << endl;
	//cerr << "****************" << endl;
      }
      //sum->SetXTitle("Seconds");
      //sum->SetYTitle("Volts");
      sum->Draw();
      sprintf(filename,"plots/2216/%d/event_%d_peaks_%d.gif",(int)outFolder,(int)event->Get_eventNum(),(int)num_peaks);
      c1->Print(filename);
      delete sum;
    }
  }

  // print delta_t vs. s2 graph here
  //dt_s2->Draw("BOX");
  //sprintf(filename, "plots/dt_s2.gif");
  //c1->Print(filename);
  

  //list.close();
  

  // draw peak_shapes
  //peak_shapes->SetXTitle("Sigma");
  //peak_shapes->SetYTitle("Constant");
  //peak_shapes->Draw("BOX");
  //sprintf(filename, "plots/box_peak_shapes.gif");
  //c1->Print(filename);
  //delete peak_shapes;
  

} // int main()

TH1D* sum_of_channels(Lux_EVT_Event* event) {  // Sum of all pulses from all channels

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
        chan_sums->AddBinContent(i + offset, -(pdata->GetBinContent(i))*voltRes*gains[pulse->Get_channelNumber()]);
	// Note that gain values are different for each data set.  Make sure the array has the right values
    }
  }
  return chan_sums;
}

Int_t count_peaks(TH1D* data) { // return number of peaks in an event 
  Int_t n_bins = data->GetNbinsX();
  Int_t lower_limit, upper_limit, bin_threshold;
  if (n_bins < 60) {
   bin_threshold = n_bins / 3;
  } else  bin_threshold = 20;
  TSpectrum* s = new TSpectrum(10);
  Int_t num_s_peaks = s->Search(data,2,"nomarkov",0.015);
  Int_t num_peaks = 0;
  TList* functions = data->GetListOfFunctions();
  TPolyMarker* pm = (TPolyMarker*)functions->FindObject("TPolyMarker");
  TPolyMarker* real_peaks = new TPolyMarker();
  Double_t bin_width = data->GetBinWidth(1);
  Int_t offset = -data->GetBinCenter(0)/bin_width;
  for (Int_t i = 0; i < num_s_peaks; i++) {
    // Get correct bin corresponding to peak
    Int_t correct_bin = 0;
    Int_t peak_bin =(Int_t)((pm->GetX()[i] / bin_width)+0.5) + offset;
    Double_t peak_x = pm->GetX()[i];
    do {
      Double_t bin_center = data->GetBinCenter(peak_bin);
      Double_t bin_left = bin_center - (bin_width/2);
      Double_t bin_right = bin_center + (bin_width/2);
      //cerr << "peak " << i << ": " << peak_x << ", " << pm->GetY()[i] << endl;
      //cerr << "Found Bin range: " << bin_left << " : " << bin_right << endl;
      if (peak_x >= bin_left && peak_x <= bin_right) {
        //cerr << "Correct bin found" << endl;
	correct_bin = 1;
      }
      else if (peak_x < bin_left){
        //cerr << "Bin too high\n";
	peak_bin--;
      }
      else {
        //cerr << "Bin too low\n";
	peak_bin++;
      }
    } while(!correct_bin);
    get_limits(n_bins, peak_bin, bin_threshold, lower_limit, upper_limit);
   // Check local area
    Double_t local_max = -FLT_MAX;
    Double_t left_total = 0.0;
    Double_t right_total = 0.0;
    Int_t local_max_bin = peak_bin;
    for (Int_t j = lower_limit; j < upper_limit; j++) {
      // Get real local max; Spectrum search is sometimes off.
      Double_t content = data->GetBinContent(j);
      if (content > local_max) {
        local_max = content;
	local_max_bin = j;
      } 
    }
    get_limits(n_bins, local_max_bin, bin_threshold, lower_limit, upper_limit);
    for (Int_t j = lower_limit; j < upper_limit; j++) {
      // get totals
      if (j < local_max_bin)
        left_total += data->GetBinContent(j);
      else if (j > local_max_bin)
        right_total += data->GetBinContent(j);
    }
    if (local_max_bin != lower_limit) {
      Double_t left_avg = left_total / (local_max_bin - lower_limit);
      Double_t right_avg = right_total / (upper_limit - local_max_bin);
      if (local_max > 0.002 && local_max > left_avg && local_max > right_avg) {
        // It's a peak.  Make sure it hasn't already been counted
	Int_t repeat = 0;
	for (Int_t p = 0; p < real_peaks->GetN(); p++) {
	  if (local_max == real_peaks->GetY()[p]) {
	    repeat = 1;
	    continue;
	  }
	}
	if (!repeat) {
	  real_peaks->SetNextPoint(data->GetBinCenter(local_max_bin), local_max);
	  num_peaks++;
	}
      }
    }
  }
  data->GetListOfFunctions()->Remove(pm);
  delete pm;
  data->GetListOfFunctions()->Add(real_peaks);
  real_peaks->SetMarkerStyle(23);
  real_peaks->SetMarkerColor(kRed);
  real_peaks->SetMarkerSize(1.1);

  return num_peaks;
}

// helper function for count_peaks
void get_limits(Int_t n_bins, Int_t center, Int_t bin_threshold, Int_t &lower, Int_t &upper) {
    if (center >= (n_bins - bin_threshold)) {
      lower = center - bin_threshold;
      upper = n_bins;
    }
    else if (center >= bin_threshold) {
      lower = center - bin_threshold;
      upper = center + bin_threshold;
    }
    else {
      lower = 0;
      upper = center + bin_threshold;
    }
}


TH1D* first_deriv(TH1D* data, Double_t dx) { // numerically approximate first derivative end point data will be no good
  TH1D* deriv = new TH1D(*data);
  Double_t five_const = 1.0/(12.0 * dx);
  Double_t three_const = 0.5 / dx;
  
  Int_t N = data->GetNbinsX();
  Double_t point; 
  
  // Use three point deriv formula for points 1 and N - 2
  point = three_const * (data->GetBinContent(2) - data->GetBinContent(0));
  deriv->SetBinContent(1, point);
  point = three_const * (data->GetBinContent(N-1) - data->GetBinContent(N-3));
  deriv->SetBinContent((N-2), point);

  // Use five point formula for other points
  for (Int_t i = 2; i < (N-2); i++) {
   point = five_const * (data->GetBinContent(i-2) - 8*data->GetBinContent(i-1)
                          + 8*data->GetBinContent(i+1) - data->GetBinContent(i+2));

   deriv->SetBinContent(i, point);
  }

  // really awful approximation for end points
  point = data->GetBinContent(1) - data->GetBinContent(0);
  deriv->SetBinContent(0, point);
  point = data->GetBinContent(N-1) - data->GetBinContent(N-2);
  deriv->SetBinContent(N-1, point);

  return deriv;
}
