#include <stdlib.h>

#include "Riostream.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TH1D.h"

#include "LUX_rq1_header.h"
#include "LUX_rq1_event.h"
#include "LUX_rq1_pulse.h"
#include "LUX_rq1_channel.h"

#include <iostream>
#include <math.h>

using namespace std;

int main(int argc, char* argv[]){   

  if (argc==2){
    
  TFile f(argv[1]);			       
  TTree *T = (TTree*)f.Get("LUX_rq1_tree");	
  TCanvas *c1 = new TCanvas("c1","c1");

  LUX_rq1_header *header = new LUX_rq1_header();	
  LUX_rq1_event *event = new LUX_rq1_event();
  LUX_rq1_pulse *pulse = new LUX_rq1_pulse();
  LUX_rq1_channel *channel = new LUX_rq1_channel();

  Char_t title[50];

  T->SetBranchAddress("Header", &header);
  T->SetBranchAddress("Event", &event);

  //accessing the header information
  T->GetEntry(0);
  printf("There are %i events in the file.\n",header->Get_nb_evts_in_file());
  sprintf(title,"Run %s",header->Get_dataset_name().c_str());
  printf("%s\n",title);

  //access the event values individually
  T->GetEntry(4);
  printf("The mean of pulse 2 in event 4 is %3.1f\n",event->Get_pulse(2)->Get_t_mean());
  T->GetEntry(8);
  printf("The mean of pulse 1 in event 8 is %3.1f\n",event->Get_pulse(1)->Get_t_mean());
  //access the pulse values
  T->GetEntry(15);
  printf("The head value of channel 1 in pulse 2 of event 15 is %3.1f\n",event->Get_channel(2,1)->Get_head());

  //show your options for just what is stored in the tree
  T->GetEntry(10);
  T->Show();

  //now make some interesting plots
  //first if you want to use ROOT's awesome plotting tools (highly recommended)
  c1->Divide(1,3);
  c1->cd(1);
  T->Draw("evt_mean","evt_mean>0");
  c1->cd(2);
  T->Draw("evt_mean:evt_std","evt_mean>0");
  //but if you have to do it by hand for some reason
  //that is possible too, but not as recommended
  TH1D *nonsense_hist = new TH1D("nonsense_hist",title,1000,0,40);
  for (UInt_t i=1;i<header->Get_nb_evts_in_file();i++){
    T->GetEvent(i);
    for (UInt_t j=0;j<event->Get_num_pulses();j++){
      nonsense_hist->Fill(event->Get_pulse(j)->Get_evt_mean()/event->Get_pulse(j)->Get_evt_std()*event->Get_channel(j,1)->Get_head());
    }
  }
  nonsense_hist->SetTitle("Ridiculous Histogram;Random Value;");
  c1->cd(3);
  gPad->SetLogy(1);
  nonsense_hist->Draw();
  
  c1->Print("test_plot.png");
  }
  else printf("Figure this out and come back\n");
}
