#include "TFile.h"
#include "TChain.h"
#include "TLine.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TMath.h"

#include <iostream>

using namespace std;

void makeTritiumNtuple(const char* input, const char* outfile)
{
  gROOT->Reset();
 
  // Create a histogram
  TH1D* _hS1 = new TH1D("_hS1","S1_{c}", 60, 0., 40.);
  TH1D* _hS2 = new TH1D("_hS2","S2", 200, -4., 4.);
  TH2D* _hS2S1 = new TH2D("_hS2S1","log_{10}(S2/S1) vs S1", 60, 0., 60., 200, 1.0, 3.0);

  _hS2->Sumw2();
  
  // Create an ntuple
  TTree t1("t1","a simple Tree with simple variables");
  double S1, S2, logS2S1, s2radius, drift_time, r, z;
  t1.Branch("S1",&S1,"S1/D");
  t1.Branch("logS2S1",&logS2S1,"logS2S1/D");
  t1.Branch("r",&r,"r/D");
  t1.Branch("z",&z,"z/D");
  
  double mean = 0;
  double sigma = 0;
  const double ermA = 2.57766;
  const double ermB = -1.39154e-01;
  const double ersA = 2.72469;
  const double ersB = -1.36292e-01;

  ifstream fileIn(input);
  if(fileIn.is_open()) {
    while(fileIn >> S1 >> S2 >> drift_time >> s2radius) {
      logS2S1 = TMath::Log10(S2/S1);
      r =s2radius;
      z=5.6 + 49. - 0.1508*drift_time;
      
      cout << "S1 = " << S1 << ", logS2S1 = " << logS2S1 << ", r = " << r << ", z = " << z << endl;
      t1.Fill();
      
      mean = ermA*pow(S1,ermB);
      sigma = ersA*pow(S1,ersB) - ermA*pow(S1,ermB);

      // Fill the histograms
      _hS1->Fill(S1);
      if(r<18. && z< 47. && z > 7. && S1 > 2 && S1 <= 3)
	_hS2->Fill( (logS2S1-mean)/sigma );
      _hS2S1->Fill(S1,logS2S1);
    }
  }
 
  TFile* fileOut = TFile::Open(outfile,"recreate");
  t1.Write();
  _hS1->Write();
  _hS2->Write();
  _hS2S1->Write();
  fileOut->Close();
}
