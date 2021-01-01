//////////////////////////////////////////////////////////////////////////
// Based on KRG's ProfileLikelihood.C
// AC 2013-07-16 -- Put in simple linear yields of initial quanta and Gaussian resolutions on S1 and S2. 
// Reverted to a beta=0 version of WIMP spectrum until normalisation of full version is fixed.
// KRG - Changing parameterization of NEST to S2/S1 vs S1 and using S1-based energy scale 
//////////////////////////////////////////////////////////////////////////
//
// Setting up an extended maximum likelihood fit
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "TROOT.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TStopwatch.h"
#include "RooAddPdf.h"
#include "RooCategory.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooEfficiency.h"
#include "RooExtendPdf.h"
#include "RooFFTConvPdf.h"
#include "RooGaussian.h"
#include "RooNumGenConfig.h"
#include "RooPlot.h"
#include "RooProfileLL.h"
#include "RooProjectedPdf.h"
#include "RooRandom.h"
#include "RooSimultaneous.h"
#include "RooWimpSpectrum.h"
#include "RooWorkspace.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodIntervalPlot.h"

#include "RooStats/FrequentistCalculator.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/HypoTestPlot.h"

#include "RooStats/NumEventsTestStat.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/SimpleLikelihoodRatioTestStat.h"
#include "RooStats/RatioOfProfiledLikelihoodsTestStat.h"
#include "RooStats/MaxLikelihoodEstimateTestStat.h"

#include "RooStats/HypoTestInverter.h"
#include "RooStats/HypoTestInverterResult.h"
#include "RooStats/HypoTestInverterPlot.h"

#include "ProfileLikelihood.h"

#include "RooWimpSpec_ctsPerkeVkgdaypb.h"

#include <string>

using namespace std;

using namespace RooFit ;
using namespace RooStats;

void PLRun3_tritium(const char* fname=NULL)
{
    
  // to time the macro
  TStopwatch t;
  t.Start();

  // to randomise seed
  RooRandom::randomGenerator()->SetSeed(0);
  
  // Create a local workspace to build the full model
  RooWorkspace *w = new RooWorkspace("w");

//  if(!fname) {
//    cout << "Need to input data file." << endl;
//    return;
//  }

  // Read in Wimp search data
  RooDataSet* searchData;
  TFile* f = TFile::Open(fname);
  TTree *tree = (TTree*)f->Get("t1");
  searchData = new RooDataSet("searchData","searchSata",RooArgSet(S1,logS2S1,r,z),Import(*tree)) ;

  TH1D *histS1 = (TH1D*)f->Get("_hS1");
  RooDataHist dhS1("dhS1","dhS1",S1,Import(*histS1)) ;
  RooHistPdf tritiumSpectrum("tritiumSpectrum","tritiumSpectrum", S1, dhS1);
  w->import(tritiumSpectrum);

  //////////////////////////////////////////////////////////////////////////////////
  // Build the model
  //////////////////////////////////////////////////////////////////////////////////

  //w->import(S1);
  w->import(logS2S1);
  w->import(r);
  w->import(z);
  
  //build ER functions for median log10(S2/S1) (aka logS2S1) as function of S1. 
  //The numerical values are from fits to NEST sims after applying 14%LC, 30%SPEres, 60%CC, 20%SEres.
  w->factory("cexpr::meanEr('ermA*pow(S1,ermB)+ermC*S1',ermA[2.6,2.5,2.7],S1,ermB[-0.143,-0.16,-0.12],ermC[0])");
  //w->factory("cexpr::meanEr('ermA*pow(S1,ermB)+ermC*S1',ermA[2.61,2.60,2.62],S1,ermB[-0.143,-0.145,-0.141],ermC[0])");
  //w->factory("cexpr::meanEr('ermA*pow(S1,ermB)+ermC*S1',ermA[2.61],S1,ermB[-0.143],ermC[0])");
  
  //build ER functions for gaussian sigma in log10(S2/S1) as function of S1 (via parametrisation of mu+1sigma)
  w->factory("cexpr::sigmaEr('ersA*pow(S1,ersB)-ermA*pow(S1,ermB)+(ersC-ermC)*S1',ermA, ermB, ermC,ersA[2.71,2.50,2.90],S1,ersB[-0.141,-0.16,-0.12],ersC[0])");
  //w->factory("cexpr::sigmaEr('ersA*pow(S1,ersB)-ermA*pow(S1,ermB)+(ersC-ermC)*S1',ermA, ermB, ermC,ersA[2.71,2.70,2.71],S1,ersB[-0.141,-0.143,-0.139],ersC[0])");
  //w->factory("cexpr::sigmaEr('ersA*pow(S1,ersB)-ermA*pow(S1,ermB)+(ersC-ermC)*S1',ermA, ermB, ermC,ersA[2.71],S1,ersB[-0.141],ersC[0])");
  
  // Create pdfs P(log10(S2/S1)|S1)
  w->factory("Gaussian::bandEr(logS2S1,meanEr,sigmaEr)");

  // Declare P(r) and P(z) for WIMPs  - only need to do this once for all other uniform pdfs
  w->factory("CEXPR::uniformSpatial('TMath::TwoPi()*r',r,z)");

  // ER background from detector materials
  w->factory("PROD::modelUnNorm(bandEr|S1, tritiumSpectrum, uniformSpatial)"); // P(log10(S2/S1))P(S1)P(r,z)

  w->factory("nSig[1000,100,10000]"); // number of signal events 
  w->factory("ExtendPdf::model(modelUnNorm, nSig)");

  w->Print();

  // Define some useful sets for later  
  w->defineSet("obs","S1,logS2S1,r,z");

  //Save the MC truth model, to be loaded for each pseudoexperiment
  w->saveSnapshot("initialConditions",w->allVars());
  
  //////////////////////////////////////////////////////////////////////////////////
  // Make fit projections if desired
  //////////////////////////////////////////////////////////////////////////////////
  
  TH1* hh_searchData = searchData->createHistogram("S1,logS2S1",100,100) ;
  
//  // Draw the data generated
//  TCanvas* c0 = new TCanvas("c0","DataHistCanvas",600,400) ;
//  gPad->SetLeftMargin(0.15) ; hh_searchData->GetZaxis()->SetTitleOffset(1.4) ; hh_searchData->Draw("col") ;

  w->pdf("model")->fitTo(*searchData);

  // Create 2D histogram of data and binned dataset
  w->factory("ExtendPdf::modelTest(modelUnNorm, nGen[1000000])");
  
  RooDataSet* genData = w->pdf("modelTest")->generate(*w->set("obs"), RooFit::Name("genData"));
  TH1* hh_genData = genData->createHistogram("S1,logS2S1",60,200) ;

  // Draw the data and P(log10(S2/S1)|S1) 
  TCanvas* c0 = new TCanvas("c0","DataHistCanvas",600,400) ;
  gPad->SetLeftMargin(0.15) ; hh_genData->GetZaxis()->SetTitleOffset(1.4) ; hh_genData->Draw("COLZ") ;
  gPad->SetLeftMargin(0.15) ; hh_searchData->SetMarkerStyle(7) ; hh_searchData->Draw("same") ;
  
  // Plot data and PDF overlaid
  RooPlot* xframe = S1.frame(Title("S1 fit projection")) ;
  searchData->plotOn(xframe,Binning(30)) ; 
  //w->pdf("modelWithConstraints")->plotOn(xframe);
  w->pdf("model")->plotOn(xframe);
  
  RooPlot* yframe = logS2S1.frame(Title("log(S2/S1) fit projection")) ;
  searchData->plotOn(yframe,Binning(90)) ;
  //w->pdf("modelWithConstraints")->plotOn(yframe);
  w->pdf("model")->plotOn(yframe);
  
  RooPlot* rframe = r.frame(Title("r fit projection")) ;
  searchData->plotOn(rframe,Binning(30)) ;
  //w->pdf("modelWithConstraints")->plotOn(rframe);
  w->pdf("model")->plotOn(rframe);
  
  RooPlot* zframe = z.frame(Title("z fit projection")) ;
  searchData->plotOn(zframe,Binning(30)) ;
  //w->pdf("modelWithConstraints")->plotOn(zframe);
  w->pdf("model")->plotOn(zframe);
  
  // Draw the frames on the canvas
  TCanvas* c1 = new TCanvas("c1","FitProjCanvas",800,600) ;
  c1->Divide(2,2);
  c1->cd(1);  gPad->SetLeftMargin(0.15) ; xframe->GetYaxis()->SetTitleOffset(1.4) ; xframe->Draw() ;
  c1->cd(2);  gPad->SetLeftMargin(0.15) ; yframe->GetYaxis()->SetTitleOffset(1.4) ; yframe->Draw() ;
  
  c1->cd(3);  gPad->SetLeftMargin(0.15) ; rframe->GetYaxis()->SetTitleOffset(1.4) ; rframe->Draw() ;
  c1->cd(4);  gPad->SetLeftMargin(0.15) ; zframe->GetYaxis()->SetTitleOffset(1.4) ; zframe->Draw() ;

  // Print out timing info  
  t.Stop();
  t.Print();
}
