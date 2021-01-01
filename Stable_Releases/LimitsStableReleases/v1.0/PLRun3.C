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
#include "TProfile.h"
#include "TNtuple.h"
#include "TStopwatch.h"
#include "RooAddPdf.h"
#include "RooCachedPdf.h"
#include "RooCategory.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooEfficiency.h"
#include "RooExtendPdf.h"
#include "RooFFTConvPdf.h"
#include "RooGaussian.h"
#include "RooHist.h"
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

#include "RooMsgService.h"

#include "ProfileLikelihood.h"

#include "RooWimpSpec_ctsPerkeVkgdaypb.h"

#include <string>
#include <fstream>
#include <iostream>
#include <cstdlib>
using namespace std;

using namespace RooFit ;
using namespace RooStats;

//void buildModel(RooWorkspace* w, RooWorkspace* wSig)
void buildModel(RooWorkspace* w)
{
  
  ////////////////////////////////////////////////////////////////////////////////////////
  // The Model building stage
  ////////////////////////////////////////////////////////////////////////////////////////

  //S1.setBins("fft",10000); // Sampling frequency for FFT

//  // Read in signal model workspace
//  w->import(*wSig->var("mWimp"));
//  w->import(*wSig->var("events_per_pbkgday"));
//  w->import(*wSig->var("nrmA"));
//  w->import(*wSig->var("nrmB"));
//  w->import(*wSig->var("nrsA"));
//  w->import(*wSig->var("nrsB"));
//
//  cout << "nrmA = " << w->var("nrmA")->getVal() << ", nrmB = " << w->var("nrmB")->getVal() << endl;
//  cout << "nrsA = " << w->var("nrsA")->getVal() << ", nrsB = " << w->var("nrsB")->getVal() << endl;

  // Import parameters
  w->import(sig_xs);
  //w->import(mWimp);
  //w->import(S1);
  //w->import(logS2S1);

  w->import(targetA);
  w->import(targetZ);
  w->import(T);
  w->import(r);
  w->import(z);
  w->import(rho);
  w->import(dayOfYear);
  w->import(vEsc);
  w->import(v0);
  w->import(vE);
  w->import(ff_c);
  w->import(ff_r0); //not used in Lewin&Smith Helm FF
  w->import(ff_a);
  w->import(ff_s);
  w->import(idm_delta);
  w->import(halo_beta); //not used for SHM
  
  w->import(rateCompton);
  w->import(ra);
  w->import(rd);
  w->import(za);
  w->import(zb);
  w->import(zd);
  
  //w->import(rateActivatedXe);
  w->import(rateActivatedXeA);
  w->import(rateActivatedXeB);

  w->import(rateRn);

  w->import(nCompton);
  w->import(nActivatedXe);
  w->import(nRn);
  
  w->import(relUncertaintyBg);
  w->import(nNrCalib);
  w->import(nErCalib);

  w->import(xeDensity);
  w->import(zLiquid);
  w->import(driftVelocity);
  w->import(eTau);
  
  w->import(thresNe);
  w->import(singleElec);

  // Parameterize the recoil energies in terms of S1
  w->factory("cexpr::Eee('aEee+bEee*S1+cEee*pow(S1,2)+dEee*pow(S1,3)+eEee*pow(S1,4)+fEee*pow(S1,5)',S1,aEee[0.135369],bEee[0.360545],cEee[-0.0186279],dEee[0.00074638],eEee[-1.42213e-05],fEee[1.01752e-07])"); // [keVee]
  w->factory("cexpr::Enr('aEnr+bEnr*S1+cEnr*pow(S1,2)+dEnr*pow(S1,3)+eEnr*pow(S1,4)+fEnr*pow(S1,5)',S1,aEnr[0.135369],bEnr[0.360545],cEnr[-0.0186279],dEnr[0.00074638],eEnr[-1.42213e-05],fEnr[1.01752e-07])"); // [keVnr]
  //w->factory("cexpr::Enr('A*pow(Eee,B)',Eee,A[7.07935],B[0.87101])"); // [keVnr]

  // Terms to transform integral from Int[f(E)dE] to Int[g(E(S1)) dE/dS1 dS1]
  w->factory("cexpr::dEee_dS1('bEee+cEee*S1+dEee*pow(S1,2)+eEee*pow(S1,3)+fEee*pow(S1,4)',S1,bEee,cEee,dEee,eEee,fEee)");
  w->factory("cexpr::dEnr_dS1('bEnr+cEnr*S1+dEnr*pow(S1,2)+eEnr*pow(S1,3)+fEnr*pow(S1,4)',S1,bEnr,cEnr,dEnr,eEnr,fEnr)"); 

  //The numerical values are from fits to NEST sims after applying 14%LC, 30%SPEres, 60%CC, 20%SEres.
  //build NR functions for median log10(S2/S1) (aka logS2S1) as function of S1. 
  //build NR functions for gaussian sigma in log10(S2/S1) as function of S1 (via parametrisation of mu+1sigma)
  //w->factory("cexpr::meanNr('nrmA*pow(S1,nrmB)+nrmC*S1-0.086',nrmA,S1,nrmB,nrmC[0])");
  //w->factory("cexpr::sigmaNr('nrsA*pow(S1,nrsB)-nrmA*pow(S1,nrmB)+(nrsC-nrmC)*S1-0.086',nrmA, nrmB, nrmC,nrsA,S1,nrsB,nrsC[0])");

//  w->factory("cexpr::meanNr('nrmA*pow(S1,nrmB)-0.086',nrmA,S1,nrmB)");
//  w->factory("cexpr::sigmaNr('nrsA*pow(S1,nrsB)-meanNr',nrsA,S1,nrsB,meanNr)");
 
  //build ER functions for median log10(S2/S1) (aka logS2S1) as function of S1. 
  //build ER functions for gaussian sigma in log10(S2/S1) as function of S1 (via parametrisation of mu+1sigma)

  if(!bandsFixed) {
    // krg fit to H-3
    w->factory("cexpr::meanEr('ermA*pow(S1,ermB)+ermC*S1',ermA[2.57766,2.566493,2.5888271],S1,ermB[-1.39154e-01,-1.408161e-01,-1.3749184e-01],ermC[0])");
    w->factory("cexpr::sigmaEr('ersA*pow(S1,ersB)-ermA*pow(S1,ermB)+(ersC-ermC)*S1',ermA, ermB, ermC,ersA[2.72469,2.711174,2.738206],S1,ersB[-1.36292e-01,-1.3819031e-01,-1.3439369e-01],ersC[0])");

    // fit from AD
    //w->factory("cexpr::meanEr('ermA*pow(S1,ermB)+ermC*S1',ermA[2.61,2.60,2.62],S1,ermB[-0.143,-0.145,-0.141],ermC[0])");
    //w->factory("cexpr::sigmaEr('ersA*pow(S1,ersB)-ermA*pow(S1,ermB)+(ersC-ermC)*S1',ermA, ermB, ermC,ersA[2.71,2.70,2.71],S1,ersB[-0.141,-0.143,-0.139],ersC[0])");
  }
  else {
    // krg fit to H-3
    w->factory("cexpr::meanEr('ermA*pow(S1,ermB)+ermC*S1',ermA[2.57766],S1,ermB[-1.39154e-01],ermC[0])");
    w->factory("cexpr::sigmaEr('ersA*pow(S1,ersB)-ermA*pow(S1,ermB)+(ersC-ermC)*S1',ermA, ermB, ermC,ersA[2.72469],S1,ersB[-1.36292e-01],ersC[0])");

    // fit from AD
    //w->factory("cexpr::meanEr('ermA*pow(S1,ermB)+ermC*S1',ermA[2.61],S1,ermB[-0.143],ermC[0])");
    //w->factory("cexpr::sigmaEr('ersA*pow(S1,ersB)-ermA*pow(S1,ermB)+(ersC-ermC)*S1',ermA, ermB, ermC,ersA[2.71],S1,ersB[-0.141],ersC[0])");
  }

  // Create pdfs P(log10(S2/S1)|S1)
  //w->factory("Gaussian::bandNr(logS2S1,meanNr,sigmaNr)");
  w->factory("Gaussian::bandEr(logS2S1,meanEr,sigmaEr)");

  // Define some useful sets for later  
  w->defineSet("obs","S1,logS2S1,r,z");
  //w->defineSet("obs","S1,logS2S1");
  //w->defineSet("nuis","nCompton,nActivatedXe,nRn,nrmA,nrmB,nrmC,nrsA,nrsB,nrsC,ermA,ermB,ermC,ersA,ersB,ersC");
  if(!bandsFixed)
    w->defineSet("nuis","nCompton,nActivatedXe,nRn,ermA,ermB,ersA,ersB");
  else
    w->defineSet("nuis","nCompton,nActivatedXe,nRn");

  ////////////////////////////////////////////////////////////////////////////////////////
  // Efficiency model
  ////////////////////////////////////////////////////////////////////////////////////////

  //S1 efficiency
  w->factory("CEXPR::effS1('0.5*TMath::Erf((S1-aEff)/(sqrt(2.)*bEff))+0.5', S1, aEff[1.48805],bEff[0.484580])");
  //w->factory("effS1[1.]");

  ////////////////////////////////////////////////////////////////////////////////////////
  // Build up model with terms
  // P(log10(S2/S1)|S1)P(S1)P(r,z)
  // for both signal & bg
  ////////////////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////////////////
  // Signal model
  ////////////////////////////////////////////////////////////////////////////////////////
  
  // Declare pdf of WIMP energy spectrum P(S1)
//  gROOT->ProcessLine(".L RooWimpSpectrum.cxx+");//load a custom library
//  w->importClassCode("RooWimpSpectrum", kTRUE);//import custom class into workspace
//  w->factory("WimpSpectrum::dRdE(Enr,mWimp, v0, vE, vEsc, rho, ff_c, ff_r0, ff_a, ff_s, idm_delta, halo_beta)");
//  //w->factory("CEXPR::dRdS1('dRdE*dEnr_dS1*LeffThreshold',dRdE,dEnr_dS1,LeffThreshold)");
//  w->factory("CEXPR::dRdS1('dRdE*dEnr_dS1',dRdE,dEnr_dS1)");

  // Declare P(r) and P(z) for WIMPs  - only need to do this once for all other uniform pdfs
  w->factory("CEXPR::UniformSpatial('TMath::TwoPi()*r',r,z)");

  // Multiply P(log10(S2/S1)|S1)P(S1)P(r)P(z)
//  w->factory("PROD::sigModelUnNorm(bandNr|S1, dRdS1, UniformSpatial)");
//  //w->factory("PROD::sigModelUnNorm(bandNr|S1, dRdS1)");
//  w->factory("PROD::sigEModel(dRdS1, UniformSpatial)");

  w->factory("PROD::sigModelUnNorm(signalHistPdf, UniformSpatial)");

  // Make a variable corresponding to the number of signal events
  double fiducialMass = (w->pdf("UniformSpatial")->createIntegral( RooArgSet(*w->var("r"),*w->var("z")) )->getVal() 
			 * w->var("xeDensity")->getVal());

  cout << "Fiducial mass = " << fiducialMass << endl;
  w->factory("fidMass[1.]");
  w->var("fidMass")->setVal(fiducialMass);
  //w->factory("prod::nSig(sig_xs, T, fidMass, nSigPerPbKgDay)"); // number of signal events 
  //w->var("nSigCubicCmPerPbDay")->setVal(nSigCubicCmPerPbDay);

//  double nSigCubicCmPerKgPbDay =w->pdf("sigEModel")->createIntegral( RooArgSet(*w->var("S1"),*w->var("r"),*w->var("z")) )->getVal(); // [day-1 pb-1 cm3]
//  cout << "double nSigCubicCmPerKgPbDay = " << nSigCubicCmPerKgPbDay << endl;
//
//  w->factory("nSigCubicCmPerKgPbDay[1]");
//  w->var("nSigCubicCmPerKgPbDay")->setVal(nSigCubicCmPerKgPbDay);
//  w->factory("prod::nSig(sig_xs,T, xeDensity, nSigCubicCmPerKgPbDay)"); // number of signal events 
  
  w->factory("nSig[0,0,100]"); // number of signal events 

  w->factory("ExtendPdf::sigModel(sigModelUnNorm, nSig)");

  ////////////////////////////////////////////////////////////////////////////////////////
  // Background model
  ////////////////////////////////////////////////////////////////////////////////////////
  
  // ER background from detector materials
  w->factory("CEXPR::comptonSpectrum('dEee_dS1*effS1',dEee_dS1,effS1)"); // P(S1)
  w->factory("CEXPR::comptonSpatial('TMath::TwoPi()*r*rateCompton*exp(pow(r/ra, 2.)+pow((z-za)/zb, 2.))', r,rateCompton,ra,z ,za,zb)"); // P(r,z)
  w->factory("PROD::comptonModelUnNorm(bandEr|S1, comptonSpectrum, comptonSpatial)"); // P(log10(S2/S1))P(S1)P(r,z)
//  w->factory("CachedPdf::comptonModelUnNorm(comptonModelUnNormUnCached)"); // P(log10(S2/S1))P(S1)P(r,z)

//  w->factory("PROD::comptonModelUnNorm(bandEr, comptonSpatial)"); // P(log10(S2/S1))P(S1)P(r,z)
  w->factory("PROD::comptonEModel(comptonSpectrum, comptonSpatial)"); // P(log10(S2/S1))P(S1)P(r,z)
  
  // Calculate number of expected Compton events due to materials
  RooAbsReal* nComptonCubicCmPerDay = w->pdf("comptonEModel")->createIntegral(RooArgSet(*w->var("S1"),*w->var("r"),*w->var("z"))); // [day-1]
  nComptonCubicCmPerDay->SetName("nComptonCubicCmPerDay");
  w->import(*nComptonCubicCmPerDay);
  w->factory("prod::nComptonPred(T, xeDensity, nComptonCubicCmPerDay)"); // number of Compton events predicted by bg model
  w->factory("prod::sigmaCompton(nComptonPred,relUncertaintyBg)"); // error on number of Compton events predicted by bg model

  //w->factory("ExtendPdf::comptonModel(comptonModelUnNorm, nCompton)");
  
  // ER background from Xe-127
  // Use 2D hist pdf from sims
  w->factory("CEXPR::activatedXeSpatial('TMath::TwoPi()*r*0.47*(exp(-(0.)/52.5)*exp(fabs(r/rd)))+(rateActivatedXeB/rateActivatedXeA)*exp(-0./52.5)*exp(fabs(z-za)/zd)',r,rateActivatedXeA,rd,rateActivatedXeB,z,za,zd)"); // P(r,z)
  w->factory("PROD::activatedXeModelUnNorm(xenonHistPdf, activatedXeSpatial)"); // P(log10(S2/S1))P(S1)P(r,z)

  //w->factory("CEXPR::activatedXeSpatial('TMath::TwoPi()*r*rateActivatedXe',r,rateActivatedXe,z)"); // P(r,z)
  // Include time-dep for later, use average value between 4/22-8/08 for time-averaged rate
  //w->factory("CEXPR::activatedXeSpectrum('(actA1*exp(-0.5*pow((Eee-actE1)/actSig1,2.))+actA2*exp(-0.5*pow((Eee-actE2)/actSig2,2.)))*dEee_dS1*effS1',actA1[2.17],Eee,actE1[0.970],actSig1[0.290],actA2[1.74],actE2[4.91],actSig2[1.473],dEee_dS1,effS1)"); // P(S1)

  //w->factory("PROD::activatedXeModelUnNorm(bandEr|S1, activatedXeSpectrum, activatedXeSpatial)"); // P(log10(S2/S1))P(S1)P(r,z)
  //w->factory("PROD::activatedXeEModel(activatedXeSpectrum, activatedXeSpatial)"); // P(log10(S2/S1))P(S1)P(r,z)
  
  // Calculate number of expected Xe-127 events
  //RooAbsReal* nActivatedXeCubicCmPerDay = w->pdf("activatedXeEModel")->createIntegral(RooArgSet(*w->var("S1"),*w->var("r"),*w->var("z"))); // [day-1]
  //nActivatedXeCubicCmPerDay->SetName("nActivatedXeCubicCmPerDay");
  //w->import(*nActivatedXeCubicCmPerDay);
  //w->factory("prod::nActivatedXePred(T, xeDensity, nActivatedXeCubicCmPerDay)"); // number of Xe-127 events predicted by bg model
  //w->factory("prod::nActivatedXePred(T, fidMass, xenon_per_kgday)"); // number of Xe-127 events predicted by bg model
  w->factory("prod::sigmaActivatedXe(nActivatedXePred,relUncertaintyBg)"); // error on number of Xe-127 events predicted by bg model

  //w->factory("ExtendPdf::activatedXeModel(activatedXeModelUnNorm, nActivatedXe)");
  
  // ER background from Rn-222 and Kr-85
  w->factory("CEXPR::rnSpectrum('dEee_dS1*effS1',dEee_dS1,effS1)"); // P(S1)
  w->factory("CEXPR::rnSpatial('TMath::TwoPi()*r*rateRn',r,rateRn,z)"); // P(r,z)
  w->factory("PROD::rnModelUnNorm(bandEr|S1, rnSpectrum, rnSpatial)"); // P(log10(S2/S1))P(S1)P(r,z)
//  w->factory("CachedPdf::rnModelUnNorm(rnModelUnNormUnCached)"); // P(log10(S2/S1))P(S1)P(r,z)

//  w->factory("PROD::rnModelUnNorm(bandEr, rnSpatial)"); // P(log10(S2/S1))P(S1)P(r,z)
  w->factory("PROD::rnEModel(rnSpectrum, rnSpatial)"); // P(log10(S2/S1))P(S1)P(r,z)
  
  // Calculate number of expected Rn-222 events
  RooAbsReal* nRnCubicCmPerDay = w->pdf("rnEModel")->createIntegral(RooArgSet(*w->var("S1"),*w->var("r"),*w->var("z"))); // [day-1]
  nRnCubicCmPerDay->SetName("nRnCubicCmPerDay");
  w->import(*nRnCubicCmPerDay);
  w->factory("prod::nRnPred(T, xeDensity, nRnCubicCmPerDay)"); // number of Xe-127 events predicted by bg model
  w->factory("prod::sigmaRn(nRnPred,relUncertaintyBg)"); // error on number of Xe-127 events predicted by bg model
  
  //w->factory("ExtendPdf::rnModel(rnModelUnNorm, nRn)");

  // Set the initial number of backgrounds to the prediction
  double nComptonInit =  w->function("nComptonPred")->getVal();
  double nActivatedXeInit = w->function("nActivatedXePred")->getVal();
  double nRnInit = w->function("nRnPred")->getVal();
  
  cout << "nActivatedXeInit = " << nActivatedXeInit << endl;

  w->var("nCompton")->setVal(nComptonInit);
  w->var("nActivatedXe")->setVal(nActivatedXeInit);
  w->var("nRn")->setVal(nRnInit);
  
  // Fix number of background events to predictions
  if(bgNormFixed) {
    w->var("nCompton")->setConstant(kTRUE);
    w->var("nActivatedXe")->setConstant(kTRUE);
    w->var("nRn")->setConstant(kTRUE);
  }

  ////////////////////////////////////////////////////////////////////////////////////////
  // Full model with constraints
  ////////////////////////////////////////////////////////////////////////////////////////
  
  // Put together complete likelihood
  w->factory("SUM::fullModel(nSig*sigModelUnNorm,nCompton*comptonModelUnNorm,nActivatedXe*activatedXeModelUnNorm,nRn*rnModelUnNorm)");
  w->factory("SUM::bgModel(nCompton*comptonModelUnNorm,nActivatedXe*activatedXeModelUnNorm,nRn*rnModelUnNorm)");
  //w->factory("ExtendPdf::bgModel(activatedXeModelUnNorm,nActivatedXe)");

  // Add constraints on nuisance parameters
  w->factory("Gaussian::comptonConstraint(nComptonPred, nCompton, sigmaCompton)");
  w->factory("Gaussian::activatedXeConstraint(nActivatedXePred, nActivatedXe, sigmaActivatedXe)");
  w->factory("Gaussian::rnConstraint(nRnPred, nRn, sigmaRn)");
  
  w->factory("PROD::modelWithConstraints(fullModel,comptonConstraint,activatedXeConstraint,rnConstraint)"); // product of terms

  // Use brute-force accept/reject method
  const RooNumGenConfig* myGen = w->pdf("modelWithConstraints")->getGeneratorConfig();
  myGen->Print();

  w->importClassCode();  
}

void PLRun3(const char* fname=NULL, const int thisWimpMass=0, bool plcToyMC = kFALSE, 
	    bool showProjections = kFALSE, bool hypoTestInv = kFALSE, bool hypoTestNull = kFALSE)
{
    
  // to time the macro
  TStopwatch t;
  t.Start();
  
  // Create a local workspace to build the full model
  RooWorkspace *w = new RooWorkspace("w");

  // Read in Wimp search data
  RooDataSet* searchData;
  if(fname) {
    // Read in search, calib data
    TFile* f = TFile::Open(fname);
    TTree *tree = (TTree*)f->Get("t1");
    searchData = new RooDataSet("searchData","searchSata",RooArgSet(S1,logS2S1,r,z),Import(*tree)) ;
  }

  TFile* f1 = TFile::Open(TString::Format("MCHistInput/wimpDensityFile_%05iGeV_S1_2to30phe.root",thisWimpMass));
  //TFile* f1 = TFile::Open(TString::Format("MCHistInput/wimpDensityFile_%05iGeV_NoEThresh.root",thisWimpMass));
  
  mWimp.setVal(thisWimpMass);
  mWimp.setConstant(kTRUE);
  w->import(mWimp);

  gROOT->cd();

  cout << "thisWimpMass = " << thisWimpMass << ", mWimp = " << w->var("mWimp")->getVal() << endl;

  TH2D *sigHist = (TH2D*)f1->Get("wimpEventDensity");
  double rateFactor = sigHist->Integral(sigHist->FindBin(S1.getMin(),0),sigHist->FindBin(S1.getMax(),0)-1,1,sigHist->GetNbinsY());
  cout << "Rate factor = " << rateFactor << endl;
  RooRealVar events_per_pbkgday("events_per_pbkgday","events_per_pbkgday",rateFactor);
  w->import(events_per_pbkgday);
  RooDataHist sigDH("sigDH","sigDH",RooArgSet(S1,logS2S1),sigHist);
  RooHistPdf signalHistPdf("signalHistPdf","signalHistPdf",RooArgSet(S1,logS2S1),sigDH);
  w->import(signalHistPdf);

  // Read in ER band
//  TFile* f2 = TFile::Open("wSmoothedTritium.root");
//  RooWorkspace* wER = (RooWorkspace*) f2->Get("wSmoothedTritium"); // file with fit to WIMP bands
//  RooHistPdf* bandEr = (RooHistPdf*)wER->pdf("smoothedTritiumHistPdf");
//  bandEr->SetName("bandEr");
//  w->import(*bandEr);

  //TFile* f2 = TFile::Open("smoothedTritiumTH2.root");
//  TFile* f2 = TFile::Open("smoothedTritiumTH2_allbins.root");
//  TH2D *tritiumHist = (TH2D*)f2->Get("kdeTritiumPdf");
//  RooDataHist tritiumDH("tritiumDH","tritiumDH",RooArgSet(S1,logS2S1),tritiumHist);
//  RooHistPdf bandEr("bandEr","bandEr",RooArgSet(S1,logS2S1),tritiumDH);
//  w->import(bandEr);

  // Read in Xe-127 hist pdf
  TFile* f3 = TFile::Open("Xe127ActivationPdf/smoothedNEST_xenon127_eventsPerBinIn117.6kgRun03.root");
  //TFile* f3 = TFile::Open("NEST_xenon127_eventsPerBinIn117.6kgRun03.root");

  TH2D *xenonHist = (TH2D*)f3->Get("xenonHist");
  double rateFactor_actXe = xenonHist->Integral(xenonHist->FindBin(S1.getMin(),0),xenonHist->FindBin(S1.getMax(),0)-1,1,xenonHist->GetNbinsY());
  cout << "Rate factor Xe-127 = " << rateFactor_actXe << endl;

  gROOT->cd();

  RooRealVar nActivatedXePred("nActivatedXePred","nActivatedXePred",rateFactor_actXe);
  w->import(nActivatedXePred);
  RooDataHist xenonDH("xenonDH","xenonDH",RooArgSet(S1,logS2S1),xenonHist);
  RooHistPdf xenonHistPdf("xenonHistPdf","xenonHistPdf",RooArgSet(S1,logS2S1),xenonDH);
  w->import(xenonHistPdf);

  TCanvas* c00 = new TCanvas("c00","SigHistCanvas",1200,400) ;
  c00->Divide(3,1);
  c00->cd(1); w->pdf("signalHistPdf")->createHistogram("S1,logS2S1")->Draw("colz");
  //c00->cd(2); w->pdf("bandEr")->createHistogram("S1,logS2S1")->Draw("colz");
  c00->cd(3); w->pdf("xenonHistPdf")->createHistogram("S1,logS2S1")->Draw("colz");
  c00->Print(TString::Format("sigHist_%05iGeV.pdf",thisWimpMass));

  //set generation to accept/reject to avoid problems with crashes in RooFoamGenerator
  RooAbsPdf::defaultGeneratorConfig()->method1D(kFALSE,kFALSE).setLabel("RooAcceptReject") ;
  RooAbsPdf::defaultGeneratorConfig()->method2D(kFALSE,kFALSE).setLabel("RooAcceptReject") ;
  RooAbsPdf::defaultGeneratorConfig()->methodND(kFALSE,kFALSE).setLabel("RooAcceptReject") ;

  //rott - Need to change *all* the generator defaults
  RooAbsPdf::defaultGeneratorConfig()->method1D(kTRUE,kFALSE).setLabel("RooAcceptReject",kTRUE) ;
  RooAbsPdf::defaultGeneratorConfig()->method2D(kTRUE,kFALSE).setLabel("RooAcceptReject",kTRUE) ;
  RooAbsPdf::defaultGeneratorConfig()->methodND(kTRUE,kFALSE).setLabel("RooAcceptReject",kTRUE) ;

  RooAbsPdf::defaultGeneratorConfig()->method1D(kFALSE,kTRUE).setLabel("RooAcceptReject",kTRUE) ;
  RooAbsPdf::defaultGeneratorConfig()->method2D(kFALSE,kTRUE).setLabel("RooAcceptReject",kTRUE) ;
  RooAbsPdf::defaultGeneratorConfig()->methodND(kFALSE,kTRUE).setLabel("RooAcceptReject",kTRUE) ;

  // to randomise seed
  RooRandom::randomGenerator()->SetSeed(0);
  
  //////////////////////////////////////////////////////////////////////////////////
  // Build the model
  //////////////////////////////////////////////////////////////////////////////////
  
  //buildModel(w,wSig);
  buildModel(w);
  
  //////////////////////////////////////////////////////////////////////////////////
  // Configure model and save to workspace
  //////////////////////////////////////////////////////////////////////////////////
  
  ModelConfig* model = new ModelConfig("model", w);
  model->SetPdf(*w->pdf("modelWithConstraints"));
  model->SetObservables(*w->set("obs"));
  //model->SetParametersOfInterest(*w->var("sig_xs"));
  model->SetParametersOfInterest(*w->var("nSig"));
  model->SetNuisanceParameters(*w->set("nuis"));
  
  w->import(*model);
  
  RooRealVar* poi = (RooRealVar*) model->GetParametersOfInterest()->first();
  
  ModelConfig*  bModel = (ModelConfig*) model->Clone();
  bModel->SetName("model_with_poi_0");
  double oldval = poi->getVal();
  poi->setVal(0);
  bModel->SetSnapshot( *poi  );
  poi->setVal(oldval);
  
  w->import(*bModel);
  
  //Save the MC truth model, to be loaded for each pseudoexperiment
  w->saveSnapshot("initialConditions",w->allVars());
  
  //////////////////////////////////////////////////////////////////////////////////
  // Make fit projections if desired
  //////////////////////////////////////////////////////////////////////////////////

  //w->pdf("modelWithConstraints")->graphVizTree("model.dot");
  w->Print();

  if(showProjections) {
  
    // Generate data
    if(!fname)
      searchData = w->pdf("bgModel")->generate(*w->set("obs"), RooFit::Name("searchData"));
    
    TH1* hh_searchData = searchData->createHistogram("S1,logS2S1",60,200) ;
    
    // Fit model with constraints
    w->pdf("modelWithConstraints")->fitTo(*searchData, Constrain(*w->set("nuis")));

//    // Create 2D histogram of data and binned dataset
//    //RooDataSet* genData = w->pdf("modelWithConstraints")->generate(*w->set("obs"), RooFit::Name("genData"));
//    w->factory("ExtendPdf::modelErTest(comptonModelUnNorm, nGen[1000000])");
//    RooDataSet* genErData = w->pdf("modelErTest")->generate(*w->set("obs"), RooFit::Name("genErData"));
//    TH1* hh_genErData = genErData->createHistogram("S1,logS2S1",60,200) ;
//
//    w->factory("ExtendPdf::modelSigTest(sigModelUnNorm, nGen[1000000])");    
//    RooDataSet* genSigData = w->pdf("modelSigTest")->generate(*w->set("obs"), RooFit::Name("genSigData"));
//    TH1* hh_genSigData = genSigData->createHistogram("S1,logS2S1",60,200) ;
//    
//    // Draw the data and P(log10(S2/S1)|S1) 
//    TCanvas* c0 = new TCanvas("c0","DataHistCanvas",800,400) ;
//    c0->Divide(2,1);
//    c0->cd(1);
//    gPad->SetLeftMargin(0.15) ; hh_genSigData->GetZaxis()->SetTitleOffset(1.4) ; hh_genSigData->Draw("COLZ") ;
//    gPad->SetLeftMargin(0.15) ; hh_searchData->SetMarkerStyle(7) ; hh_searchData->Draw("same") ;
//    c0->cd(2);
//    gPad->SetLeftMargin(0.15) ; hh_genErData->GetZaxis()->SetTitleOffset(1.4) ; hh_genErData->Draw("COLZ") ;
//    gPad->SetLeftMargin(0.15) ; hh_searchData->SetMarkerStyle(7) ; hh_searchData->Draw("same") ;
//    c0->Print(TString::Format("genDataHist_%05iGeV.pdf",thisWimpMass));

    /*
    // Plot data and PDF overlaid
    RooPlot* S1frame = S1.frame(Title("S1 fit projection")) ;
    searchData->plotOn(S1frame,Binning(30)) ; 
    //w->pdf("modelWithConstraints")->plotOn(S1frame);
    w->pdf("bgModel")->plotOn(S1frame);
    w->pdf("sigModel")->plotOn(S1frame,LineStyle(kDashed),LineColor(kRed));

    // Construct histograms with the residuals and pulls of the data w.r.t. the curve
    RooHist* hS1resid = S1frame->residHist() ;    
    RooHist* hS1pull = S1frame->pullHist() ;

    // Create a new frame to draw the residual distribution and add the distribution to the frame
    RooPlot* S1frame2 = S1.frame(Title("Residual Distribution")) ;
    S1frame2->addPlotable(hS1resid,"P") ;

    // Create a new frame to draw the pull distribution and add the distribution to the frame  
    RooPlot* S1frame3 = S1.frame(Title("Pull Distribution")) ;
    S1frame3->addPlotable(hS1pull,"P") ;
    
    TCanvas* c = new TCanvas("rf109_chi2residpull","rf109_chi2residpull",900,300) ;
    c->Divide(3) ;
    c->cd(1) ; gPad->SetLeftMargin(0.15) ; S1frame->GetYaxis()->SetTitleOffset(1.6) ; S1frame->Draw() ;
    c->cd(2) ; gPad->SetLeftMargin(0.15) ; S1frame2->GetYaxis()->SetTitleOffset(1.6) ; S1frame2->Draw() ;
    c->cd(3) ; gPad->SetLeftMargin(0.15) ; S1frame3->GetYaxis()->SetTitleOffset(1.6) ; S1frame3->Draw() ;

    
    RooPlot* S2S1frame = logS2S1.frame(Title("log(S2/S1) fit projection")) ;
    //RooPlot S2S1frame(RooArgSet(S1,logS2S1),Title("log(S2/S1) fit projection")) ;
    searchData->plotOn(S2S1frame,Binning(30)) ;
    w->pdf("bgModel")->plotOn(S2S1frame);
    w->pdf("sigModel")->plotOn(S2S1frame,LineStyle(kDashed),LineColor(kRed));

    // Construct histograms with the residuals and pulls of the data w.r.t. the curve
    RooHist* hS2S1resid = S2S1frame->residHist() ;
    RooHist* hS2S1pull = S2S1frame->pullHist() ;

    // Create a new frame to draw the residual distribution and add the distribution to the frame
    RooPlot* S2S1frame2 = logS2S1.frame(Title("Residual Distribution")) ;
    S2S1frame2->addPlotable(hS2S1resid,"P") ;

    // Create a new frame to draw the pull distribution and add the distribution to the frame  
    RooPlot* S2S1frame3 = logS2S1.frame(Title("Pull Distribution")) ;
    S2S1frame3->addPlotable(hS2S1pull,"P") ;
    
    TCanvas* d = new TCanvas("rf109_chi2residpull","rf109_chi2residpull",900,300) ;
    d->Divide(3) ;
    d->cd(1) ; gPad->SetLeftMargin(0.15) ; S2S1frame->GetYaxis()->SetTitleOffset(1.6) ; S2S1frame->Draw() ;
    d->cd(2) ; gPad->SetLeftMargin(0.15) ; S2S1frame2->GetYaxis()->SetTitleOffset(1.6) ; S2S1frame2->Draw() ;
    d->cd(3) ; gPad->SetLeftMargin(0.15) ; S2S1frame3->GetYaxis()->SetTitleOffset(1.6) ; S2S1frame3->Draw() ;

    
    RooPlot* rframe = r.frame(Title("r fit projection")) ;
    searchData->plotOn(rframe,Binning(30)) ;
    w->pdf("bgModel")->plotOn(rframe);
    w->pdf("sigModel")->plotOn(rframe,LineStyle(kDashed),LineColor(kRed));
    
    RooPlot* zframe = z.frame(Title("z fit projection")) ;
    searchData->plotOn(zframe,Binning(30)) ;
    w->pdf("bgModel")->plotOn(zframe);
    w->pdf("sigModel")->plotOn(zframe,LineStyle(kDashed),LineColor(kRed));
    
    // Draw the frames on the canvas
    TCanvas* c1 = new TCanvas("c1","FitProjCanvas",800,600) ;
    c1->Divide(2,2);
    c1->cd(1);  gPad->SetLeftMargin(0.15) ; S1frame->GetYaxis()->SetTitleOffset(1.4) ; S1frame->Draw() ;
    c1->cd(2);  gPad->SetLeftMargin(0.15) ; S2S1frame->GetYaxis()->SetTitleOffset(1.4) ; S2S1frame->Draw() ;
    
    c1->cd(3);  gPad->SetLeftMargin(0.15) ; rframe->GetYaxis()->SetTitleOffset(1.4) ; rframe->Draw() ;
    c1->cd(4);  gPad->SetLeftMargin(0.15) ; zframe->GetYaxis()->SetTitleOffset(1.4) ; zframe->Draw() ;
    c1->Print(TString::Format("fitProjHist_%05iGeV.pdf",thisWimpMass));
    */

    // Evaluate test statistic in data
    RooAbsReal* nll = w->pdf("modelWithConstraints")->createNLL(*searchData, CloneData(kTRUE), Constrain(*w->set("nuis")) );
    double nLL = nll->getVal();
    cout << "Value of nLL = " << nLL << endl;

    // Draw profile likelihod by hand
    TProfile* pfx = new TProfile("pfx","Profile", 100, 0., 10.);
    for (int i = 0; i< 100; i++ ) {
      double n_i = (double)i/10.;
      w->var("nSig")->setVal(n_i);
      w->var("nSig")->setConstant(kTRUE);

      // Calculate null p-value "by-hand"
      w->pdf("modelWithConstraints")->fitTo(*searchData, Constrain(*w->set("nuis")));
      
      // Evaluate test statistic in data
      RooAbsReal* nll_i = w->pdf("modelWithConstraints")->createNLL(*searchData, CloneData(kTRUE), Constrain(*w->set("nuis")) );
      double nLL_i = nll_i->getVal();
      pfx->Fill(n_i, 2*(nLL_i-nLL));
    }

    TCanvas* c2 = new TCanvas("c2","DataHistCanvas",600,400) ;
    pfx->Draw();

    w->var("nSig")->setVal(0);
    w->var("nSig")->setConstant(kTRUE);

    w->pdf("modelWithConstraints")->fitTo(*searchData, Constrain(*w->set("nuis")));

    RooAbsReal* nll_null = w->pdf("modelWithConstraints")->createNLL(*searchData, CloneData(kTRUE), Constrain(*w->set("nuis")) );
    double nLL_null = nll_null->getVal();
    double q0_obs = 2*(nLL_null-nLL);
    cout << "Value of null nLL = " << nLL_null << endl;
    cout << "Value of q0 = " << q0_obs << endl;

    w->var("nSig")->setConstant(kFALSE);

    // Fit model with constraints
    w->pdf("modelWithConstraints")->fitTo(*searchData, Constrain(*w->set("nuis")));

    ProfileLikelihoodCalculator plc(*searchData, *model);
    // ProfileLikelihoodCalculator makes central intervals
    // to get one-sided upper-limit with desired CL, one needs to use twice larger alpha test size
    if(setOneSided)
      plc.SetTestSize(2.*(1.-limitCL));
    else
      plc.SetConfidenceLevel(limitCL);
    
    LikelihoodInterval* lrint = plc.GetInterval();  // that was easy.

//    cout << "Likelihood ratio = " << lrint->GetLikelihoodRatio()->getVal() << endl;
//
//    // Let's make a plot
//    TCanvas* dataCanvas = new TCanvas("dataCanvas");
//    LikelihoodIntervalPlot plotInt((LikelihoodInterval*)lrint);
//    plotInt.SetTitle("Profile Likelihood Ratio and Posterior for nSig");
//    plotInt.Draw();

    double ul=lrint->UpperLimit(*w->var("nSig"));

    //double ul=lrint->UpperLimit(*w->var("sig_xs"));
//    cout << "Profile likelihood limit: " << ul << ", corresponding to number of signal event = " 
//	 << ul * w->var("T")->getVal() * w->var("xeDensity")->getVal() * w->var("nSigCubicCmPerKgPbDay")->getVal() << endl;

    cout << "Profile likelihood limit: " << ul << ", corresponding to signal cross section = " 
	 << ul / w->var("T")->getVal() / w->var("fidMass")->getVal() / w->var("events_per_pbkgday")->getVal() << endl;

  }

  //////////////////////////////////////////////////////////////////////////////////
  // Try ProfileLikelihoodCalculator (Wilk's theorem)
  //////////////////////////////////////////////////////////////////////////////////
  
  if(plcToyMC) {
    // Make ntuple to store values of UL
    TNtuple* mcTrials = new TNtuple("mcTrials","PLR limits from MC trials","upperLimit:searchEvents");

    // Generate pseudoexperiments for each trial and fit with ProfileLikelihoodCalculator
    for(int i = 0; i<ntrials;i++){
      w->loadSnapshot("initialConditions");
      //generate  one pseudoexperiment
      RooDataSet *searchMcData = w->pdf("bgModel")->generate(*w->set("obs"), RooFit::Extended(),RooFit::Name("searchMcData"));
      //RooDataSet *searchMcData = w->pdf("fullModel")->generate(*w->set("obs"), RooFit::Extended(),RooFit::Name("searchMcData"));
      cout << "*** generated search data, numEntries = " << searchMcData->numEntries() << endl;
      
      ProfileLikelihoodCalculator plc(*searchMcData, *model);
      // ProfileLikelihoodCalculator makes central intervals
      // to get one-sided upper-limit with desired CL, one needs to use twice larger alpha test size
      if(setOneSided)
	plc.SetTestSize(2.*(1.-limitCL));
      else
	plc.SetConfidenceLevel(limitCL);
      
      LikelihoodInterval* lrint = plc.GetInterval();  // that was easy.
      //double ul=lrint->UpperLimit(*w->var("sig_xs"));
      double ul=lrint->UpperLimit(*w->var("nSig"));
      cout << "Profile likelihood limit: " << ul <<endl;
      
      //Fill the tree, tidy
      //later could save three TH2Ds for each trial?
      mcTrials->Fill(ul,searchMcData->numEntries());
      searchMcData->Delete();
    }
    
    //plot and save results
    //mcTrials->Draw("upperLimit>>hist(100,1e-10,1e-08)","","");
    mcTrials->Draw("upperLimit>>histUL(100,0,25)","","");
    mcTrials->Draw("searchEvents>>hist(100,0,100)","","");
    mcTrials->SaveAs( TString::Format("plr_reducedParams_%.0fgev_%0.fd_r%.0fcm.root",
				      w->var("mWimp")->getVal(),w->var("T")->getVal(),w->var("r")->getMax()) );
  }

  //////////////////////////////////////////////////////////////////////////////////
  // Hypothesis test inversion
  //////////////////////////////////////////////////////////////////////////////////

  if(hypoTestInv) {

    // If no input data, generate some
    if(!fname)
      searchData = w->pdf("bgModel")->generate(*w->set("obs"), RooFit::Name("searchData"));

//    // Use Wilk's theorem to get starting place
//    ProfileLikelihoodCalculator plc(*searchData, *model);
//    // ProfileLikelihoodCalculator makes central intervals
//    // to get one-sided upper-limit with desired CL, one needs to use twice larger alpha test size
//    if(setOneSided)
//      plc.SetTestSize(2.*(1.-limitCL));
//    else
//      plc.SetConfidenceLevel(limitCL);
//    
//    LikelihoodInterval* lrint = plc.GetInterval();  // that was easy.
//    //double ul=lrint->UpperLimit(*w->var("sig_xs"));
//    double ul=lrint->UpperLimit(*w->var("nSig"));
//
//    cout << "Profile likelihood limit using Wilk's theorem: " << ul << ", corresponding to signal cross section = " 
//	 << ul / w->var("T")->getVal() / w->var("fidMass")->getVal() / w->var("events_per_pbkgday")->getVal() << endl;
//    
//    // Set range to scan
//    poimin = ul-1.;
//    poimax = ul+5;

    // Set up frequentist calculator    
    FrequentistCalculator *  fc  = new FrequentistCalculator(*searchData, *bModel, *model); // inputs: data, alt model, null model
    fc->SetToys(ntrials,ntrials/2.);    // 1000 for null (S+B) , 500 for alt (B)
    
    // create hypotest inverter 
    HypoTestInverter calc(*fc);
    
    // set confidence level (e.g. 95% upper limits)
    calc.SetConfidenceLevel(limitCL);
    //calc.SetParameter("ReuseAltToys", kTRUE);
   
    // for CLS
    bool useCLs = false;
    calc.UseCLs(useCLs);
    //calc.SetVerbose(false);
    calc.SetVerbose(0);
    
    ToyMCSampler *toymcs = (ToyMCSampler*)calc.GetHypoTestCalculator()->GetTestStatSampler();

    // profile likelihood test statistics 
    ProfileLikelihoodTestStat profll(*model->GetPdf());
    // Use one-sided profile likelihood?
    if(setOneSided)
      profll.SetOneSided(true);
    
    // ratio of profile likelihood - need to pass snapshot for the alt 
    // RatioOfProfiledLikelihoodsTestStat ropl(*sbModel->GetPdf(), *bModel->GetPdf(), bModel->GetSnapshot());
    
    // set the test statistic to use 
    toymcs->SetTestStatistic(&profll);
    
    std::cout << "Doing a fixed scan  in interval : " << poimin << " , " << poimax << std::endl;
    calc.SetFixedScan(npointsHypoTest,poimin,poimax);

    //rott - this suppresses some of the warnings from RooAbsPdf
    cout << "Printing messages:\n";
    RooMsgService::instance().setStreamStatus(0,kFALSE);
    RooMsgService::instance().Print();

//    std::cout << "Doing auto scan" << std::endl;
//    calc.SetAutoScan();
    
    HypoTestInverterResult * result = calc.GetInterval();

    // Write HypoTestInverterResult to workspace
    w->import(*result);
    
    double upperLimit = result->UpperLimit();
    double norm = w->var("T")->getVal() * w->var("fidMass")->getVal() * w->var("events_per_pbkgday")->getVal();
    double ulError = result->UpperLimitEstimatedError();
    cout << "The computed upper limit is: " << upperLimit << " +/- " << ulError 
	 << ", corresponding to signal cross section = " 
	 << upperLimit / norm
	 << endl;

    // compute expected limit
    cout << "Expected upper limits, using the B (alternate) model : " << endl;
    cout << " expected limit (median) " << result->GetExpectedUpperLimit(0) / norm << endl;
    cout << " expected limit (-1 sig) " << result->GetExpectedUpperLimit(-1) / norm << endl;
    cout << " expected limit (+1 sig) " << result->GetExpectedUpperLimit(1) / norm << endl;
    cout << " expected limit (-2 sig) " << result->GetExpectedUpperLimit(-2) / norm << endl;
    cout << " expected limit (+2 sig) " << result->GetExpectedUpperLimit(2) / norm << endl;
    

    ofstream resultsFile;
    resultsFile.open(TString::Format("results_%05iGeV.txt",thisWimpMass));
    resultsFile << thisWimpMass << " " << upperLimit / norm 
		<< " " << result->GetExpectedUpperLimit(-2) / norm 
		<< " " << result->GetExpectedUpperLimit(-1) / norm 
		<< " " << result->GetExpectedUpperLimit(0) / norm 
		<< " " << result->GetExpectedUpperLimit(1) / norm 
		<< " " << result->GetExpectedUpperLimit(2) / norm 
		<< endl;
    
    resultsFile.close();

    // plot now the result of the scan 
    HypoTestInverterPlot *plot = new HypoTestInverterPlot("HTI_Result_Plot","",result);
    
    // plot in a new canvas with style
    TCanvas * c1 = new TCanvas("HypoTestInverter Scan"); 
    c1->SetLogy(false);
    
    //plot->Draw("CLb 2CL");  // plot also CLb and CLs+b 
    plot->Draw("CLb");  // plot also CLb and CLs+b 
    //plot->Draw("OBS");  // plot only observed p-value
    c1->Print(TString::Format("hypoTest_%05iGeV.pdf",thisWimpMass));
    
    // plot test statistics distributions for the two hypothesis
    // when distribution is generated (case of FrequentistCalculators)
    const int n = result->ArraySize();
    if (n> 0 &&  result->GetResult(0)->GetNullDistribution() ) { 
      TCanvas * c2 = new TCanvas("Test Statistic Distributions","",2);
      if (n > 1) {
         int ny = TMath::CeilNint( sqrt(n) );
         int nx = TMath::CeilNint(double(n)/ny);
         c2->Divide( nx,ny);
      }
      for (int i=0; i<n; i++) {
         if (n > 1) c2->cd(i+1);
         SamplingDistPlot * pl = plot->MakeTestStatPlot(i);
         pl->SetLogYaxis(true);
         pl->Draw();
	 
         // Write sampling distributions to workspace
         SamplingDistribution * nullSampDist = (SamplingDistribution*)result->GetNullTestStatDist(i);
         w->import(*nullSampDist);
         SamplingDistribution * altSampDist = (SamplingDistribution*)result->GetAltTestStatDist(i);
         w->import(*altSampDist);
      }
      c2->Print(TString::Format("hypoTestDists_%05iGeV.pdf",thisWimpMass));
    }

    // Write workspace to file
    w->writeToFile(TString::Format("ws_model_%05iGeV.root",thisWimpMass));
  }
  
  //////////////////////////////////////////////////////////////////////////////////
  // Null hypothesis test (i.e. "discovery")
  //////////////////////////////////////////////////////////////////////////////////

  if(hypoTestNull) {

    // If no input data, generate some
    if(!fname)
      searchData = bModel->GetPdf()->generate(*w->set("obs"), RooFit::Name("searchData"));

    // Fit with nSig floating
    w->pdf("modelWithConstraints")->fitTo(*searchData, Constrain(*w->set("nuis")));

    // Capture fit results with nSig floating
    w->saveSnapshot("genConditions",w->allVars());

    // Evaluate denominator of test statistic in data
    RooAbsReal* nll = w->pdf("modelWithConstraints")->createNLL(*searchData, CloneData(kTRUE), Constrain(*w->set("nuis")) );
    double nLL = nll->getVal();
    cout << "Value of nLL = " << nLL << endl;

    // Set nSig to zero and get NLL
    w->var("nSig")->setVal(0);
    w->var("nSig")->setConstant(kTRUE);

    // Refit with nSig=0
    w->pdf("modelWithConstraints")->fitTo(*searchData, Constrain(*w->set("nuis")));

    // Capture fit results with nSig=0
    w->saveSnapshot("genNullConditions",w->allVars());

    // Evaluate numerator of test statistic
    RooAbsReal* nll_null = w->pdf("modelWithConstraints")->createNLL(*searchData, CloneData(kTRUE), Constrain(*w->set("nuis")) );
    double nLL_null = nll_null->getVal();

    // Determine null test stat in data
    double q0_obs = 2*(nLL_null-nLL);
    cout << "Value of null nLL = " << nLL_null << endl;
    cout << "Value of q0 = " << q0_obs << endl;

    w->var("nSig")->setConstant(kFALSE);

    TNtuple* mcTrials = new TNtuple("mcTrials","PLR limits from MC trials","testStat:searchEvents");
    double q0_i = 0;
    double count = 0;

    // Generate pseudoexperiments for each trial and fit with ProfileLikelihoodCalculator
    for(int i = 0; i<ntrials;i++){
      w->loadSnapshot("genConditions");
      //w->loadSnapshot("genNullConditions");

      w->var("nSig")->setVal(0);
      w->var("nSig")->setConstant(kTRUE);

      //generate  one pseudoexperiment
      RooDataSet *searchMcData = w->pdf("modelWithConstraints")->generate(*w->set("obs"), RooFit::Extended(),RooFit::Name("searchMcData"));
      cout << "*** generated search data, numEntries = " << searchMcData->numEntries() << endl;

      // Evaluate numerator of test statistic in pseudo-experiment
      w->pdf("modelWithConstraints")->fitTo(*searchMcData, Constrain(*w->set("nuis")));
      RooAbsReal* nll_null_i = w->pdf("modelWithConstraints")->createNLL(*searchMcData, CloneData(kTRUE), Constrain(*w->set("nuis")) );
      nLL_null = nll_null_i->getVal();
 
      // Let nSig float
      w->loadSnapshot("genConditions");
      w->var("nSig")->setConstant(kFALSE);

      // Evaluate denominator of test statistic in pseudo-experiment
      w->pdf("modelWithConstraints")->fitTo(*searchMcData, Constrain(*w->set("nuis")));
      RooAbsReal* nll_i = w->pdf("modelWithConstraints")->createNLL(*searchMcData, CloneData(kTRUE), Constrain(*w->set("nuis")) );
      nLL = nll_i->getVal();
      q0_i = 2*(nLL_null-nLL);
     
      if(q0_i > q0_obs) count++;
      
      //Fill the tree, tidy
      mcTrials->Fill(q0_i,searchMcData->numEntries());
      searchMcData->Delete();
    }

    cout << "Null p-value = " << count/ntrials << endl;
    
    //plot and save results
    mcTrials->Draw("testStat>>histTS(100,0,5)","","");
    //mcTrials->Draw("searchEvents>>histSE(100,0,200)","","");
    mcTrials->SaveAs( TString::Format("plr_reducedParams_%.0fgev_%0.fd_r%.0fcm.root",
				      w->var("mWimp")->getVal(),w->var("T")->getVal(),w->var("r")->getMax()) );
    /*
    // Test bg model as null    
    FrequentistCalculator *  fc  = new FrequentistCalculator(*searchData, *model, *bModel); // inputs: data, alt model, null model
    fc->SetToys(ntrials,ntrials); 

    // create the test statistics
    ProfileLikelihoodTestStat profll(*model->GetPdf());
    // use one-sided profile likelihood
    profll.SetOneSidedDiscovery(true);

    // configure  ToyMCSampler and set the test statistics
    ToyMCSampler *toymcs = (ToyMCSampler*)fc->GetTestStatSampler();
    toymcs->SetTestStatistic(&profll);

    // run the test
    HypoTestResult* result = fc->GetHypoTest();
    result->Print();

    // plot test statistic distributions
    HypoTestPlot * plot = new HypoTestPlot(*result);
    plot->Draw();
    
    // Write workspace to file
    w->writeToFile(TString::Format("ws_model_nullTest_%05iGeV.root",thisWimpMass));
    */    

  }

  // print timing info
  t.Stop();
  t.Print();
}


#ifndef __CINT__

int main(int argc, char *argv[]) {
  if(argc != 3) {
    cout << "Usage: runMass WIMPfile mass\nWhere mass is an int and WIMPfile is a root file\n";
    return 0;
  }

  PLRun3(argv[1],atoi(argv[2]),kFALSE,kFALSE,kTRUE);

  return 0;

}

#endif
