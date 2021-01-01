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
#include "RooDataSet.h"
#include "RooEfficiency.h"
#include "RooExtendPdf.h"
#include "RooFFTConvPdf.h"
#include "RooGaussian.h"
#include "RooProjectedPdf.h"
#include "RooPlot.h"
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

void buildModel(RooWorkspace* w)
{
  
  ////////////////////////////////////////////////////////////////////////////////////////
  // The Model building stage
  ////////////////////////////////////////////////////////////////////////////////////////

  //S1.setBins("fft",10000); // Sampling frequency for FFT

  // Import parameters
  w->import(sig_xs);
  w->import(mWimp);
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
  w->import(S1);
  w->import(logS2S1);
  
  w->import(rateCompton);
  w->import(ra);
  w->import(za);
  w->import(zb);
  
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
  w->factory("cexpr::Eee('aEee+bEee*S1+cEee*pow(S1,2)+dEee*pow(S1,3)+eEee*pow(S1,4)+fEee*pow(S1,5)',S1,aEee[0.465028],bEee[0.329153],cEee[-0.0180685],dEee[0.000993007],eEee[-2.8398e-05],fEee[3.16677e-07])"); // [keVee]
  w->factory("cexpr::Enr('aEnr+bEnr*S1+cEnr*pow(S1,2)+dEnr*pow(S1,3)+eEnr*pow(S1,4)+fEnr*pow(S1,5)',S1,aEnr[1.08087],bEnr[1.53183],cEnr[-0.0874158],dEnr[0.00527215],eEnr[-0.000163157],fEnr[1.93408e-06])"); // [keVnr]
  //w->import(Enr);

  // Terms to transform integral from Int[f(E)dE] to Int[g(E(S1)) dE/dS1 dS1]
  w->factory("cexpr::dEee_dS1('bEee+cEee*S1+dEee*pow(S1,2)+eEee*pow(S1,3)+fEee*pow(S1,4)',S1,bEee,cEee,dEee,eEee,fEee)");
  w->factory("cexpr::dEnr_dS1('bEnr+cEnr*S1+dEnr*pow(S1,2)+eEnr*pow(S1,3)+fEnr*pow(S1,4)',S1,bEnr,cEnr,dEnr,eEnr,fEnr)"); 
  //w->factory("dEnr_dS1[1.]");

  if(!bandsFixed) {
    //The numerical values are from fits to NEST sims after applying 14%LC, 30%SPEres, 60%CC, 20%SEres.
    //build NR functions for median log10(S2/S1) (aka logS2S1) as function of S1. 
    //build NR functions for gaussian sigma in log10(S2/S1) as function of S1 (via parametrisation of mu+1sigma)
    w->factory("cexpr::meanNr('nrmA*pow(S1,nrmB)+nrmC*S1',nrmA[1.12,1.0,2.0],S1,nrmB[-0.124,-0.2,0],nrmC[-0.005,-0.01,0])");
    w->factory("cexpr::sigmaNr('nrsA*pow(S1,nrsB)-nrmA*pow(S1,nrmB)+(nrsC-nrmC)*S1',nrmA, nrmB, nrmC,nrsA[1.47,1.0,3.0],S1,nrsB[-0.195,-0.2,0],nrsC[-0.0034,-0.01,0])");
    
    //build ER functions for median log10(S2/S1) (aka logS2S1) as function of S1. 
    //build ER functions for gaussian sigma in log10(S2/S1) as function of S1 (via parametrisation of mu+1sigma)
    w->factory("cexpr::meanEr('ermA*pow(S1,ermB)+ermC*S1',ermA[1.71,1.0,2.0],S1,ermB[-0.161,-0.2,0],ermC[-0.0056,-0.01,0])");
    w->factory("cexpr::sigmaEr('ersA*pow(S1,ersB)-ermA*pow(S1,ermB)+(ersC-ermC)*S1',ermA, ermB, ermC,ersA[1.95,1.0,3.0],S1,ersB[-0.175,-0.2,0],ersC[-0.0032,-0.01,0])");
  }
  else {
    // try a fixed band
    w->factory("cexpr::meanNr('nrmA*pow(S1,nrmB)+nrmC*S1',nrmA[1.12],S1,nrmB[-0.124],nrmC[-0.005])");
    w->factory("cexpr::sigmaNr('nrsA*pow(S1,nrsB)-nrmA*pow(S1,nrmB)+(nrsC-nrmC)*S1',nrmA, nrmB, nrmC,nrsA[1.47],S1,nrsB[-0.195],nrsC[-0.0034])");
    w->factory("cexpr::meanEr('ermA*pow(S1,ermB)+ermC*S1',ermA[1.71],S1,ermB[-0.161],ermC[-0.0056])");
    w->factory("cexpr::sigmaEr('ersA*pow(S1,ersB)-ermA*pow(S1,ermB)+(ersC-ermC)*S1',ermA, ermB, ermC,ersA[1.95],S1,ersB[-0.175],ersC[-0.0032])"); 
  }

  // Create pdfs P(log10(S2/S1)|S1)
  w->factory("Gaussian::bandNr(logS2S1,meanNr,sigmaNr)");
  w->factory("Gaussian::bandEr(logS2S1,meanEr,sigmaEr)");

  // Generate with shift in the mean of ER
  w->factory("cexpr::meanErSys('ermASys*pow(S1,ermBSys)+ermCSys*S1',ermASys[1.54,1.0,2.0],S1,ermBSys[-0.161,-0.2,0],ermCSys[-0.0050,-0.01,0])");
  w->factory("cexpr::sigmaErSys('ersASys*pow(S1,ersBSys)-ermASys*pow(S1,ermBSys)+(ersCSys-ermCSys)*S1',ermASys, ermBSys, ermCSys,ersASys[1.95,1.0,3.0],S1,ersBSys[-0.175,-0.2,0],ersCSys[-0.1732,-0.01,0])");

  // Generate with best fit to data
  w->factory("cexpr::meanNrData('nrmAData*pow(S1,nrmBData)+nrmCData*S1',nrmAData[1.85,1.0,2.0],S1,nrmBData[-0.099,-0.2,0],nrmCData[-0.0025,-0.01,0])");
  w->factory("cexpr::sigmaNrData('nrsAData*pow(S1,nrsBData)-nrmAData*pow(S1,nrmBData)+(nrsCData-nrmCData)*S1',nrmAData, nrmBData, nrmCData,nrsAData[2.1,1.0,3.0],S1,nrsBData[-0.108,-0.2,0],nrsCData[-1.27652e-07,-0.01,0])");
  w->factory("cexpr::meanErData('ermAData*pow(S1,ermBData)+ermCData*S1',ermAData[2.22,1.0,3.0],S1,ermBData[-0.0575,-0.2,0],ermCData[-0.0010,-0.01,0])");
  w->factory("cexpr::sigmaErData('ersAData*pow(S1,ersBData)-ermAData*pow(S1,ermBData)+(ersCData-ermCData)*S1',ermAData, ermBData, ermCData,ersAData[2.06,1.0,3.0],S1,ersBData[-0.0755,-0.2,0],ersCData[-0.0056,-0.01,0])");

  // Create pdfs P(log10(S2/S1)|S1)
  //w->factory("Gaussian::bandNrGen(logS2S1,meanNrData,sigmaNrData)");
  //w->factory("Gaussian::bandErGen(logS2S1,meanErData,sigmaErData)");

  //w->factory("Gaussian::bandNrGen(logS2S1,meanNr,sigmaNr)");
  //w->factory("Gaussian::bandErGen(logS2S1,meanErSys,sigmaErSys)");

  w->factory("Gaussian::bandNrGen(logS2S1,meanNr,sigmaNr)");
  w->factory("Gaussian::bandErGen(logS2S1,meanEr,sigmaEr)");

  // Define some useful sets for later  
  w->defineSet("obs","S1,logS2S1,r,z");
  w->defineSet("nuis","nCompton,nActivatedXe,nRn,nrmA,nrmB,nrmC,nrsA,nrsB,nrsC,ermA,ermB,ermC,ersA,ersB,ersC");

  ////////////////////////////////////////////////////////////////////////////////////////
  // Efficiency model
  ////////////////////////////////////////////////////////////////////////////////////////

  //S1 efficiency
  w->factory("CEXPR::effS1('1./( 1.+exp(-2.*kS1*(S1-thresS1)) )',S1,kS1[1.77374],thresS1[3.48613])");
  //w->factory("effS1[1]");

  // S2 efficiency - need to figure out how to incorporate this
  w->factory("cexpr::S2('pow(10,logS2S1)*S1*exp(-(zLiquid-z)/(driftVelocity*eTau))',S1,logS2S1,z,zLiquid,driftVelocity,eTau)");
  //w->factory("CEXPR::effS2('1./( 1.+exp(-2.*kS2*(S2-thresNe*singleElec)) )',S2,kS2[20.],thresNe,singleElec)");
  w->factory("CEXPR::effS2('S2>thresNe*singleElec',S2,thresNe,singleElec)");

  w->factory("CEXPR::bandNrWithS2Eff('bandNr*effS2',bandNr,effS2)");
  w->factory("CEXPR::bandErWithS2Eff('bandEr*effS2',bandEr,effS2)");

  ////////////////////////////////////////////////////////////////////////////////////////
  // Build up model with terms
  // P(log10(S2/S1)|S1)P(S1)P(r,z)
  // for both signal & bg
  ////////////////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////////////////
  // Signal model
  ////////////////////////////////////////////////////////////////////////////////////////
  
  // Declare pdf of WIMP energy spectrum P(Enr)
  //gROOT->ProcessLine(".L RooWimpSpec_ctsPerkeVkgdaypb.cxx+");//load a custom library
  //w->importClassCode("RooWimpSpec_ctsPerkeVkgdaypb",kTRUE);//import custom class into workspace
  //w->factory("WimpSpec_ctsPerkeVkgdaypb::dRdE(Enr,mWimp,targetA,targetZ,rho,dayOfYear,vEsc,v0,ff_c,ff_a,ff_s,idm_delta)");
   
  // Declare pdf of WIMP energy spectrum P(S1)
  gROOT->ProcessLine(".L RooWimpSpectrum.cxx+");//load a custom library
  w->importClassCode("RooWimpSpectrum", kTRUE);//import custom class into workspace
  w->factory("WimpSpectrum::dRdE(Enr,mWimp, v0, vE, vEsc, rho, ff_c, ff_r0, ff_a, ff_s, idm_delta, halo_beta)");
  w->factory("CEXPR::dRdS1('dRdE*dEnr_dS1*effS1',dRdE,dEnr_dS1,effS1)");

  //w->factory("CEXPR::dRdS1('exp(-S1/a)',S1,a[14.7])");

  //RooAbsReal* wimpNorm = w->pdf("dRdE")->createIntegral(*w->var("Enr"));
  RooAbsReal* wimpNorm = w->pdf("dRdS1")->createIntegral(*w->var("S1"));
  wimpNorm->SetName("wimpNorm");
  w->import(*wimpNorm);

//  // Check normalization, Enr
//  Enr.setRange("r0",4.5,25.);
//  cout << "Number in range [4.5,25] = " << (w->pdf("dRdE")->createIntegral(*w->var("Enr"),Range("r0"))->getVal() 
//					    * w->var("T")->getVal() *100. *w->var("sig_xs")->getVal()) << endl;
//
//  Enr.setRange("r1",5.,6.);
//  cout << "Number in range [5,6] = " << (w->pdf("dRdE")->createIntegral(*w->var("Enr"),Range("r1"))->getVal() 
//					    * w->var("T")->getVal() *100. *w->var("sig_xs")->getVal()) << endl;
//
//  Enr.setRange("r2",10.,11.);
//  cout << "Number in range [10,11] = " << (w->pdf("dRdE")->createIntegral(*w->var("Enr"),Range("r2"))->getVal() 
//					    * w->var("T")->getVal() *100. *w->var("sig_xs")->getVal()) << endl;
//
//  Enr.setRange("r3",24.,25.);
//  cout << "Number in range [24,25] = " << (w->pdf("dRdE")->createIntegral(*w->var("Enr"),Range("r3"))->getVal() 
//					    * w->var("T")->getVal() *100. *w->var("sig_xs")->getVal()) << endl;

  //Check normalization, S1
  S1.setRange("r0",2.,30.);
  cout << "Number in range [4.5,25] = " << (w->pdf("dRdS1")->createIntegral(*w->var("S1"),Range("r0"))->getVal() 
					    * w->var("T")->getVal() *100. *w->var("sig_xs")->getVal()) << endl;

  S1.setRange("r1",3.2,4.);
  cout << "Number in range [5,6] = " << (w->pdf("dRdS1")->createIntegral(*w->var("S1"),Range("r1"))->getVal() 
					    * w->var("T")->getVal() *100. *w->var("sig_xs")->getVal()) << endl;

  S1.setRange("r2",8.5,9.5);
  cout << "Number in range [10,11] = " << (w->pdf("dRdS1")->createIntegral(*w->var("S1"),Range("r2"))->getVal() 
					    * w->var("T")->getVal() *100. *w->var("sig_xs")->getVal()) << endl;

  S1.setRange("r3",28.7,30.);
  cout << "Number in range [24,25] = " << (w->pdf("dRdS1")->createIntegral(*w->var("S1"),Range("r3"))->getVal() 
					    * w->var("T")->getVal() *100. *w->var("sig_xs")->getVal()) << endl;

  w->factory("prod::wimpSpec1(wimpNorm,sig_xs,M[100.],T,dRdE)");
  w->factory("prod::wimpSpec2(wimpNorm,sig_xs,M[100.],T,dRdS1)");

  //RooPlot* frame = Enr.frame(Title("Wimp Spectrum")) ;
  ////w->function("wimpSpec1")->plotOn(frame);
  //w->pdf("dRdE")->plotOn(frame,LineColor(kBlue));

  RooPlot* frame2 = S1.frame(Title("Wimp Spectrum")) ;
  //w->function("wimpSpec2")->plotOn(frame2);
  w->pdf("dRdS1")->plotOn(frame2,LineColor(kBlue));

  //TH1* hh_wimpSpec1 = w->function("wimpSpec1")->createHistogram("Enr",20) ;
  TH1* hh_wimpSpec2 = w->function("wimpSpec2")->createHistogram("S1",27) ;

  TCanvas* c00 = new TCanvas("c00","Wimp Spectrum",600,400) ;
  //c00->Divide(2,1);
  //c00->cd(1);
  //frame->Draw() ;
  //c00->cd(2);
  frame2->Draw() ;

  TCanvas* c01 = new TCanvas("c01","Wimp Spectrum Hists",600,400) ;
  //hh_wimpSpec1->Draw() ;
  hh_wimpSpec2->Draw() ;


  // Add energy resolution - this should probably be incorporated directly into dR/dE
  //w->factory("prod::sigmaEnrRes(relEnrRes[0.254865],Enr)");
  //w->factory("Gaussian::resEnr(Enr,meanEnrRes[0],sigmaEnrRes)");
  //w->factory("FFTConvPdf::dRdS1ConvGaus(S1,dRdS1,resEnr)");

  // Second try- smear in S1 directly
  //w->factory("cexpr::sigmaS1Res('0.14*sqrt(S1)',S1)");
  //w->factory("Gaussian::resS1(S1,meanS1Res[0],sigmaS1Res)");
  //w->factory("FFTConvPdf::dRdS1Gaus(S1,dRdS1,resS1)");
  //w->factory("PROD::dRdS1ConvGaus(dRdS1Gaus,effS1)");

  // Declare P(r) and P(z) for WIMPs  - only need to do this once for all other uniform pdfs
  w->factory("CEXPR::UniformSpatial('TMath::TwoPi()*r',r,z)");
  
  // Multiply P(log10(S2/S1)|S1)P(S1)P(r)P(z)
  //w->factory("PROD::sigModel(bandNr|S1, dRdS1, UniformSpatial)");
  w->factory("PROD::sigModelUnNorm(bandNr|S1, dRdS1, UniformSpatial)");

  w->factory("CEXPR::sigModelTest('bandNr*dRdS1',bandNr,dRdS1)");
  w->factory("PROD::sigModelGen(bandNrGen|S1, dRdS1, UniformSpatial)");
  w->factory("PROD::sigEModel(dRdS1, UniformSpatial)");

  double sigNormInt = w->pdf("sigModelUnNorm")->createIntegral( *w->set("obs") )->getVal(); // [day-1 pb-1 cm3]

  double bandNrNormInt = w->pdf("bandNr")->createIntegral( RooArgSet(*w->var("S1"),*w->var("logS2S1")) )->getVal(); // [day-1 pb-1 cm3]
  double dRdS1NormInt = w->pdf("dRdS1")->createIntegral( RooArgSet(*w->var("S1")) )->getVal(); // [day-1 pb-1 cm3]
  double uniformSpatialNormInt = w->pdf("UniformSpatial")->createIntegral( RooArgSet(*w->var("r"), *w->var("z")) )->getVal(); // [day-1 pb-1 cm3]

  double sigTestNormInt = w->pdf("sigModelTest")->createIntegral( RooArgSet(*w->var("S1"),*w->var("logS2S1")) )->getVal(); // [day-1 pb-1 cm3]

  cout << "sigNormInt = " << sigNormInt << endl;
  cout << "bandNrNormInt = " << bandNrNormInt << endl;
  cout << "dRdS1NormInt = " << dRdS1NormInt << endl;
  cout << "uniformSpatialNormInt = " << uniformSpatialNormInt << endl;
  cout << "********************" << endl;
  cout << "sigTestNormInt = " << sigTestNormInt << endl;
  cout << "********************" << endl;

  // Make a variable corresponding to the number of signal events
  ////RooAbsReal* nSigCubicCmPerKgPbDay =w->pdf("sigModel")->createIntegral( *w->set("obs") ); // [day-1 pb-1 cm3]
  //RooAbsReal* nSigCubicCmPerKgPbDay =w->pdf("sigEModel")->createIntegral( RooArgSet(*w->var("S1"),*w->var("r"),*w->var("z")) ); // [day-1 pb-1 cm3]
  //nSigCubicCmPerKgPbDay->SetName("nSigCubicCmPerKgPbDay");
  //w->import(*nSigCubicCmPerKgPbDay);
  
  double nSigCubicCmPerKgPbDay =w->pdf("sigEModel")->createIntegral( RooArgSet(*w->var("S1"),*w->var("r"),*w->var("z")) )->getVal(); // [day-1 pb-1 cm3]
  cout << "double nSigCubicCmPerKgPbDay = " << nSigCubicCmPerKgPbDay << endl;

  //w->factory("nSigCubicCmPerKgPbDay[1]");
  //w->var("nSigCubicCmPerKgPbDay")->setVal(nSigCubicCmPerKgPbDay);
  //w->factory("prod::nSig(sig_xs,T, xeDensity, nSigCubicCmPerKgPbDay)"); // number of signal events 

  //w->factory("nSig[1, 0,100]"); // number of signal events 
  w->factory("nSig[2, 0,100]"); // number of signal events 

  w->factory("ExtendPdf::sigModel(sigModelUnNorm, nSig)");
  
  ////////////////////////////////////////////////////////////////////////////////////////
  // Background model
  ////////////////////////////////////////////////////////////////////////////////////////
  
  // ER background from detector materials
  w->factory("CEXPR::comptonSpectrum('dEee_dS1*effS1',dEee_dS1,effS1)"); // P(S1)
  w->factory("CEXPR::comptonSpatial('TMath::TwoPi()*r*rateCompton*exp(pow(r/ra, 2.)+pow((z-za)/zb, 2.))', r,rateCompton,ra,z ,za,zb)"); // P(r,z)
  w->factory("PROD::comptonModelUnNorm(bandEr|S1, comptonSpectrum, comptonSpatial)"); // P(log10(S2/S1))P(S1)P(r,z)
  w->factory("PROD::comptonModelGen(bandErGen|S1, comptonSpectrum, comptonSpatial)"); // P(log10(S2/S1))P(S1)P(r,z)
  //w->factory("PROD::comptonModel(bandErWithS2Eff|S1, comptonSpectrum, comptonSpatial)"); // P(log10(S2/S1))P(S1)P(r,z)
  w->factory("PROD::comptonEModel(comptonSpectrum, comptonSpatial)"); // P(log10(S2/S1))P(S1)P(r,z)
  
  // Calculate number of expected Compton events due to materials
  //RooAbsReal* nComptonCubicCmPerKgDay = w->pdf("comptonModel")->createIntegral(*w->set("obs")); // [day-1]
  RooAbsReal* nComptonCubicCmPerKgDay = w->pdf("comptonEModel")->createIntegral(RooArgSet(*w->var("S1"),*w->var("r"),*w->var("z"))); // [day-1]
  nComptonCubicCmPerKgDay->SetName("nComptonCubicCmPerKgDay");
  w->import(*nComptonCubicCmPerKgDay);
  w->factory("prod::nComptonPred(T, xeDensity, nComptonCubicCmPerKgDay)"); // number of Compton events predicted by bg model
  //w->factory("nComptonPred[2,0,100]"); // number of Compton events predicted by bg model
  w->factory("prod::sigmaCompton(nComptonPred,relUncertaintyBg)"); // error on number of Compton events predicted by bg model

  w->factory("ExtendPdf::comptonModel(comptonModelUnNorm, nCompton)");
  
  // ER background from Xe-127
  w->factory("CEXPR::activatedXeSpectrum('dEee_dS1*effS1',dEee_dS1,effS1)"); // P(S1)
  //w->factory("CEXPR::activatedXeSpatial('TMath::TwoPi()*r*rateActivatedXe',r,rateActivatedXe,z)"); // P(r,z)
  w->factory("CEXPR::activatedXeSpatial('TMath::TwoPi()*r*(rateActivatedXeA*exp(-(T/2.)/52.5)*exp(fabs(r/ra)))+rateActivatedXeB*exp(-0.5)*exp(fabs(z-za)/zb)',r,rateActivatedXeA,T,ra,rateActivatedXeB,z,za,zb)"); // P(r,z)
  w->factory("PROD::activatedXeModelUnNorm(bandEr|S1, activatedXeSpectrum, activatedXeSpatial)"); // P(log10(S2/S1))P(S1)P(r,z)
  w->factory("PROD::activatedXeModelGen(bandErGen|S1, activatedXeSpectrum, activatedXeSpatial)"); // P(log10(S2/S1))P(S1)P(r,z)
  //w->factory("PROD::activatedXeModel(bandErWithS2Eff|S1, activatedXeSpectrum, activatedXeSpatial)"); // P(log10(S2/S1))P(S1)P(r,z)
  w->factory("PROD::activatedXeEModel(activatedXeSpectrum, activatedXeSpatial)"); // P(log10(S2/S1))P(S1)P(r,z)
  
  // Calculate number of expected Xe-127 events
  //RooAbsReal* nActivatedXeCubicCmPerKgDay = w->pdf("activatedXeModel")->createIntegral(*w->set("obs")); // [day-1]
  RooAbsReal* nActivatedXeCubicCmPerKgDay = w->pdf("activatedXeEModel")->createIntegral(RooArgSet(*w->var("S1"),*w->var("r"),*w->var("z"))); // [day-1]
  nActivatedXeCubicCmPerKgDay->SetName("nActivatedXeCubicCmPerKgDay");
  w->import(*nActivatedXeCubicCmPerKgDay);
  w->factory("prod::nActivatedXePred(T, xeDensity, nActivatedXeCubicCmPerKgDay)"); // number of Xe-127 events predicted by bg model
  //w->factory("nActivatedXePred[1,0,100]"); // number of Xe-127 events predicted by bg model
  w->factory("prod::sigmaActivatedXe(nActivatedXePred,relUncertaintyBg)"); // error on number of Xe-127 events predicted by bg model
  
  w->factory("ExtendPdf::activatedXeModel(activatedXeModelUnNorm, nActivatedXe)");

  // ER background from Rn-222 and Kr-85
  w->factory("CEXPR::rnSpectrum('dEee_dS1*effS1',dEee_dS1,effS1)"); // P(S1)
  w->factory("CEXPR::rnSpatial('TMath::TwoPi()*r*rateRn',r,rateRn,z)"); // P(r,z)
  w->factory("PROD::rnModelUnNorm(bandEr|S1, rnSpectrum, rnSpatial)"); // P(log10(S2/S1))P(S1)P(r,z)
  w->factory("PROD::rnModelGen(bandErGen|S1, rnSpectrum, rnSpatial)"); // P(log10(S2/S1))P(S1)P(r,z)
  //w->factory("PROD::rnModel(bandErWithS2Eff|S1, rnSpectrum, rnSpatial)"); // P(log10(S2/S1))P(S1)P(r,z)
  w->factory("PROD::rnEModel(rnSpectrum, rnSpatial)"); // P(log10(S2/S1))P(S1)P(r,z)
  
  // Calculate number of expected Rn-222 events
  //RooAbsReal* nRnCubicCmPerKgDay = w->pdf("rnModel")->createIntegral(*w->set("obs")); // [day-1]
  RooAbsReal* nRnCubicCmPerKgDay = w->pdf("rnEModel")->createIntegral(RooArgSet(*w->var("S1"),*w->var("r"),*w->var("z"))); // [day-1]
  nRnCubicCmPerKgDay->SetName("nRnCubicCmPerKgDay");
  w->import(*nRnCubicCmPerKgDay);
  w->factory("prod::nRnPred(T, xeDensity, nRnCubicCmPerKgDay)"); // number of Xe-127 events predicted by bg model
  //w->factory("nRnPred[1,0,100]"); // number of Xe-127 events predicted by bg model
  w->factory("prod::sigmaRn(nRnPred,relUncertaintyBg)"); // error on number of Xe-127 events predicted by bg model

  w->factory("ExtendPdf::rnModel(rnModelUnNorm, nRn)");
  
  // Set the initial number of backgrounds to the prediction
  double nComptonInit =  w->function("nComptonPred")->getVal();
  double nActivatedXeInit = w->function("nActivatedXePred")->getVal();
  double nRnInit = w->function("nRnPred")->getVal();

  w->var("nCompton")->setVal(nComptonInit);
  w->var("nActivatedXe")->setVal(nActivatedXeInit);
  w->var("nRn")->setVal(nRnInit);

  //w->var("nCompton")->setVal(0);
  //w->var("nActivatedXe")->setVal(0);
  //w->var("nRn")->setVal(0);
  
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
  //w->factory("SUM::fullModel(nSig*sigModel,nCompton*comptonModel,nActivatedXe*activatedXeModel,nRn*rnModel)");
  //w->factory("SUM::bgModel(nCompton*comptonModel,nActivatedXe*activatedXeModel,nRn*rnModel)");
  //w->factory("ExtendPdf::fullModel(sigModelUnNorm,nSig)");
  //w->factory("ExtendPdf::bgModel(comptonModelUnNorm,nCompton)");

  //w->factory("SUM::fullModel(sigModel,comptonModel,activatedXeModel,rnModel)");
  //w->factory("SUM::fullModel(sigModel,comptonModel)");
  //w->factory("SUM::bgModel(comptonModel,activatedXeModel,rnModel)");

  RooAddPdf fullModel("fullModel","fullModel",RooArgList(*w->pdf("sigModel"),*w->pdf("comptonModel"),*w->pdf("activatedXeModel"),*w->pdf("rnModel")) );
  w->import(fullModel);
  RooAddPdf bgModel("bgModel","bgModel",RooArgList(*w->pdf("comptonModel"),*w->pdf("activatedXeModel"),*w->pdf("rnModel")) );
  w->import(bgModel);

  w->factory("SUM::bgModelGen(nCompton*comptonModelGen,nActivatedXe*activatedXeModelGen,nRn*rnModelGen)");
  
  // Add constraints on nuisance parameters
  w->factory("Gaussian::comptonConstraint(nComptonPred, nCompton, sigmaCompton)");
  w->factory("Gaussian::activatedXeConstraint(nActivatedXePred, nActivatedXe, sigmaActivatedXe)");
  w->factory("Gaussian::rnConstraint(nRnPred, nRn, sigmaRn)");
  
  w->factory("PROD::modelWithConstraints(fullModel,comptonConstraint,activatedXeConstraint,rnConstraint)"); // product of terms
  w->importClassCode();
  
  
  ////////////////////////////////////////////////////////////////////////////////////////
  // Make simultaneous pdf
  ////////////////////////////////////////////////////////////////////////////////////////
  /*
  // Decalre extended pdfs for the calibration samples
  w->factory("CEXPR::fallingNrSpectrum('exp(-S1/14.7)*effS1*dEnr_dS1',S1,dEnr_dS1,effS1)"); // P(S1)
  w->factory("PROD::fallingNrModel(bandNr|S1,fallingNrSpectrum)"); // P(log10(S2/S1))P(S1)
  w->factory("PROD::fallingNrModelGen(bandNrGen|S1,fallingNrSpectrum)"); // P(log10(S2/S1))P(S1)
  w->factory("ExtendPdf::calibNrModel(fallingNrModel,nNrCalib)");  
  w->factory("ExtendPdf::calibNrModelGen(fallingNrModelGen,nNrCalib)");  

  w->factory("PROD::flatErModel(bandEr|S1,comptonSpectrum)"); // P(log10(S2/S1))P(S1)
  w->factory("PROD::flatErModelGen(bandErGen|S1,comptonSpectrum)"); // P(log10(S2/S1))P(S1)
  //w->factory("PROD::flatErModel(bandErWithS2Eff|S1,comptonSpectrum)"); // P(log10(S2/S1))P(S1)
  w->factory("ExtendPdf::calibErModel(flatErModel,nErCalib)");
  w->factory("ExtendPdf::calibErModelGen(flatErModelGen,nErCalib)");
  
  // Decalre RooCategory    
  w->factory("sample[nr,er,search]");
  
  //make a RooSimultaneous in order to calculate the joint likelihood of search and calibration data.
  RooSimultaneous *combinedModel = new RooSimultaneous("combinedModel","simultaneous pdf for the three sample types",*w->cat("sample")) ;
  combinedModel->addPdf(*w->pdf("calibNrModel"),"nr");
  combinedModel->addPdf(*w->pdf("calibErModel"),"er");
  combinedModel->addPdf(*w->pdf("modelWithConstraints"),"search");  
  w->import(*combinedModel);
  */
}

void ProfileLikelihood_ReducedParams(const char* fname=NULL, bool plcToyMC = kTRUE, bool showProjections = kFALSE, bool hypoTestInv = kFALSE)
{
    
  // to time the macro
  TStopwatch t;
  t.Start();
  
  // Create a local workspace to build the full model
  RooWorkspace *w = new RooWorkspace("w");

  RooDataSet* searchData;
  if(fname) {
    // Read in search, calib data
    TFile* f = TFile::Open(fname);
    TTree *tree = (TTree*)f->Get("t1");
    searchData = new RooDataSet("searchData","searchData",RooArgSet(S1,logS2S1,r,z),Import(*tree)) ;
  }

  // to randomise seed
  RooRandom::randomGenerator()->SetSeed(519);
  
  //////////////////////////////////////////////////////////////////////////////////
  // Build the model
  //////////////////////////////////////////////////////////////////////////////////
  
  buildModel(w);
  
  //////////////////////////////////////////////////////////////////////////////////
  // Configure model and save to workspace
  //////////////////////////////////////////////////////////////////////////////////
  
  ModelConfig* model = new ModelConfig("model", w);
  model->SetPdf(*w->pdf("modelWithConstraints"));
  //model->SetPdf(*w->pdf("combinedModel"));
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
  
  if(showProjections) {
    w->Print();

    // Generate data
    if(!fname)
      searchData = w->pdf("bgModel")->generate(*w->set("obs"), RooFit::Name("searchData"));

    //RooDataSet *sigData = w->pdf("sigModel")->generate(*w->set("obs"),RooFit::Name("sigData"));

    //RooDataSet *erData = w->pdf("calibErModel")->generate(*w->set("obs"), RooFit::Name("erDdata"));
    //RooDataSet *nrData = w->pdf("calibNrModel")->generate(*w->set("obs"), RooFit::Name("nrData"));

    // generate binned - not sure how to implement this with both search data (unbinned) and binned calib data
    //RooDataHist *erData = w->pdf("calibErModel")->generateBinned(*w->set("obs"),1000, RooFit::Name("erDdata"));
    //RooDataHist *nrData = w->pdf("calibNrModel")->generateBinned(*w->set("obs"),1000, RooFit::Name("nrData"));

//    RooDataSet *combinedData = new RooDataSet("combinedData","combined data",* w->set("obs"),Index(*w->cat("sample")),
//					      Import("er",*erData),Import("nr",*nrData),Import("search",*searchData)); 
    
    w->import(*searchData);
    //w->import(*combinedData);
    
    //w->pdf("combinedModel")->fitTo(*combinedData, ConditionalObservables(*w->var("S1")),RooFit::Constrain(*w->var("nCompton")),
//    w->pdf("combinedModel")->fitTo(*combinedData, RooFit::Constrain(*w->var("nCompton")),
//				   RooFit::Constrain(*w->var("nActivatedXe")), RooFit::Constrain(*w->var("nRn")) );
    
    TH1* hh_searchData = searchData->createHistogram("S1,logS2S1",100,100) ;
    //TH1* hh_sigData = sigData->createHistogram("S1",30) ;
    //TH1* hh_erData = erData->createHistogram("S1,logS2S1",100,100) ;
    //TH1* hh_nrData = nrData->createHistogram("S1,logS2S1",100,100) ;

    //hh_sigData->Draw();
    
    // Draw the data generated
    TCanvas* c0 = new TCanvas("c0","DataHistCanvas",1200,400) ;
    c0->Divide(3,1);
    c0->cd(1); gPad->SetLeftMargin(0.15) ; hh_searchData->GetZaxis()->SetTitleOffset(1.4) ; hh_searchData->Draw("col") ;
    //c0->cd(2); gPad->SetLeftMargin(0.15) ; hh_erData->GetZaxis()->SetTitleOffset(1.4) ; hh_erData->Draw("col") ;
    //c0->cd(3); gPad->SetLeftMargin(0.15) ; hh_nrData->GetZaxis()->SetTitleOffset(1.4) ; hh_nrData->Draw("col") ;
    
    // Plot data and PDF overlaid
    RooPlot* xframe = S1.frame(Title("S1 fit projection")) ;
    searchData->plotOn(xframe) ; 
    w->pdf("modelWithConstraints")->plotOn(xframe);
    w->pdf("sigModel")->plotOn(xframe,LineStyle(kDashed),LineColor(kRed));
    
    RooPlot* yframe = (*w->var("logS2S1")).frame(Title("log(S2/S1) fit projection")) ;
    searchData->plotOn(yframe) ;
    w->pdf("modelWithConstraints")->plotOn(yframe);
    w->pdf("sigModel")->plotOn(yframe,LineStyle(kDashed),LineColor(kRed));
    
    RooPlot* rframe = r.frame(Title("r fit projection")) ;
    searchData->plotOn(rframe) ;
    w->pdf("modelWithConstraints")->plotOn(rframe);
    w->pdf("sigModel")->plotOn(rframe,LineStyle(kDashed),LineColor(kRed));
    
    RooPlot* zframe = z.frame(Title("z fit projection")) ;
    searchData->plotOn(zframe) ;
    w->pdf("modelWithConstraints")->plotOn(zframe);
    w->pdf("sigModel")->plotOn(zframe,LineStyle(kDashed),LineColor(kRed));
    
    // Draw the frames on the canvas
    TCanvas* c1 = new TCanvas("c1","FitProjCanvas",800,600) ;
    c1->Divide(2,2);
    c1->cd(1);  gPad->SetLeftMargin(0.15) ; xframe->GetYaxis()->SetTitleOffset(1.4) ; xframe->Draw() ;
    c1->cd(2);  gPad->SetLeftMargin(0.15) ; yframe->GetYaxis()->SetTitleOffset(1.4) ; yframe->Draw() ;
    
    c1->cd(3);  gPad->SetLeftMargin(0.15) ; rframe->GetYaxis()->SetTitleOffset(1.4) ; rframe->Draw() ;
    c1->cd(4);  gPad->SetLeftMargin(0.15) ; zframe->GetYaxis()->SetTitleOffset(1.4) ; zframe->Draw() ;


//    ProfileLikelihoodCalculator plc(*searchData, *model);
//    //ProfileLikelihoodCalculator plc(*combinedData, *model);
//    
//    // ProfileLikelihoodCalculator makes central intervals
//    // to get one-sided upper-limit with desired CL, one needs to use twice larger alpha test size
//    if(setOneSided)
//      plc.SetTestSize(2.*(1.-limitCL));
//    else
//      plc.SetConfidenceLevel(limitCL);
//    
//    LikelihoodInterval* lrint = plc.GetInterval();  // that was easy.
//    //double ll=lrint->LowerLimit(*w->var("sig_xs"));
//    //double ul=lrint->UpperLimit(*w->var("sig_xs"));
//
//    double ul=lrint->UpperLimit(*w->var("nSig"));
//    //cout << "Profile likelihood limits: " << ll << "--" << ul <<endl;
//    cout << "Profile likelihood limit: " << ul <<endl;
//
//    //cout << "Signal xsec factor = " << w->var("nSigCubicCmPerKgPbDay")->getVal() * w->var("T")->getVal() * w->var("xeDensity")->getVal() << endl;
//
//    //erData->Delete();
//    //nrData->Delete();
//    searchData->Delete();
//    //combinedData->Delete();
  }

  // Write workspace to file
  w->writeToFile("ws_model.root");
  
  //////////////////////////////////////////////////////////////////////////////////
  // Try ProfileLikelihoodCalculator (Wilk's theorem)
  //////////////////////////////////////////////////////////////////////////////////
  
  if(plcToyMC) {
    w->Print();

    // Make ntuple to store values of UL
    TNtuple* mcTrials = new TNtuple("mcTrials","PLR limits from MC trials","upperLimit:searchEvents");

    // Generate pseudoexperiments for each trial and fit with ProfileLikelihoodCalculator
    for(int i = 0; i<ntrials;i++){
      w->loadSnapshot("initialConditions");
      //generate  one pseudoexperiment
      //RooDataSet *searchMcData = w->pdf("bgModel")->generate(*w->set("obs"), RooFit::Extended(),RooFit::Name("searchMcData"));
      RooDataSet *searchMcData = w->pdf("fullModel")->generate(*w->set("obs"), RooFit::Extended(),RooFit::Name("searchMcData"));
      //RooDataSet *erMcData = w->pdf("calibErModel")->generate(*w->set("obs"), RooFit::Extended(),RooFit::Name("erMcDdata"));
      //RooDataSet *nrMcData = w->pdf("calibNrModel")->generate(*w->set("obs"), RooFit::Extended(),RooFit::Name("nrMcData"));


      //RooDataSet *searchMcData = w->pdf("bgModelGen")->generate(*w->set("obs"), RooFit::Extended(),RooFit::Name("searchMcData"));
      //RooDataSet *erMcData = w->pdf("calibErModelGen")->generate(*w->set("obs"), RooFit::Extended(),RooFit::Name("erMcDdata"));
      //RooDataSet *nrMcData = w->pdf("calibNrModelGen")->generate(*w->set("obs"), RooFit::Extended(),RooFit::Name("nrMcData"));

//      RooDataSet *combinedMcData = new RooDataSet("combinedData","combined data",* w->set("obs"),Index(*w->cat("sample")),
//      					  Import("er",*erMcData),Import("nr",*nrMcData),Import("search",*searchMcData)); 
      
      cout << "*** generated search data, numEntries = " << searchMcData->numEntries() << endl;
      //cout << "*** generated ER calib data, numEntries = " << erMcData->numEntries() << endl;
      //cout << "*** generated NR calib data, numEntries = " << nrMcData->numEntries() << endl;


      ProfileLikelihoodCalculator plc(*searchMcData, *model);
      //ProfileLikelihoodCalculator plc(*combinedMcData, *model);
      
      // ProfileLikelihoodCalculator makes central intervals
      // to get one-sided upper-limit with desired CL, one needs to use twice larger alpha test size
      if(setOneSided)
	plc.SetTestSize(2.*(1.-limitCL));
      else
	plc.SetConfidenceLevel(limitCL);
      
      LikelihoodInterval* lrint = plc.GetInterval();  // that was easy.
      //double ll=lrint->LowerLimit(*w->var("sig_xs"));
      //double ul=lrint->UpperLimit(*w->var("sig_xs"));

      //double ul=lrint->UpperLimit(*w->var("sig_xs"))*(w->var("nSigCubicCmPerKgPbDay")->getVal() * w->var("T")->getVal() * w->var("xeDensity")->getVal());
      
      double ul=lrint->UpperLimit(*w->var("nSig"));
      //cout << "Profile likelihood limits: " << ll << "--" << ul <<endl;
      cout << "Profile likelihood limit: " << ul <<endl;

      //cout << "Signal xsec factor = " << w->var("nSigCubicCmPerKgPbDay")->getVal() * w->var("T")->getVal() * w->var("xeDensity")->getVal() << endl;
 
      //Fill the tree, tidy
      //later could save three TH2Ds for each trial?
      //mcTrials->Fill(ll,ul,searchMcData->numEntries());
      mcTrials->Fill(ul,searchMcData->numEntries());
      //erMcData->Delete();
      //nrMcData->Delete();
      searchMcData->Delete();
      //combinedMcData->Delete();
      
    }
    
    //plot and save results
    //mcTrials->Draw("upperLimit:lowerLimit>>hist(10,-0.5,9.5,20,-0.5,19.5)","","colz");
    //mcTrials->Draw("upperLimit>>hist(100,1e-10,1e-08)","","");
    mcTrials->Draw("upperLimit>>histUL(100,0,25)","","");
    mcTrials->Draw("searchEvents>>hist(100,0,100)","","");
    mcTrials->SaveAs( TString::Format("plr_reducedParams_%.0fgev_%0.fd_r%0.fcm.root",
				      w->var("mWimp")->getVal(),w->var("T")->getVal(),w->var("r")->getMax()) );
				      //w->var("nErCalib")->getVal(),w->var("nNrCalib")->getVal()) );
  }

  //////////////////////////////////////////////////////////////////////////////////
  // Hypothesis test inversion
  //////////////////////////////////////////////////////////////////////////////////

  // Note - doesn't work with simultaneous pdf right now! 
  if(hypoTestInv) {

    // Generate data
    searchData = w->pdf("bgModel")->generate(*w->set("obs"), RooFit::Name("searchData"));
    //RooDataSet *erData = w->pdf("calibErModel")->generate(*w->set("obs"), RooFit::Name("erDdata"));
    //RooDataSet *nrData = w->pdf("calibNrModel")->generate(*w->set("obs"), RooFit::Name("nrData"));
//    RooDataSet *combinedData = new RooDataSet("combinedData","combined data",* w->set("obs"),Index(*w->cat("sample")),
//					      Import("er",*erData),Import("nr",*nrData),Import("search",*searchData)); 
    
    //FrequentistCalculator *  fc  = new FrequentistCalculator(*combinedData, *bModel, *model);
    FrequentistCalculator *  fc  = new FrequentistCalculator(*searchData, *bModel, *model);
    fc->SetToys(ntrials,ntrials/2.);    // 1000 for null (S+B) , 50 for alt (B)
    
    // create hypotest inverter 
    HypoTestInverter calc(*fc);
    
    // set confidence level (e.g. 95% upper limits)
    calc.SetConfidenceLevel(limitCL);
   
    // for CLS
    bool useCLs = false;
    calc.UseCLs(useCLs);
    calc.SetVerbose(false);
    
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
    
    //cout << "Doing an auto scan." << endl;
    //calc.SetAutoScan();
  
    HypoTestInverterResult * result = calc.GetInterval();
    
    double upperLimit = result->UpperLimit();
    double ulError = result->UpperLimitEstimatedError();
    cout << "The computed upper limit is: " << upperLimit << " +/- " << ulError << endl;
  
    // compute expected limit
    cout << "Expected upper limits, using the B (alternate) model : " << endl;
    cout << " expected limit (median) " << result->GetExpectedUpperLimit(0) << endl;
    cout << " expected limit (-1 sig) " << result->GetExpectedUpperLimit(-1) << endl;
    cout << " expected limit (+1 sig) " << result->GetExpectedUpperLimit(1) << endl;
    
    // plot now the result of the scan 
    HypoTestInverterPlot *plot = new HypoTestInverterPlot("HTI_Result_Plot","Feldman-Cousins Interval",result);
    
    // plot in a new canvas with style
    TCanvas * c1 = new TCanvas("HypoTestInverter Scan"); 
    c1->SetLogy(false);
    
    plot->Draw("CLb 2CL");  // plot also CLb and CLs+b 
    //plot->Draw("OBS");  // plot only observed p-value
    

    // plot also in a new canvas the test statistics distributions 
    
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
      }
   }
   
   //erData->Delete();
   //nrData->Delete();
   searchData->Delete();
   //combinedData->Delete();
  }
  
  // print timing info
  t.Stop();
  t.Print();
}
