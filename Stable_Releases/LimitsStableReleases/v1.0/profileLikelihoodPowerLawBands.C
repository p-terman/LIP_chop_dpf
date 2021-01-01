//A program to investigate the effect of MC sample size on PLR limits
//Basic idea is to work in log10(S2/S1)-vs-S1 and have 3 parameters each for the S1-dependence of mu and sigma for ER and NE (A*Exp(B*S1)+C*S1 . 

#include "TStopwatch.h"
#include "TCanvas.h"
#include "TNtuple.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"

#include <iostream>

#include "RooRandom.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "RooStats/ProfileLikelihoodCalculator.h"

void profileLikelihoodPowerLawBands(const double muCalib=200., const double muSearch=60. ,const int ntrials=100)
{
  using namespace RooFit;
  using namespace RooStats;
  // to randomise seed
  RooRandom::randomGenerator()->SetSeed(0);
  RooWorkspace *w = new RooWorkspace();
   
//   TFile *_file0 = TFile::Open("pp60keVnr_LCp010-Cath10.root");
//   TFile *_file1 = TFile::Open("pp15keVgr_LCp010-Cath10.root");
//   
//   TTree*nr = (TTree*)_file0->Get("observedTree");
//   TTree*er = (TTree*)_file1->Get("observedTree");
//   
//   RooRealVar fS1_phe("fS1_phe","fS1_phe",0.1,50);
//   RooRealVar fS2_e("fS2_e","S2",0.1,1000);
//   
//   RooDataSet nestNr("nestNr","nestNr",nr,RooArgSet(fS1_phe,fS2_e));
//   RooDataSet nestEr("nestEr","nestEr",er,RooArgSet(fS1_phe,fS2_e));
// 
//   w->import(nestEr);
//   w->import(nestNr);
  
  //w->factory("cexpr::scattY('log10(fS2_e/fS1_phe)',fS1_phe,fS2_e)"); uncomment to use S2 rather than log10(S2/S1) as the RooRealVar
  w->factory("scattY[-2,2]");
  w->factory("fS1_phe[0.1,50]");
  
  //build functions for median log10(S2/S1) (aka scattY) as function of S1. The numerical values are from fits to NEST sims after applying 14%LC, 30%SPEres, 60%CC, 20%SEres.
  w->factory("cexpr::meanEr('ermA*pow(fS1_phe,ermB)+ermC*fS1_phe',ermA[1.71,1.0,2.0],fS1_phe,ermB[-0.161,-0.2,0],ermC[-0.0056,-0.01,0])");
  //w->factory("cexpr::meanEr('ermA*pow(fS1_phe,ermB)+ermC*fS1_phe',ermA[1.71],fS1_phe,ermB[-0.161],ermC[-0.0056])"); //try a fixed band
  w->factory("cexpr::meanNr('nrmA*pow(fS1_phe,nrmB)+nrmC*fS1_phe',nrmA[1.12,1.0,2.0],fS1_phe,nrmB[-0.124,-0.2,0],nrmC[-0.005,-0.01,0])");
  //w->factory("cexpr::meanNr('nrmA*pow(fS1_phe,nrmB)+nrmC*fS1_phe',nrmA[1.12],fS1_phe,nrmB[-0.124],nrmC[-0.005])"); //try a fixed band
  
  //build functions for gaussian sigma in log10(S2/S1) as function of S1 (via parametrisation of mu+1sigma)
  w->factory("cexpr::sigmaEr('ersA*pow(fS1_phe,ersB)-ermA*pow(fS1_phe,ermB)+(ersC-ermC)*fS1_phe',ermA, ermB, ermC,ersA[1.95,1.0,3.0],fS1_phe,ersB[-0.175,-0.2,0],ersC[-0.0032,-0.01,0])");
  //w->factory("cexpr::sigmaEr('ersA*pow(fS1_phe,ersB)-ermA*pow(fS1_phe,ermB)+(ersC-ermC)*fS1_phe',ermA, ermB, ermC,ersA[1.95],fS1_phe,ersB[-0.175],ersC[-0.0032])"); // try with fixed band
  w->factory("cexpr::sigmaNr('nrsA*pow(fS1_phe,nrsB)-nrmA*pow(fS1_phe,nrmB)+(nrsC-nrmC)*fS1_phe',nrmA, nrmB, nrmC,nrsA[1.47,1.0,3.0],fS1_phe,nrsB[-0.195,-0.2,0],nrsC[-0.0034,-0.01,0])");
  //w->factory("cexpr::sigmaNr('nrsA*pow(fS1_phe,nrsB)-nrmA*pow(fS1_phe,nrmB)+(nrsC-nrmC)*fS1_phe',nrmA, nrmB, nrmC,nrsA[1.47],fS1_phe,nrsB[-0.195],nrsC[-0.0034])"); //try with fixed band
  
  w->factory("Gaussian::scattYEr(scattY,meanEr,sigmaEr)");
  w->factory("Gaussian::scattYNr(scattY,meanNr,sigmaNr)");
  
  
  //build distributions in log10(S2/S1) -v- S1 for NR and ER assuming flat in S1
  w->factory("PROD::pdfNr(Exponential::almostLikeAWIMP(fS1_phe,decayConst[-0.1]),scattYNr|fS1_phe)");
  w->factory("PROD::pdfEr(Uniform::s1specEr(fS1_phe),scattYEr|fS1_phe)");
  
  //extend the distributions
  w->factory("ExtendPdf::densityNr(pdfNr,calibNr[10000])");
  w->var("calibNr")->setVal(muCalib);
  w->factory("ExtendPdf::densityEr(pdfEr,calibEr[10000])");
  w->var("calibEr")->setVal(muCalib);
  
  w->factory("SUM::densitySearch(nNrInSearch[0,0,1000]*pdfNr,nErInSearch[100,0,1000]*pdfEr)");
  w->var("nErInSearch")->setVal(muSearch);
  w->factory("sample[nr,er,search]");
  
  //Save the MC truth model, to be loaded for each pseudoexperiment
  w->saveSnapshot("initialConditions",w->allVars());

  //make a RooSimultaneous in order to calculate the joint likelihood of search and calibration data.
  RooSimultaneous *combPdf = new RooSimultaneous("combPdf","simultaneous pdf for the three sample types",*w->cat("sample")) ;
  combPdf->addPdf(* w->pdf("densityEr"),"er");
  combPdf->addPdf(* w->pdf("densityNr"),"nr");
  combPdf->addPdf(* w->pdf("densitySearch"),"search");  
  w->import(*combPdf);
  
  //configure the model's observable, interesting and nuisance parameters
  w->defineSet("observables","fS1_phe,scattY");
  w->defineSet("poi","nNrInSearch");
  w->defineSet("np","nErInSearch,ermA,ermB,ermC,ersA,ersB,ersC,nrmA,nrmB,nrmC,nrsA,nrsB,nrsC");
  //w->defineSet("np","nErInSearch"); //try a fixed band
  ModelConfig *model=new ModelConfig("model",w);
  model->SetParametersOfInterest(* w->set("poi"));
  model->SetObservables(* w->set("observables"));
  model->SetPdf(* w->pdf("combPdf"));
  w->import(*model);
  
  TNtuple* mcTrials = new TNtuple("mcTrials","PLR limits from MC trials","lowerLimit:upperLimit:searchEvents");

  for(int i = 0; i<ntrials;i++){
    w->loadSnapshot("initialConditions");
    //generate  one pseudoexperiment
    RooDataSet *er_mcdata = w->pdf("densityEr")->generate(*w->set("observables"),RooFit::Extended(), RooFit::Name("er_mcdata"));
    RooDataSet *nr_mcdata = w->pdf("densityNr")->generate(*w->set("observables"),RooFit::Extended(), RooFit::Name("nr_mcdata"));
    RooDataSet *search_mcdata = w->pdf("densitySearch")->generate(*w->set("observables"),RooFit::Extended(), RooFit::Name("search_mcdata"));
    cout<<"*** generated search data, numEntries = "<<search_mcdata->numEntries()<<endl;
    RooDataSet *combData = new RooDataSet("combData","combined data",* w->set("observables"),Index(*w->cat("sample")),Import("er",*er_mcdata),Import("nr",*nr_mcdata),Import("search",*search_mcdata)); 

    //calculate a Wilks' theorem confidence interval for the number of signal events
    ProfileLikelihoodCalculator plc(*combData,*model);
    plc.SetConfidenceLevel(0.9);
    LikelihoodInterval* lrint = plc.GetInterval();
    double ll=lrint->LowerLimit(* w->var("nNrInSearch"));
    double ul=lrint->UpperLimit(* w->var("nNrInSearch"));
    cout<<ll<<"--"<<ul<<endl;
    
    //Fill the tree, tidy
    //later could save three TH2Ds for each trial?
    mcTrials->Fill(ll,ul,search_mcdata->numEntries());
    er_mcdata->Delete();
    nr_mcdata->Delete();
    search_mcdata->Delete();
    combData->Delete();

  }
  //plot and save results
  mcTrials->Draw("upperLimit:lowerLimit>>hist(10,-0.5,9.5,20,-0.5,19.5)","","colz");
  mcTrials->SaveAs(TString::Format("powerlaw%.0fcal%.0fobs.root",muCalib,muSearch));
}
