//////////////////////////////////////////////////////////////////////////
//Based on KRG's ProfileLikelihood.C
//AC 2013-07-16 -- Put in simple linear yields of initial quanta and Gaussian resolutions on S1 and S2. Reverted to a beta=0 version of WIMP spectrum until normalisation of full version is fixed.
//////////////////////////////////////////////////////////////////////////
//
// Setting up an extended maximum likelihood fit
//
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooProjectedPdf.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TFile.h"
#include "TStopwatch.h"
#include "RooPlot.h"
#include "RooWorkspace.h"
#include "RooWimpSpectrum.h"
#include "RooS1Response.h"
#include "ProfileLikelihood.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodIntervalPlot.h"
#include "TROOT.h"

#include "RooWimpSpec_ctsPerkeVkgdaypb.h"

#include <string>

using namespace std;

using namespace RooFit ;
using namespace RooStats;

void buildModel(RooWorkspace* w)
{

   /////////////////////////////////////////
   // The Model building stage
   /////////////////////////////////////////

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
   w->import(Enr);
   w->import(Eee);
   w->import(nPhot);
   w->import(nE);
   w->import(S1);
   w->import(S2);

   w->import(rateCompton);
   w->import(ra);
   w->import(za);
   w->import(zb);

   w->import(rateActivatedXe);
   w->import(rateRn);
   w->import(xeDensity);

   //w->import(nSig);
   w->import(nCompton);
   w->import(nActivatedXe);
   w->import(nRn);

   w->import(relUncertaintyBg);

   w->defineSet("obs","S1,S2,r,z");
   //w->defineSet("obs","Enr,nPhot,nE,r,z");
   w->defineSet("projNR","Enr,nPhot,nE");
   w->defineSet("projER","Eee,nPhot,nE");
   w->defineSet("s1Cond","nPhot,r,z");
   w->defineSet("s2Cond","nE,r,z");
   w->defineSet("nuis","nCompton,nActivatedXe,nRn");

   ////////////////////////////////////////////////////////////////////////////////////////
   // Build up model with terms
   // P(E,r,z)P(nPhot,nE|E)P(S1|nPhot,r,z)P(S2|nE,r,z,t)
   // for both signal & bg
   ////////////////////////////////////////////////////////////////////////////////////////

   /////////////////////////////////////
   // Signal model
   /////////////////////////////////////

   // Declare pdf of WIMP energy spectrum P(Enr)
   //w->factory("WimpSpectrum::dRdE(Enr, mWimp, v0, vE, vEsc, rho, c, r, a, s, delta, beta)");
   //gROOT->ProcessLine(".L RooWimpSpec_ctsPerkeVkgdaypb.cxx+");//load a custom library
   //w->importClassCode("RooWimpSpec_ctsPerkeVkgdaypb",kTRUE);//import custom class into workspace
   //w->factory("WimpSpec_ctsPerkeVkgdaypb::dRdE(Enr,mWimp,targetA,targetZ,rho,dayOfYear,vEsc,v0,ff_c,ff_a,ff_s,idm_delta)");

   gROOT->ProcessLine(".L RooWimpSpectrum.cxx+");//load a custom library
   w->importClassCode("RooWimpSpectrum", kTRUE);//import custom class into workspace
   w->factory("WimpSpectrum::dRdE(Enr,mWimp, v0, vE, vEsc, rho, ff_c, ff_r0, ff_a, ff_s, idm_delta, halo_beta)");

   // Create pdf P(nPhot,nE|Enr)
   w->factory("cexpr::toyMeanPhotNr('toyPhotNrYield*Enr',toyPhotNrYield[8.0],Enr)");
   w->factory("cexpr::toyMeanElecNr('toyElecNrYield*Enr',toyElecNrYield[4.0],Enr)");
   w->factory("cexpr::sigmaPhotNr('sqrt(nPhot)',nPhot)");
   w->factory("cexpr::sigmaElecNr('sqrt(nE)',nE)");
   w->factory("PROD::NrQuantaModel(Poisson::photNr(nPhot,toyMeanPhotNr),Poisson::elecNr(nE,toyMeanElecNr))");
   //w->factory("PROD::NrQuantaModel(Gaussian::photNr(nPhot,toyMeanPhotNr,sigmaPhotNr),Gaussian::elecNr(nE,toyMeanElecNr,sigmaElecNr))");
   //To do: replace toy above with a real model below...
   //w->factory("QuantaNR::NrQuantaModel(Enr, nPhot, nE)");

   // Declare P(r) and P(z) for WIMPs  - only need to do this once for all other uniform pdfs
   w->factory("CEXPR::UniformR('TMath::TwoPi()*r',r)");
   w->factory("Uniform::UniformZ(z)");

   // Multiply P(nPhot,nE|Enr)P(Enr)P(r)P(z)
   w->factory("PROD::SigQuantaModel(NrQuantaModel|Enr, dRdE, UniformR, UniformZ)");

   // Create pdf P(S1|nPhot,r,z) - only need to do this once for all future pdfs
   w->factory("cexpr::toyS1res('0.1+sqrt(nPhot/toyPheYield)',toyPheYield[0.14],nPhot)");
   w->factory("Gaussian::S1Res(S1,nPhot,toyS1res)");
   //To do: replace toy above with a real model below...
   // w->factory("S1Response::S1Res(S1,nPhot,r,z)");

   // Create pdf P(S2|nE,r,z,t) - only need to do this once for all future pdfs
   w->factory("cexpr::toyS2res('0.1+sqrt(nE/toySEYield)',toySEYield[0.60],nE)");
   w->factory("Gaussian::S2Res(S2,nE,toyS2res)"); //pure gaussian S2 resolution
   //To do: replace toy above with a real model below...
   //w->factory("S2Response::S2Res(S2,nPhot,r,z)");

   // Multiply to get full signal model
   w->factory("PROD::SigModel(S1Res|nPhot, S2Res|nE, SigQuantaModel)");
   //To do: replace toy above once S1Res and S2Res have position-dependence...
   //w->factory("PROD::SigModel(S1Res|s1Cond, S2Res|s2Cond, SigQuantaModel)");

   // Marginalize over Enr, nPhot, nE
   RooProjectedPdf* SigModelProj = (RooProjectedPdf*) w->pdf("SigModel")->createProjection(*w->set("projNR"));
   SigModelProj->SetName("SigModelProj");
   w->import(*SigModelProj);

   //Make a variable corresponding to the number of signal events
   RooAbsReal* nSigPerPbDay =SigModelProj->createIntegral( *w->set("obs") ); // [day-1 pb-1 cm3]
   nSigPerPbDay->SetName("nSigPerPbDay");
   w->import(*nSigPerPbDay);
   w->factory("prod::nSig(sig_xs,T, xeDensity, nSigPerPbDay)"); // number of signal events 


   /////////////////////////////////////
   // Background model
   /////////////////////////////////////

   // ER background from detector materials
   w->factory("Uniform::ComptonSpectrum(Eee)");
   w->factory("cexpr::toyMeanPhotEr('toyPhotErYield*Eee',toyPhotErYield[30.0],Eee)");
   w->factory("cexpr::toyMeanElecEr('toyElecErYield*Eee',toyElecErYield[30.0],Eee)");
   w->factory("cexpr::sigmaPhotEr('sqrt(nPhot)',nPhot)");
   w->factory("cexpr::sigmaElecEr('sqrt(nE)',nE)");

   w->factory("PROD::ErQuantaModel(Poisson::photEr(nPhot,toyMeanPhotEr),Poisson::elecEr(nE,toyMeanElecEr))");
   //w->factory("PROD::ErQuantaModel(Gaussian::photEr(nPhot,toyMeanPhotNr,sigmaPhotEr),Gaussian::elecEr(nE,toyMeanElecEr,sigmaElecEr))");
   //w->factory("QuantaER::ErQuantaModel(Eee, nPhot, nE)");

   // Dave's model of spatial dependence of Compton scatters from materials
   w->factory("CEXPR::ComptonSpatial('TMath::TwoPi()*r*rateCompton*exp(pow(r/ra, 2.)+pow((z-za)/zb, 2.))', r, rateCompton, ra, z , za, zb)");

   w->factory("PROD::ComptonQuantaModel(ErQuantaModel|Eee,ComptonSpectrum,ComptonSpatial)");
   //w->factory("PROD::ComptonQuantaModel(ErQuantaModel|Eee,ComptonSpectrum)");

//   // Multiply to get full Compton model
//   w->factory("PROD::ComptonModel(S1Res|s1Cond,S2Res|s2Cond,ComptonQuantaModel)");
//
//   // Marginalize over Eee, nPhot, nE
//   RooProjectedPdf* ComptonModelProj = (RooProjectedPdf*) w->pdf("ComptonModel")->createProjection(*w->set("projER"));
//   ComptonModelProj->SetName("ComptonModelProj");
//   w->import(*ComptonModelProj);
//
//   // Calculate number of expected Compton events due to materials
//   RooRealVar nComptonPerDay = ComptonModelProj->createIntegral(*w->set("obs")); // [day-1]
//
//   nComptonPerDay.SetName("nComptonPerDay");
//   w->import(nComptonPerDay);
//   w->factory("prod::nComptonPred(T, xeDensity, nComptonPerDay)"); // number of Compton events predicted by bg model

   // Calculate number of expected Compton events due to materials
   RooAbsReal* nComptonPerDay = w->pdf("ComptonQuantaModel")->createIntegral(*w->set("obs")); // [day-1]
   nComptonPerDay->SetName("nComptonPerDay");
   w->import(*nComptonPerDay);
   w->factory("prod::nComptonPred(T, xeDensity, nComptonPerDay)"); // number of Compton events predicted by bg model

   w->factory("prod::sigmaCompton(nComptonPred,relUncertaintyBg)"); // error on number of Compton events predicted by bg model

   // ER background from Xe-127
   w->factory("Uniform::ActivatedXeSpectrum(Eee)");
   w->factory("CEXPR::ActivatedXeSpatial('TMath::TwoPi()*r*rateActivatedXe',r,rateActivatedXe,z)");
   w->factory("PROD::ActivatedXeQuantaModel(ErQuantaModel|Eee,ActivatedXeSpectrum,ActivatedXeSpatial)");

//   // Multiply to get full Xe-127 model
//   w->factory("PROD::ActivatedXeModel(S1Res|s1Cond,S2Res|s2Cond,ActivatedXeQuantaModel)");
//
//   // Marginalize over Eee, nPhot, nE
//   RooProjectedPdf* ActivatedXeModelProj = (RooProjectedPdf*) w->pdf("ActivatedXeModel")->createProjection(*w->set("projER"));
//   ActivatedXeModelProj->SetName("ActivatedXeModelProj");
//   w->import(*ActivatedXeModelProj);
//
//   // Calculate number of expected Xe-127 events
//   RooRealVar nActivatedXePerDay = ActivatedXeModelProj->createIntegral(*w->set("obs")); // [day-1]
//
//   nActivatedXePerDay.SetName("nActivatedXePerDay");
//   w->import(nActivatedXePerDay);
//   w->factory("prod::nActivatedXePred(T, xeDensity, nActivatedXePerDay)"); // number of Xe-127 events predicted by bg model

   // Calculate number of expected Xe-127 events
   RooAbsReal* nActivatedXePerDay = w->pdf("ActivatedXeQuantaModel")->createIntegral(*w->set("obs")); // [day-1]
   nActivatedXePerDay->SetName("nActivatedXePerDay");
   w->import(*nActivatedXePerDay);
   w->factory("prod::nActivatedXePred(T, xeDensity, nActivatedXePerDay)"); // number of Xe-127 events predicted by bg model

   w->factory("prod::sigmaActivatedXe(nActivatedXePred,relUncertaintyBg)"); // error on number of Xe-127 events predicted by bg model

   // ER background from Rn-222 and Kr-85
   w->factory("Uniform::RnSpectrum(Eee)");
   w->factory("CEXPR::RnSpatial('TMath::TwoPi()*r*rateRn',r,rateRn,z)");
   w->factory("PROD::RnQuantaModel(ErQuantaModel|Eee,RnSpectrum,RnSpatial)");

//   // Multiply to get full Rn model
//   w->factory("PROD::RnModel(S1Res|s1Cond,S2Res|s2Cond,RnQuantaModel)");
//
//   // Marginalize over Eee, nPhot, nE
//   RooProjectedPdf* RnModelProj = (RooProjectedPdf*) w->pdf("RnModel")->createProjection(*w->set("projER"));
//   RnModelProj->SetName("RnModelProj");
//   w->import(*RnModelProj);
//
//   // Calculate number of expected Xe-127 events
//   RooRealVar nRnPerDay = RnModelProj->createIntegral(*w->set("obs")); // [day-1]
//
//   nRnPerDay.SetName("nRnPerDay");
//   w->import(nRnPerDay);
//   w->factory("prod::nRnPred(T, xeDensity, nRnPerDay)"); // number of Xe-127 events predicted by bg model

   // Calculate number of expected Xe-127 events
   RooAbsReal* nRnPerDay = w->pdf("RnQuantaModel")->createIntegral(*w->set("obs")); // [day-1]
   nRnPerDay->SetName("nRnPerDay");
   w->import(*nRnPerDay);
   w->factory("prod::nRnPred(T, xeDensity, nRnPerDay)"); // number of Xe-127 events predicted by bg model

   w->factory("prod::sigmaRn(nRnPred,relUncertaintyBg)"); // error on number of Xe-127 events predicted by bg model

   // Put together complete likelihood
   //w->factory("SUM::fullModel(nSig*SigModelProj,nCompton*ComptonModelProj,nActivatedXe*ActivatedXeModelProj,nRn*RnModelProj)");
   w->factory("SUM::fullModel(nSig*SigQuantaModel,nCompton*ComptonQuantaModel,nActivatedXe*ActivatedXeQuantaModel,nRn*RnQuantaModel)");
   w->factory("SUM::bgModel(nCompton*ComptonQuantaModel,nActivatedXe*ActivatedXeQuantaModel,nRn*RnQuantaModel)");

   // Add constraints on nuisance parameters
   w->factory("Gaussian::comptonConstraint(nComptonPred, nCompton, sigmaCompton)");
   w->factory("Gaussian::activatedXeConstraint(nActivatedXePred, nActivatedXe, sigmaActivatedXe)");
   w->factory("Gaussian::rnConstraint(nRnPred, nRn, sigmaRn)");

   // Full model with constraints
   w->factory("PROD::modelWithConstraints(fullModel,comptonConstraint,activatedXeConstraint,rnConstraint)"); // product of terms

   w->importClassCode();
   
}

void ProfileLikelihood(const char* fileName, bool showProjections = kFALSE)
{
    
    // to time the macro
    TStopwatch t;
    t.Start();
    
    // Create a local workspace to build the full model
    RooWorkspace *w = new RooWorkspace("w");
    
    // Build the model
    buildModel(w);
    
    // Read models back to local space
    RooAbsPdf* fullModel   = w->pdf("fullModel");
    RooAbsPdf* modelWithConstraints   = w->pdf("modelWithConstraints");
    //RooAbsPdf* sigModel = w->pdf("SigModelProj");
    RooAbsPdf* sigModel = w->pdf("SigQuantaModel");
    RooAbsPdf* bgModel   = w->pdf("bgModel");
    
    
    // Generate data from full model
    RooDataSet* data = bgModel->generate(*w->set("obs"),Name("data"));
    w->import(*data);
    
    // Print structure of composite p.d.f.
    //modelWithConstraints->Print("t");
    
    if(showProjections) {
        //fullModel->fitTo(*data);
      modelWithConstraints->fitTo(*data, RooFit::Constrain(*w->var("nCompton")),RooFit::Constrain(*w->var("nActivatedXe")),RooFit::Constrain(*w->var("nRn")) );
        
      // Plot data and PDF overlaid
      //RooPlot* xframe = S1.frame(Title("S1 fit projection")) ;
      RooPlot* xframe = nPhot.frame(Title("S1 fit projection")) ;
      data->plotOn(xframe) ;
      sigModel->plotOn(xframe,LineStyle(kDashed),LineColor(kRed));
      //fullModel->plotOn(xframe);
      modelWithConstraints->plotOn(xframe);
        
      //RooPlot* yframe = S2.frame(Title("S2 fit projection")) ;
      //RooPlot* yframe = logS2S1.frame(Title("S2 fit projection")) ;
      RooPlot* yframe = nE.frame(Title("S2 fit projection")) ;
      data->plotOn(yframe) ;
      sigModel->plotOn(yframe,LineStyle(kDashed),LineColor(kRed));
      //fullModel->plotOn(yframe);
      modelWithConstraints->plotOn(yframe);
      
      RooPlot* rframe = r.frame(Title("r fit projection")) ;
      data->plotOn(rframe) ;
      sigModel->plotOn(rframe,LineStyle(kDashed),LineColor(kRed));
      //fullModel->plotOn(rframe);
      modelWithConstraints->plotOn(rframe);
      
      RooPlot* zframe = z.frame(Title("z fit projection")) ;
      data->plotOn(zframe) ;
      sigModel->plotOn(zframe,LineStyle(kDashed),LineColor(kRed));
      //fullModel->plotOn(zframe);
      modelWithConstraints->plotOn(zframe);
      
//      // Draw the data generated
//      TCanvas* c0 = new TCanvas("c0","DataHistCanvas",800,600) ;
//      gPad->SetLeftMargin(0.15) ; hh_data->GetZaxis()->SetTitleOffset(1.4) ; hh_data->Draw("col") ;
//      //gPad->SetLeftMargin(0.15) ; hh_Edata->GetZaxis()->SetTitleOffset(1.4) ; hh_Edata->Draw("col") ;
      
      // Draw the frames on the canvas
      TCanvas* c1 = new TCanvas("c1","FitProjCanvas",800,600) ;
      c1->Divide(2,2);
      c1->cd(1);  gPad->SetLeftMargin(0.15) ; xframe->GetYaxis()->SetTitleOffset(1.4) ; xframe->Draw() ;
      c1->cd(2);  gPad->SetLeftMargin(0.15) ; yframe->GetYaxis()->SetTitleOffset(1.4) ; yframe->Draw() ;
      
      c1->cd(3);  gPad->SetLeftMargin(0.15) ; rframe->GetYaxis()->SetTitleOffset(1.4) ; rframe->Draw() ;
      c1->cd(4);  gPad->SetLeftMargin(0.15) ; zframe->GetYaxis()->SetTitleOffset(1.4) ; zframe->Draw() ;
    }
    
    
    //////////////////////////////////////////////////////////////////////////////////
    // Configure model and save to workspace
    //////////////////////////////////////////////////////////////////////////////////
    
    ModelConfig* myModel = new ModelConfig("myModel", w);
    //myModel.SetPdf(*fullModel);
    myModel->SetPdf(*modelWithConstraints);
    myModel->SetObservables(*w->set("obs"));
    myModel->SetParametersOfInterest(*w->var("sig_xs"));
    myModel->SetNuisanceParameters(*w->set("nuis"));
    
    w->import(*myModel);
    
    RooRealVar* poi = (RooRealVar*) myModel->GetParametersOfInterest()->first();
    
    ModelConfig*  bModel = (ModelConfig*) myModel->Clone();
    bModel->SetName("myModel_with_poi_0");
    double oldval = poi->getVal();
    poi->setVal(0);
    bModel->SetSnapshot( *poi  );
    poi->setVal(oldval);
    
    w->import(*bModel);
    
    //w->Print();
    w->writeToFile("ws_test.root");
    
    // First, let's use a Calculator based on the Profile Likelihood Ratio
    ProfileLikelihoodCalculator plc(*data, *myModel);
    
    // ProfileLikelihoodCalculator makes central intervals
    // to get one-sided upper-limit with desired CL, one needs to use twice larger alpha test size
    //plc.SetTestSize(2.*(1.-limitCL));
    plc.SetConfidenceLevel(limitCL);
    
    RooArgSet* nullParams = (RooArgSet*) myModel->GetParametersOfInterest();
    nullParams->setRealValue("sig_xs",0); 
    plc.SetNullParameters( *nullParams);
    
    ConfInterval* lrint = plc.GetInterval();  // that was easy.
    
    // Let's make a plot
    TCanvas* dataCanvas = new TCanvas("dataCanvas");
    LikelihoodIntervalPlot plotInt((LikelihoodInterval*)lrint);
    plotInt.SetTitle("Profile Likelihood Ratio and Posterior for sig_xs");
    plotInt.Draw();
    
    // Get upper limit from Profile Calculator
    cout << endl;
    cout << "Profile upper limit on s = " << ((LikelihoodInterval*) lrint)->UpperLimit(sig_xs) << endl;
    
    /// print timing info
    t.Stop();
    t.Print();
}
