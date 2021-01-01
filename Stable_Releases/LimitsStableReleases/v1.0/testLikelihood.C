#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "TH2F.h"
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
#include "RooLightCollection.h"

#include <string>

using namespace std;

using namespace RooFit ;
using namespace RooStats;

void testLikelihood(){
   RooWorkspace *w = new RooWorkspace();
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
   w->import(ff_c);
   //w->import(ff_r0);
   w->import(ff_a);
   w->import(ff_s);
   w->import(idm_delta);
   //w->import(halo_beta);
   w->import(Enr);
   w->import(Eee);
   w->import(nPhot);
   w->import(nE);
   w->import(S1);
   w->import(S2);

   w->import(rCompton);
   w->import(ra);
   w->import(za);
   w->import(zb);

   w->import(nSig);
   w->import(nCompton);
   w->import(nActivatedXe);
   w->import(nRn);

   w->defineSet("obs","S1,S2,r,z");
   w->defineSet("projNR","Enr,nPhot,nE");
   w->defineSet("projER","Eee,nPhot,nE");
   w->defineSet("s1Cond","nPhot,r,z");
   w->defineSet("s2Cond","nE,r,z");

   gROOT->ProcessLine(".L RooWimpSpec_ctsPerkeVkgdaypb.cxx+");//load a custom library
   w->importClassCode("RooWimpSpec_ctsPerkeVkgdaypb",kTRUE);//import custom class into workspace
   w->factory("WimpSpec_ctsPerkeVkgdaypb::dRdE(Enr,mWimp,targetA,targetZ,rho,dayOfYear,vEsc,v0,ff_c,ff_a,ff_s,idm_delta)");
   w->factory("cexpr::toyMeanPhotNr('toyPhotNrYield*Enr',toyPhotNrYield[8.0],Enr)");
   w->factory("cexpr::toyMeanElecNr('toyElecNrYield*Enr',toyElecNrYield[4.0],Enr)");
   w->factory("PROD::NrQuantaModel(Poisson::photNr(nPhot,toyMeanPhotNr),Poisson::elecNr(nE,toyMeanElecNr))");

   // Declare P(r) and P(z) for WIMPs  - only need to do this once for all other uniform pdfs
   w->factory("CEXPR::UniformR('2*TMath::Pi()*r',r)");
   w->factory("Uniform::UniformZ(z)");

   // Multiply P(nPhot,nE|Enr)P(Enr)P(r)P(z)
   //w->factory("PROD::SigQuantaModel(NrQuantaModel|Enr, dRdE, UniformR, UniformZ)");
   //for now, skip position discrimination
   w->factory("PROD::SigQuantaModel(NrQuantaModel|Enr, dRdE)");

//    w->factory("cexpr::toyMeanPhe('toyLightYield*nPhot',toyLightYield[0.14],nPhot)");
//    w->factory("Poisson::toyPhe(nPhe[0,1200],toyMeanPhe)");
//   //could insert coincidence efficiency here
//    w->factory("cexpr::toyS1Sigma('sqrt(toyS1ResPar*toyPhe)/toyLightYield',toyS1ResPar[0.4],toyPhe,toyLightYield)");
//    w->factory("cexpr::toyS1Mean('toyPhe/toyLightYield',toyPhe,toyLightYield)");
//    w->factory("Gaussian::S1Res(S1,toyS1Mean,toyS1Sigma)");
   w->factory("cexpr::toyS1res('0.1+sqrt(nPhot/toyPheYield)',toyPheYield[0.14],nPhot)");
   w->factory("Gaussian::S1Res(S1,nPhot,toyS1res)");  

   w->factory("cexpr::toyS2res('0.1+sqrt(nE/toySEYield)',toySEYield[0.60],nE)");
   w->factory("Gaussian::S2Res(S2,nE,toyS2res)");


   
   //OK, can I marginalise over anything?
   TCanvas *c1 = new TCanvas("c1","c1");
   c1->Divide(2,1);
   c1->cd(1);
   RooArgSet quantaSet(Enr);
   RooProjectedPdf *SigModelIdeal = (RooProjectedPdf *) w->pdf("SigQuantaModel")->createProjection(quantaSet);
   SigModelIdeal->SetName("SigModelIdeal");
   w->import(*SigModelIdeal);
   SigModelIdeal->generateBinned(RooArgSet(nPhot,nE),10000)->createHistogram("nPhot,nE")->Draw("colz");
   c1->Print("quanta.eps");
   
   // Multiply to get full signal model
   w->factory("PROD::SigModel(S1Res|nPhot, S2Res|nE, SigModelIdeal)");
   c1->cd(2);
   RooArgSet simpleProjSet(nPhot,nE);
   RooProjectedPdf *SigModelSimpleProj = (RooProjectedPdf *) w->pdf("SigModel")->createProjection(simpleProjSet);
   SigModelSimpleProj->SetName("SigModelSimpleProj");
   w->import(*SigModelSimpleProj);
   w->pdf("SigModelSimpleProj")->generateBinned(RooArgSet(S1,S2),100)->createHistogram("S1,S2")->Draw("colz");
   c1->Print("observed.eps");
   //w->pdf("SigModel")->generateBinned(RooArgSet(S1,S2),100)->createHistogram("S1,S2")->Draw("colz");
   // Marginalize over E, nPhot, nE
//    RooProjectedPdf* SigModelProj = (RooProjectedPdf*)(*w->pdf("SigModel")).createProjection(*w->set("projNR"));
//    SigModelProj->SetName("SigModelProj");
//    TCanvas *c1 = new TCanvas("c1","c1");
//    SigModelProj->generateBinned(RooArgSet(S1,S2),1000)->createHistogram("S1,S2")->Draw("colz");
//    w->import(*SigModelProj);
   w->importClassCode();
   w->SetName("w");
   gDirectory->Add(w);
   w->writeToFile("w2.root");
}