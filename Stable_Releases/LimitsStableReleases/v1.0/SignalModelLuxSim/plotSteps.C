void importVars(RooWorkspace *w){
    // Import parameters
  gROOT->ProcessLine(".L signalParams.h"); // you might want to change this header (model params)
  
  w->import(targetA);
  w->import(targetZ);
  w->import(rho);
  w->import(vEsc);
  w->import(v0);
  w->import(vE);
  w->import(ff_c);
  w->import(ff_a);
  w->import(ff_s);
  w->import(ff_r0);
  w->import(idm_delta);
  w->import(halo_beta);

}

void plotSteps(){
  Double_t minRawS2_phe=200.;
  Int_t minEnergy_eV=3000;
  const double maxEnergy_keV=50.;
  Int_t minS1_phe=2;
  Int_t maxS1_phe=30;
  
//  bool regenerate=1;
  const int nm=20;
  double masses[nm]={
    6.0,
    7.0,
    8.0,
    9.0,
    10.0,
    12.0,
    14.0,
    17.0,
    21.0,
    27.0,
    33.0,
    40.0,
    50.0,
    60.0,
    70.0,
    80.0,
    100.0,
    300.0,
    1000.0,
    2000.
  };
  double eventsPerPbKgDay[nm];
  double integralAbove0PerPbKgDay[nm];
  double integralAboveCutofferPbKgDay[nm];
  
  gROOT->ProcessLine(".L signalParams.h");
  Double_t minEnergy_keV=double(minEnergy_eV)/1000.;
  gSystem->Load("libRooStats.so");
  RooWorkspace *w = new RooWorkspace("w");
  importVars(w);
  w->factory(Form("mWimp[%.5g]",50));
  gROOT->ProcessLine(".L ../RooWimpSpectrum.cxx+");//load a custom library
  w->importClassCode("RooWimpSpectrum", kTRUE);//import custom class into workspace
  w->factory(Form("Enr[%.1f,50.]",0.0));
  w->var("Enr")->setRange("fullRange",0,50);
  w->var("Enr")->setRange("cutoffRange",3,50);
  w->factory("WimpSpectrum::dRdE(Enr,mWimp, v0, vE, vEsc, rho, ff_c, ff_r0, ff_a, ff_s, idm_delta, halo_beta)");
  RooPlot *specframe = w->var("Enr")->frame(Title("recoil spectra per pb.kg.day"),Range(0.0,10.0));
  int colorCode[20];
  colorCode[0]=kRed-2;
  colorCode[1]=kRed;
  colorCode[13]=kRed+1;
  colorCode[19]=kBlack;
  for(int i=0;i<nm;i++){
    w->var("mWimp")->setVal(masses[i]);
    TFile wimpEventDensityFile(Form("wimpDensityFile_%05iGeV_S1_%uto%uphe.root",int(masses[i]),minS1_phe,maxS1_phe));
    TH2D *theHist = (TH2D*)wimpEventDensityFile.Get("wimpEventDensity");
    eventsPerPbKgDay[i]=theHist->Integral(1,theHist->GetNbinsX(),1,theHist->GetNbinsY());
    integralAbove0PerPbKgDay[i]=w->function("dRdE")->createIntegral(w->argSet("Enr"),"fullRange")->getVal();
    integralAboveCutofferPbKgDay[i]=w->function("dRdE")->createIntegral(w->argSet("Enr"),"cutoffRange")->getVal();
    if(i==0||i==1||i==13||i==19){
      w->function("dRdE")->plotOn(specframe,Normalization(integralAbove0PerPbKgDay[i],RooAbsReal::Raw),Name(Form("%.f GeV",masses[i])),LineColor(colorCode[i]));
    }
  }
  
  for(int i=0;i<nm;i++){
    cout<<Form("%.3g\t%.3g\t%.3g\t%.3g",masses[i],integralAbove0PerPbKgDay[i],integralAboveCutofferPbKgDay[i],eventsPerPbKgDay[i])<<endl;
  }
  specframe->Draw();

//   TGraph *gr = new TGraph(nm,masses,limitFromEmpty3000_cm2);
//   gr->SetTitle(Form("0-BG limit from net 3000kg.days, %u-%u phe S1c, modeling 3-50keV recoils, >%.1f S2r;WIMP mass / GeV.c^{-2};#sigma_{0}/cm^{-2}",minS1_phe,maxS1_phe,minRawS2_phe));
//   gr->Draw("alp");
//   gr->SetName(Form("simpleLimits__S1_%uto%uphe",minS1_phe,maxS1_phe));
//   gr->SaveAs(Form("simpleLimits__S1_%uto%uphe.root",minS1_phe,maxS1_phe));
//   gROOT->ProcessLine(".x drawOtherLimits.C");
}
