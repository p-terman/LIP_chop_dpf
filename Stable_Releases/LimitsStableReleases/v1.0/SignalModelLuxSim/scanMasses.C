{
  Double_t minRawS2_phe=200.;
  Int_t minEnergy_eV=3000;
  Int_t minS1_phe=2;
  Int_t maxS1_phe=30;
  
  bool regenerate=1;
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
  double limitFromEmpty3000_cm2[nm];
  if(regenerate){
    gROOT->ProcessLine(".L wimpSignal_systematicNR.C");
    //gROOT->ProcessLine(".L wimpSignal.C");
    for(int i=0;i<nm;i++){
      eventsPerPbKgDay[i]=wimpSignal(masses[i],minRawS2_phe,minEnergy_eV,minS1_phe,maxS1_phe);
      limitFromEmpty3000_cm2[i]=2.3E-36/(3000*eventsPerPbKgDay[i]);
    }
  }
  else{
    for(int i=0;i<nm;i++){
      TFile wimpEventDensityFile(Form("wimpDensityFile_%05iGeV_S1_%uto%uphe.root",int(masses[i]),minS1_phe,maxS1_phe));
      TH2D *theHist = (TH2D*)wimpEventDensityFile.Get("wimpEventDensity");
      eventsPerPbKgDay[i]=theHist->Integral(1,theHist->GetNbinsX(),1,theHist->GetNbinsY());
      limitFromEmpty3000_cm2[i]=2.3E-36/(3000*eventsPerPbKgDay[i]);
      cout<<masses[i]<<"\t"<<eventsPerPbKgDay[i]<<endl;
    }
  }
  TGraph *gr = new TGraph(nm,masses,limitFromEmpty3000_cm2);
  gr->SetTitle(Form("0-BG limit from net 3000kg.days, %u-%u phe S1c, modeling 3-50keV recoils, >%.1f S2r;WIMP mass / GeV.c^{-2};#sigma_{0}/cm^{-2}",minS1_phe,maxS1_phe,minRawS2_phe));
  gr->Draw("alp");
  gr->SetName(Form("simpleLimits__S1_%uto%uphe",minS1_phe,maxS1_phe));
  gr->SaveAs(Form("simpleLimits__S1_%uto%uphe.root",minS1_phe,maxS1_phe));
  gROOT->ProcessLine(".x drawOtherLimits.C");
}
