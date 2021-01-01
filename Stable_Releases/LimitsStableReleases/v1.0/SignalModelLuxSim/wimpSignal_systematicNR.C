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

double wimpSignal(Double_t mWimp, Double_t minRawS2_phe=200., Int_t minEnergy_eV=3000, Int_t minS1_phe=2, Int_t maxS1_phe=30){
  const double shiftInLog10S2=-0.07;
  
  const double maxEnergy_keV=50.;
  if (minEnergy_eV%250!=0){
    cout<<"sorry, minEnergy_eV needs to be n*250, n=1..200"<<endl;
    return 0.;
  }
  Double_t minEnergy_keV=double(minEnergy_eV)/1000.;
  gSystem->Load("libRooStats.so");
  RooWorkspace *w = new RooWorkspace("w");
  importVars(w);
  w->factory(Form("mWimp[%.5g]",mWimp));
  gROOT->ProcessLine(".L ../RooWimpSpectrum.cxx+");//load a custom library
  w->importClassCode("RooWimpSpectrum", kTRUE);//import custom class into workspace
  
  w->factory(Form("Enr[%.1f,50.]",minEnergy_keV));
  w->factory("WimpSpectrum::dRdE(Enr,mWimp, v0, vE, vEsc, rho, ff_c, ff_r0, ff_a, ff_s, idm_delta, halo_beta)");
  
  
//   RooPlot *eframe = w->var("Enr")->frame(Title("spectrum in Enr"));
//   eframe->SetName("eframe");
//   w->function("dRdE")->plotOn(eframe);
//   gDirectory->Add(w);
//   gDirectory->Add(eframe);
//   
//   RooArgSet enAS = w->argSet("Enr");
//   RooAbsReal *wimpSigIntergal = (RooAbsReal*) w->function("dRdE")->createRunningIntegral(enAS);
//   wimpSigIntergal->plotOn(eframe);
//   eframe->Draw();
  
  char* simRqPrefix="/media/hdv0/TestData/rq/luxsm_20130914T0100_cp05431/rootfiles/luxsm_20130914T0100_f";
  char* simRqSuffix = "_cp05431.rq.root";
  TFile wimpEventDensityFile(TString::Format("wimpDensityFile_NRUp%0.2f_%05iGeV_S1_%uto%uphe.root",shiftInLog10S2,int(mWimp),minS1_phe,maxS1_phe),"recreate");
  TH2D *wimpEventDensity = new TH2D("wimpEventDensity",TString::Format("%i GeV WIMP signal per pb.kg.day.bin;S1c / phe;log_{10}(#frac{S2b}{S1c})",int(mWimp)),2*(maxS1_phe-minS1_phe),minS1_phe,maxS1_phe,200,0.5,2.5); //you might want to change this (logS2/S1 range)
  wimpEventDensity->Sumw2();
  TCanvas *cWimps = new TCanvas("cWimps","cWimps");
  cWimps->Divide(2,2);
  TTree*  events;
  //cWimps.Print("cWimps.pdf[");
  Double_t primariesPerkeV = 1000.*4.; //1000 throws per file, 4 files per keV
  for(int e_eV=minEnergy_eV;e_eV<(1000.*maxEnergy_keV+1);e_eV+=250){ // you might want to change this (low energy cutoff)
    gDirectory->cd();
    TString filename = TString::Format("%s%09i%s",simRqPrefix,e_eV,simRqSuffix);
    TFile simRqFile(TString::Format("%s%09i%s",simRqPrefix,e_eV,simRqSuffix));
    events=(TTree*)simRqFile.Get("events");
    w->var("Enr")->setVal(double(e_eV/1000.));
    if(e_eV==minEnergy_eV){
      w->var("Enr")->setVal(double((e_eV+62.5)/1000.)); //make the first bin a half-width bin
    }
    double weighting=w->function("dRdE")->getVal()/primariesPerkeV;
    if(e_eV==minEnergy_eV){
      weighting *=0.5; //make the first bin a half-width bin
    }
    TH2D *theHist = (TH2D*)wimpEventDensityFile.Get("wimpEventDensity");
    //cout<<theHist->GetName()<<endl;
    //cout<<e_eV<<":"<<weighting<<endl;
    events->SetAlias("S2b","MaxIf$(xyz_corrected_pulse_area_bot_phe,pulse_classification==2)");
    events->SetAlias("S1c","MaxIf$(xyz_corrected_pulse_area_all_phe,pulse_classification==1)");
    events->SetAlias("S2r","MaxIf$(pulse_area_phe,pulse_classification==2)");
    wimpEventDensityFile.cd();
    events->Draw(Form("log10(S2b/S1c)+%0.2f:S1c>>+wimpEventDensity",shiftInLog10S2),TString::Format("(S2r>%.1f&&golden[0])*%G",minRawS2_phe/pow(10,shiftInLog10S2),weighting),"goff"); //****apply shift in S2|S1 here****
    //wimpEventDensity->Draw();
    //cWimps.Print("cWimps.pdf");
  }
  //cWimps.Print("cWimps.pdf]");
  cWimps->SetLogz();
  wimpEventDensity->GetYaxis()->SetTitleOffset(1.5);
  wimpEventDensity->GetZaxis()->SetRangeUser(1E-5*wimpEventDensity->GetMaximum(),1.5*wimpEventDensity->GetMaximum());
  Double_t events_per_pbkgday = wimpEventDensity->Integral(1,wimpEventDensity->GetNbinsX(),1,wimpEventDensity->GetNbinsY());
  wimpEventDensity->SetTitle(Form("#splitline{%s}{(%.2E events per pb.kg.day,%u-%uphe&&%.2f-%.2fkeV&&S2raw>%.1f)}",wimpEventDensity->GetTitle(),events_per_pbkgday,minS1_phe,maxS1_phe,minEnergy_keV,maxEnergy_keV,minRawS2_phe));
//   wimpEventDensity->FitSlicesY();
//   wimpEventDensity_1->Draw("same");
  //cWimps->cd(1)->SetLogz();
  wimpEventDensity->Draw("colz");
  wimpEventDensity->ProfileX();
  wimpEventDensity_pfx->Draw("same");
  //cWimps->cd(3)->SetLogy();
  //wimpEventDensity->ProjectionX()->Fit("expo");
  cout<<2.3/(3000*events_per_pbkgday)<<endl;
  cWimps->Print(TString::Format("wimpDensity_NRUp%0.2f_%05iGeV_S1_%uto%uphe.pdf",shiftInLog10S2,int(mWimp),minS1_phe,maxS1_phe));
  //fill the output file with convenient formats
  RooWorkspace *wSignal = new RooWorkspace("wSignal");
  wSignal->factory(Form("S1[%.5g,%.5g]",double(minS1_phe),double(maxS1_phe)));
  wSignal->factory("logS2S1[0.5,2.5]");
  RooDataHist rdh("rdh","rdh",wSignal->argSet("S1,logS2S1"),wimpEventDensity);
  RooHistPdf signalHistPdf("signalHistPdf","signalHistPdf",wSignal->argSet("S1,logS2S1"),rdh);
  wSignal->import(signalHistPdf);
  
  //save the ranges used to generate the signal model
  wSignal->factory(Form("minEnergy_keV[%.5g]",double(minEnergy_keV)));
  wSignal->factory(Form("maxEnergy_keV[%.5g]",double(maxEnergy_keV)));
  wSignal->factory(Form("minS1_phe[%.5g]",double(minS1_phe)));
  wSignal->factory(Form("maxS1_phe[%.5g]",double(maxS1_phe)));
  
  //save the WIMP Mass and the scaling from xs to events in histPdf
  wSignal->factory(Form("mWimp[%.5g]",mWimp));
  wSignal->factory(Form("events_per_pbkgday[%.5g]",events_per_pbkgday));
  wimpEventDensityFile.Add(wSignal);
  wimpEventDensityFile.Write();
  return events_per_pbkgday;
  //gDirectory->Add(wimpEventDensity);
}
