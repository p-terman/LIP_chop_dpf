// ProfileLikelihood.h : a header file to define all the RooRealVars for ProfileLikelihood.C
#include "RooRealVar.h"

using namespace RooFit;


//We tend to use xenon...
RooRealVar targetZ("targetZ","Protons per target nucleus",54);
RooRealVar targetA("targetA","Nucleons per target nucleus ",131);

//Halo parameters
RooRealVar rho("rho","WIMP density [GeV.c^-{2}.cm^{-3}]",0.3); // density of WIMPs in LUX [GeV.c^-2.cm^-3]
RooRealVar vE("vE","v_{earth} [km/s]",245.); // speed of the Earth in the DM halo frame [km/s] //245 corresponds to mean of LUX run
RooRealVar vEsc("vEsc","v_{esc} [km/s]",544); // local escape velocity [km/s]
RooRealVar v0("v0","v_0 of halo",220); // sqrt(2/3) times the velocity dispersion of WIMPs in halo frame

// Form factor parameters
RooRealVar ff_c("ff_c","Nuclear radius parameter",5.6384); // nuclear radius parameter [fm]
RooRealVar ff_r0("ff_r0","Nuclear radius parameter",6.0945);  // nuclear radius parameter - not currently used in f.f. [fm]
RooRealVar ff_a("ff_a","Nuclear radius parameter",0.523); // nuclear radius parameter [fm]
RooRealVar ff_s("ff_s","Skin depth parameter",1.0); // nuclear skin thickness[fm]
RooRealVar idm_delta("idm_delta","Inelastic mass splitting",0.); // inelastic mass splitting [keV]
RooRealVar halo_beta("halo_beta","Halo beta parameter",0.); // currently assumed 0 in WIMP velocity integral

