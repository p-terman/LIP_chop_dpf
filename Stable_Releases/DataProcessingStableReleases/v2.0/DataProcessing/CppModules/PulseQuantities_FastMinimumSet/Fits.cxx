#include <vector>
#include <cmath>
#include <iostream>
#include <TMinuit.h>
#include "Fits.h"

#define PI 3.14159265358979323   

using namespace std;

////////////////////////////////////////////////////////////////////////////////
// Global variables 
////////////////////////////////////////////////////////////////////////////////
FIter pulse_start;
FIter pulse_end;
TMinuit *minuit = new TMinuit(5);
double chisq = 0;
double gerror = 1;
double minimizer_small = 1e-10;
//______________________________________________________________________________
double GausFunc(double xx, double *par) {
	// [0] amplitude, [1] mean, [2] sigma 
	double x = (xx-par[1])/par[2];
	return par[0]*exp( -0.5*x*x);
}
//______________________________________________________________________________
double TwoExpFunc(double xx, double *par) {
	// [0] amplitude, [1] offset, [2] risetime, [3] falltime 
	double x = xx-par[1];
	double a1 = x/(fabs(par[3])+minimizer_small);
	double a2 = x/(fabs(par[2])+minimizer_small);
	if (x<0.0 or par[3]<=par[2]) {
		return 0;
	}else
		return par[0]*( exp(-a1) - exp(-a2) );
}
//______________________________________________________________________________
double TwoExpConvFunc(double xx, double *par) {
  // [0] amplitude, [1] offset, [2] risetime, [3] falltime, [4] sigma 
  double ss = par[4]*par[4];
  if (ss <= 0) ss = fabs(ss)+minimizer_small;
  double t = xx-par[1];
  double mu1 = t - ss/(fabs(par[2]) + minimizer_small);
  double mu2 = t - ss/(fabs(par[3]) + minimizer_small); 
  double f1 = exp(-0.5*t*t/ss + 0.5*mu1*mu1/ss);
  double f2 = exp(-0.5*t*t/ss + 0.5*mu2*mu2/ss);
  double a1 = (1-erf(-mu1/sqrt(ss*2)))/2.0;
  double a2 = (1-erf(-mu2/sqrt(ss*2)))/2.0;    
  return par[0]*(f2*a2 - f1*a1);
}
//______________________________________________________________________________
void Stats(double *stats) {
  // calculates basic stats quickly
  // [0] max
  // [1] area 
  // [2] mean 
  // [3] sigma 
  // [4] median 
  // [5] left sigma 
  // [6] right sigma 
  
  stats[0] = stats[1] = stats[2] = stats[3] = 0;
  stats[4] = stats[5] = stats[6] = 0;
  
  FVec cdf;
  double counter = 0;
  for (FIter it=pulse_start; it != pulse_end; it++) {
    if (stats[0] < *it) stats[0] = *it;
    stats[1] += *it;
    stats[2] += *it*counter;
    stats[3] += *it*counter*counter;
    counter += 1;   
    cdf.push_back(stats[1]);
  }
  if (stats[1] <= 0) stats[1] = fabs(stats[1])+minimizer_small;
  stats[2] /= stats[1];
  stats[3] -= stats[1]*stats[2]*stats[2];
  stats[3] = sqrt(fabs(stats[3])/(fabs(stats[1]-1)+minimizer_small));
  
  for (size_t i=0; i<cdf.size(); i++) {
    if (cdf[i]/stats[1] < 0.25) stats[5] = i;
    if (cdf[i]/stats[1] < 0.50) stats[4] = i;
    if (cdf[i]/stats[1] < 0.75) stats[6] = i;
  }
  stats[5] = fabs(stats[5]-stats[4])/0.6745;
  stats[6] = fabs(stats[6]-stats[4])/0.6745;
}
////////////////////////////////////////////////////////////////////////////////
// TMinuit Functions
////////////////////////////////////////////////////////////////////////////////
void FNC1(int &npar, double *grad, double &fval, double *par, int flag) {
  // This function calculates the Chisquare. It is passed to a minimizer 
  // (TMinuit) to minimize the Chisquare value. 
  // npar = number of parameters (defined in TMinuit)
  // grad = first derivitives of the function (optional). Not used. 
  // fval = Chisquare value 
  // par  = parameters. The first parameter (par[0]) is set to a constant and is
  //        used to choose which function to fit to, ExpFunc or GausFunc.
  // flag = TMinuit flag corresponding to its options 
  
  fval = 0;
  double theory = 0, error;
  size_t i = 0;
  for (FIter it=pulse_start; it != pulse_end; it++) {
    theory = TwoExpConvFunc(i, &par[0]);
    error = sqrt(theory*theory+gerror*gerror)+minimizer_small;
    fval += (theory - *it)*(theory - *it)/error;
    i++;
  }
  if (!fval) fval = 1e100;
  chisq = fval;
}
//______________________________________________________________________________
void FNC2(int &npar, double *grad, double &fval, double *par, int flag) {
  // This function calculates the Chisquare. It is passed to a minimizer 
  // (TMinuit) to minimize the Chisquare value. 
  // npar = number of parameters (defined in TMinuit)
  // grad = first derivitives of the function (optional). Not used. 
  // fval = Chisquare value 
  // par  = parameters. The first parameter (par[0]) is set to a constant and is
  //        used to choose which function to fit to, ExpFunc or GausFunc.
  // flag = TMinuit flag corresponding to its options 
  
  fval = 0;
  double theory = 0, error;
  size_t i = 0;
  for (FIter it=pulse_start; it != pulse_end; it++) {
    theory = GausFunc(i, &par[0]);
    error = sqrt(theory*theory+gerror*gerror)+minimizer_small;
    fval += (theory - *it)*(theory - *it)/error;
    i++;
  }
  if (!fval) fval = 1e100;
  chisq = fval;
}
////////////////////////////////////////////////////////////////////////////////
// Fitting 
////////////////////////////////////////////////////////////////////////////////
void Fit(FVec *wave, size_t start, size_t end, double perror, double *results) {
  // results [0] = exp_fit_amplitude_phe_per_sample 
  // results [1] = exp_fit_time_offset_samples - t0_height_samples
  // results [2] = exp_fit_tau_rise_samples 
  // results [3] = exp_fit_tau_fall_samples 
  // results [4] = exp_fit_chisq 
  // results [5] = exp_fit_dof 
  // results [6] = gaus_fit_amplitude_phe_per_sample
  // results [7] = gaus_fit_mu_samples - t0_height_samples
  // results [8] = gaus_fit_sigma_samples
  // results [9] = gaus_fit_chisq
  // results [10] = gaus_fit_dof
  // results [11] = convolution sigma 
  
  // Set values to 0
  results[0] = results[1] = results[2] = results[3] = results[4] = 0;
  results[5] = results[6] = results[7] = results[8] = results[9] = 0;
  results[10] = results[11] = 0;
  if (end <= start) {
    return;
  }
  
  // Set Globals
  gerror = perror;
  pulse_start = wave->begin()+start;
  pulse_end = wave->begin()+end+1;
  minuit->SetMaxIterations(50);
  minuit->SetPrintLevel(-1); // (-1) = quiet mode 
  
  // get guesses
  double stats[7];
  Stats(stats);
  double max = stats[0];
  double area = stats[1];
  double mean = stats[2];
  double sigma = stats[3];
  double median = stats[4];
  double lsigma = stats[5];
  double rsigma = stats[6];  
  double ampl = (max + area/sigma)/2; 
  double b0 = 0, b1 = int(end)-int(start);
  
  // 2 Exponential Fit 
  results[5] = int(end)-int(start)+1-5;  
  double error = 0; // used to get the error of the parameter
  if (results[5] >= 0) {
    minuit->SetFCN(&FNC1);
    minuit->DefineParameter(0, "Amplitude", ampl, 2, ampl/10, ampl*10);
    minuit->DefineParameter(1, "Offset", mean/2, 0.5, b0, b1);
    minuit->DefineParameter(2, "Rise", lsigma, .01, 1e-4, lsigma*10);
    minuit->DefineParameter(3, "Fall", rsigma, .01, 1e-3, rsigma*10);
    minuit->DefineParameter(4, "Sigma", 1, 0.001, 0.001, 5);
    minuit->FixParameter(4);
    minuit->Migrad();
    minuit->GetParameter(0, results[0], error);
    minuit->GetParameter(1, results[1], error);
    minuit->GetParameter(2, results[2], error);
    minuit->GetParameter(3, results[3], error);
    minuit->GetParameter(4, results[11], error);
    results[4] = chisq;
  }
  if (!isfinite(results[0])) results[0] = ampl;
  if (!isfinite(results[1])) results[1] = mean/2;
  if (!isfinite(results[2])) results[2] = lsigma;
  if (!isfinite(results[3])) results[3] = rsigma;
  if (!isfinite(results[11])) results[11] = 1;
  if (!isfinite(results[4])) results[4] = 0;
  
  // Gaussian Fit
  results[10] = int(end)-int(start)+1-3; 
  if (results[5] >= 0) {
    ampl = area/(2.51*sigma);
    minuit->SetFCN(&FNC2);
    minuit->DefineParameter(0, "Amplitude", ampl, 2., ampl/10, ampl*10);
    minuit->DefineParameter(1, "Mean", (mean+median)/2, 0.5, b0, b1);
    minuit->DefineParameter(2, "Sigma", sigma, .01, sigma/10, sigma*10);
    minuit->Migrad();
    minuit->GetParameter(0, results[6], error);
    minuit->GetParameter(1, results[7], error);
    minuit->GetParameter(2, results[8], error);
    results[9] = chisq;
  }
  if (!isfinite(results[6])) results[6] = area/(2.51*sigma);
  if (!isfinite(results[7])) results[7] = mean;
  if (!isfinite(results[8])) results[8] = sigma;
  if (!isfinite(results[9])) results[9] = 0;
}

