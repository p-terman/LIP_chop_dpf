#include <cmath>
#include <iostream>

#include "PHET.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
// Generic Functions
////////////////////////////////////////////////////////////////////////////////
double HorizontalFunction(double *x, double *p) { return p[1]; }
//______________________________________________________________________________
double VerticalFunction(double *x, double *p) { return p[0]; }
//______________________________________________________________________________
double LineFunction(double *x, double *p) {
	if (p[2]-p[0] == 0) return p[1];
	double a0 = p[3] * (*x-p[0])/(p[2]-p[0]);
	double a1 = p[1] * (*x-p[2])/(p[0]-p[2]);
	return a0 + a1;
}
//______________________________________________________________________________
double QuadraticFunction(double *x, double *p) {
	double d = (p[0]-p[2])*(p[0]-p[4])*(p[2]-p[4]);
	if (d == 0) {
		if (*x < p[2]) return LineFunction(x, &p[0]);
		else return LineFunction(x, &p[2]);
	}
	double a0 = p[1] * (*x-p[2])/(p[0]-p[2]) * (*x-p[4])/(p[0]-p[4]);
	double a1 = p[3] * (*x-p[0])/(p[2]-p[0]) * (*x-p[4])/(p[2]-p[4]);
	double a2 = p[5] * (*x-p[2])/(p[4]-p[2]) * (*x-p[0])/(p[4]-p[0]);
	return a0 + a1 + a2;
}
//______________________________________________________________________________
double SigmoidFunction(double *x, double *p) {
	double yMid   = (p[3]+p[1])/2.; //(Y[1]+Y[0])/2.;
	double yScale = (p[3]-p[1])/2.; // (Y[1]-Y[0])/2.;
	double xMid   = (p[2]+p[0])/2.; // (X[1]+X[0])/2.;
	double xScale = (p[2]-p[0])/5.; // (X[1]-X[0])/5.;
	return yMid + yScale*tanh((*x-xMid)/xScale);
}

////////////////////////////////////////////////////////////////////////////////
// Boundaries
////////////////////////////////////////////////////////////////////////////////

//------++++++------++++++------++++++------++++++------++++++------++++++------
// Settings for pulse_area_phe vs. WAW
//------++++++------++++++------++++++------++++++------++++++------++++++------
double S1LowerAreaVsWAW(double area) {
  //  points    =   { x, y }
  double pars[] = {log10(1.0), log10(0.25), 
                   log10(1e5), log10(2.0)};
  area = log10(area);
  // functions 
  if (area < pars[0]) return pow(10., HorizontalFunction(&area, &pars[0]));
  if (area < pars[2]) return pow(10., SigmoidFunction(&area, &pars[0]));
  return pow(10., HorizontalFunction(&area, &pars[2]));
}
//______________________________________________________________________________
double S1UpperAreaVsWAW(double area) {
  //  points    =   { x, y }
  double pars[] = {log10(1.00), log10(25.0), 
                   log10(3.00), log10(20.0),
                   log10(10.0), log10(10.0),
                   log10(100.), log10(20.0),
                   log10(3.e6), log10(20.0)};
  area = log10(area);
  // functions 
  if (area < pars[0]) return pow(10., HorizontalFunction(&area, &pars[0]));
  if (area < pars[2]) return pow(10., LineFunction(&area, &pars[0]));
  if (area < pars[6]) return pow(10., QuadraticFunction(&area, &pars[2]));
  if (area < pars[8]) return pow(10., LineFunction(&area, &pars[6]));
  return pow(10., HorizontalFunction(&area, &pars[8]));
}
//______________________________________________________________________________
double S2LowerAreaVsWAW(double area) {
  //  points    =   { x, y }
  double pars[] = {log10(1.00), log10(20.0), 
                   log10(3.00), log10(15.0),
                   log10(10.0), log10(7.00),
                   log10(100.), log10(15.0),
                   log10(3.e6), log10(15.0)};
  area = log10(area);
  // functions 
  if (area < pars[0]) return pow(10., HorizontalFunction(&area, &pars[0]));
  if (area < pars[2]) return pow(10., LineFunction(&area, &pars[0]));
  if (area < pars[6]) return pow(10., QuadraticFunction(&area, &pars[2]));
  if (area < pars[8]) return pow(10., LineFunction(&area, &pars[6]));
  return pow(10., HorizontalFunction(&area, &pars[8]));
}
//______________________________________________________________________________
double S2UpperAreaVsWAW(double area) {
  //  points    =   { x, y }
  double pars[] = {log10(10.0), log10(200.), 
                   log10(1.e4), log10(200.),
                   log10(1.e5), log10(5.e3)};
  area = log10(area);
  // functions 
  if (area < pars[0]) return pow(10., HorizontalFunction(&area, &pars[0]));
  if (area < pars[2]) return pow(10., LineFunction(&area, &pars[0]));
  if (area < pars[4]) return pow(10., LineFunction(&area, &pars[2]));
  return pow(10., HorizontalFunction(&area, &pars[4]));
}

//------++++++------++++++------++++++------++++++------++++++------++++++------
// Settings for tb asymmetry vs. max_peak_area/pulse_area
//------++++++------++++++------++++++------++++++------++++++------++++++------
double S1S2LowerTBvsMPA(double mpa) {
  //  points    =   { x, y }
  double pars[] = {0.5,-1.0, 
                   1.0,-0.7, 
                   1.4,-0.6};
  if (mpa < pars[0]) return HorizontalFunction(&mpa, &pars[0]);
  if (mpa < pars[4]) return QuadraticFunction(&mpa, &pars[0]);
  return HorizontalFunction(&mpa, &pars[4]);
}
//______________________________________________________________________________
double S1S2UpperTBvsMPA(double mpa) {
  //  points    =   { x, y }
  double pars[] = {0.5,1.0, 
                   1.0,0.7, 
                   1.4,0.6};
  if (mpa < pars[0]) return HorizontalFunction(&mpa, &pars[0]);
  if (mpa < pars[4]) return QuadraticFunction(&mpa, &pars[0]);
  return HorizontalFunction(&mpa, &pars[4]);
}

//------++++++------++++++------++++++------++++++------++++++------++++++------
// Settings for max_peak_area/pulse_area vs pulse_area
//------++++++------++++++------++++++------++++++------++++++------++++++------
double S1S2UpperMPAvsPA(double area) {
  //  points    =   { x, y }
  double pars[] = {log10(0.7), 0.65, 
                   log10(2.0), 0.70, 
                   log10(5.0) ,0.80, 
                   log10(20.), 0.85, 
                   log10(30.), 0.95};
  area = log10(area);
  if (area < pars[0]) return HorizontalFunction(&area, &pars[0]);
  if (area < pars[4]) return QuadraticFunction(&area, &pars[0]);
  if (area < pars[6]) return QuadraticFunction(&area, &pars[4]);
  return HorizontalFunction(&area, &pars[6]);
}

//------++++++------++++++------++++++------++++++------++++++------++++++------
// Settings for s2 filter area difference ratio vs. pulse area phase space
//------++++++------++++++------++++++------++++++------++++++------++++++------
double S2LowerFilterDiffVsArea(double area) {
  //  points    =   { x, y }
  double pars[] = {log10(5.0), log10(3.e-4), 
                   log10(3.e3),log10(0.3)};
  area = log10(area);
  if (area < pars[0]) return pow(10., HorizontalFunction(&area, &pars[0]));
  if (area < pars[2]) return pow(10., SigmoidFunction(&area, &pars[0]));
  return pow(10., HorizontalFunction(&area, &pars[2]));
}
//______________________________________________________________________________
double S1UpperFilterDiffVsArea(double area) {
  //  points    =   { x, y }
  double pars[] = { log10(1.0), log10(1.0), 
                    log10(5.0), log10(0.3), 
                    log10(15.), log10(4.e-2), 
                    log10(60.), log10(0.14), 
                    log10(1.e6), log10(0.14)};
  area = log10(area);
  if (area < pars[0]) return pow(10., HorizontalFunction(&area, &pars[0]));
  if (area < pars[2]) return pow(10., LineFunction(&area, &pars[0]));
  if (area < pars[6]) return pow(10., QuadraticFunction(&area, &pars[2]));
  if (area < pars[8]) return pow(10., LineFunction(&area, &pars[6]));
  return pow(10., HorizontalFunction(&area, &pars[8]));
}



////////////////////////////////////////////////////////////////////////////////
// Pulse Classification 
////////////////////////////////////////////////////////////////////////////////
int PulseClassify(double *values, double *parameters) {
  /// Pulse types: 0 empty, 1 S1, 2 S2, 3 SP, 4 RE, 5 Other, 6 BL, 7 Ch, 8 ET
  
  double area    = values[0];
  double waw     = values[1];
  double pRatio  = values[2];
  double coin    = values[3];
  double tbAsym  = values[4];
  double filter  = values[5];
  double filter_area_ratio = filter/area;

  // What the parameters are 
  double s1MinArea        = parameters[0]; // 1.0
  double s1MinCoincidence = parameters[1]; // 2.0
  double s2MinArea        = parameters[2]; // 50.0
  double spMaxArea        = parameters[3]; // 2.0
  
  bool s1like=true, s2like=true, tooWide=true;
  
  // S1 Cuts
  s1like = s1like && waw    > S1LowerAreaVsWAW(area);
  s1like = s1like && waw    < S1UpperAreaVsWAW(area);
  s1like = s1like && tbAsym > S1S2LowerTBvsMPA(pRatio);
  s1like = s1like && tbAsym < S1S2UpperTBvsMPA(pRatio);
  s1like = s1like && pRatio < S1S2UpperMPAvsPA(area);
  s1like = s1like && filter_area_ratio < S1UpperFilterDiffVsArea(area);

  // S2 Cuts
  s2like = s2like && waw    > S2LowerAreaVsWAW(area);
  s2like = s2like && waw    < S2UpperAreaVsWAW(area);
  s2like = s2like && tbAsym > S1S2LowerTBvsMPA(pRatio);
  s2like = s2like && tbAsym < S1S2UpperTBvsMPA(pRatio);
  s2like = s2like && pRatio < S1S2UpperMPAvsPA(area);
	s2like = s2like && filter_area_ratio > S2LowerFilterDiffVsArea(area);
  
  tooWide = waw >= S2UpperAreaVsWAW(area);
  
  // go through the KRraaaZY logic  
  if (coin < s1MinCoincidence) {
    if (s1like && area>0) {
      if (area>spMaxArea) return 3; // Single Photon
      return 7; // Cherenkov
    }
    if (area>s1MinArea) {
      return 5; // Don't know, but most likely junk
    }
    return 6; // Baseline 
  }
  
  if (s1like && s2like) return 5; // Can't decide
  
  if (s1like) {
    if (area>=s1MinArea) return 1; // S1
    return 3; // Single Photon
  }
  
  if (tooWide) {
    if (area < s1MinArea) return 6; // digitized baseline 
    if (area > s2MinArea) return 8; // e-train 
    return 4; // Single Electron 
  }
	
	if (s2like) {
    if (area>=s2MinArea) return 2; // S2
    return 4; // Single Electron
  }
  
  // Something else entirely
  return 5;
}
