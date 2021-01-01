#ifndef __PHET_H__
#define __PHET_H__

////////////////////////////////////////////////////////////////////////////////
// Generic Functions
////////////////////////////////////////////////////////////////////////////////

double HorizontalFunction(double *x, double *p);
double VerticalFunction(double *x, double *p);
double LineFunction(double *x, double *p);
double QuadraticFunction(double *x, double *p);
double SigmoidFunction(double *x, double *p);

////////////////////////////////////////////////////////////////////////////////
// Boundaries
////////////////////////////////////////////////////////////////////////////////

double S1LowerAreaVsWAW(double area);
double S1UpperAreaVsWAW(double area);
double S2LowerAreaVsWAW(double area);
double S2UpperAreaVsWAW(double area);

double S1S2LowerTBvsMPA(double mpa);
double S1S2UpperTBvsMPA(double mpa);

double S1S2UpperMPAvsPA(double area);

double S2LowerFilterDiffVsArea(double area);
double S1UpperFilterDiffVsArea(double area);

////////////////////////////////////////////////////////////////////////////////
// Pulse Classification 
////////////////////////////////////////////////////////////////////////////////
int PulseClassify(double *values, double *parameters);

#endif
