#ifndef __USEFUL_FUNCTIONS_H__ 
#define __USEFUL_FUNCTIONS_H__ 1
#include <vector>
#include <cmath>

using std::sqrt;
using std::fabs;
//using std::erf;
//using std::tgamma;
using std::exp;
using std::vector;

typedef std::vector<float> FVec;
typedef std::vector<float>::iterator FIter;

double small = 1e-6;

//______________________________________________________________________________
void LinearFit(FIter start, FIter end, double &a0, double &a1) {  
  double n, x, x2, y, xy, N;
  n = x = x2 = y = xy = 0;
  for (FIter it=start; it != end; ++it) {
    y += *it; xy += (*it)*n;
    x += n; x2 += n*n;
    n += 1; 
  }
  if (n < 2) return;
  N = n*x2 -x*x;
  if (N == 0) N = small;
  a0 = (x2*y - x*xy)/N;
  a1 = (n*xy - x*y)/N;
}
//______________________________________________________________________________
void QuadraticFit(FIter start, FIter end, double &a0, double &a1, double &a2) {  
  double n, x, x2, x3, x4, y, xy, x2y, N;
  n = x = x2 = x3 = x4 = y = xy = x2y = 0;
  for (FIter it=start; it != end; ++it) {
    y += *it; xy += (*it)*n; x2y += (*it)*n*n;
    x += n; x2 += n*n; x3 += n*n*n; x4 += n*n*n*n;
    n += 1; 
  }
  if (n < 3) return;
  N = n*(x2*x4-x3*x3) + x*(x2*x3-x*x4) + x2*(x*x3-x2*x2);
  if (N == 0) N = small;
  a0 = ((x2*x4-x3*x3)*y + (x2*x3-x4*x)*xy + (x*x3-x2*x2)*x2y)/N;
  a1 = ((x2*x3-x4*x)*y + (n*x4-x2*x2)*xy + (x*x2-n*x3)*x2y)/N;
  a2 = ((x*x3-x2*x2)*y + (x*x2-n*x3)*xy + (n*x2-x*x)*x2y)/N;
}
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
	double a1 = x/(fabs(par[3])+small);
	double a2 = x/(fabs(par[2])+small);
	if (x<0.0 or par[3]<=par[2]) {
		return 0;
	}else
		return par[0]*( exp(-a1) - exp(-a2) );
}
//______________________________________________________________________________
double TwoExpConvFunc(double xx, double *par) {
  // [0] amplitude, [1] offset, [2] risetime, [3] falltime, [4] sigma 
  double ss = par[4]*par[4];
  if (ss <= 0) ss = fabs(ss)+small;
  double t = xx-par[1];
  double mu1 = t - ss/(fabs(par[2]) + small);
  double mu2 = t - ss/(fabs(par[3]) + small); 
  double f1 = exp(-0.5*t*t/ss + 0.5*mu1*mu1/ss);
  double f2 = exp(-0.5*t*t/ss + 0.5*mu2*mu2/ss);
  double a1 = (1-erf(-mu1/sqrt(ss*2)))/2.0;
  double a2 = (1-erf(-mu2/sqrt(ss*2)))/2.0;    
  return par[0]*(f2*a2 - f1*a1);
}
//______________________________________________________________________________
double StudentT(double m1, double s1, int n1, double m2, double s2, int n2) {
  // This performs a Student T Test on two distributions with some mean and std.
  // It returns the z value which is equvilent as t at large ndof. 
  if (n1 == 1 or n2 == 1) return 0;  
  double t, ss, ndof, denom;
  // Welch's' t-test 
  ss = sqrt(fabs(s1*s1/n1 + s2*s2/n2));
  ndof = ss*ss/( (s1*s1/n1)*(s1*s1/n1)/(n1-1) +  (s2*s2/n2)*(s2*s2/n2)/(n2-1) );
  t = fabs(m1-m2)/ss;
  // Convert t distributions to a z distribution 
  // Taken from root.cern.ch TMath::Student() function
  if (ndof < 1) return 0;
  denom = sqrt(ndof*M_PI)*tgamma(0.5*ndof)*pow(1+t*t/ndof, 0.5*ndof+0.5);
  return tgamma(0.5*ndof+0.5)/denom;
}



#endif
