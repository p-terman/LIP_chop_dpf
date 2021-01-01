#include "SplitGaussian.h"
#include <complex>

using namespace std;

double sign(double x) {return (x<0)? -1.: 1.;}

SplitGaussian::SplitGaussian() :
  amplitude(0),
  mode(0),
  lsigma(1),
  rsigma(1),
  area(0),
  mean(0),
  sigma(1),
  skew(0)
{ 
  c1 = sqrt(2./M_PI);
  c2 = 1.-2./M_PI;
  c3 = 4./M_PI-1.;
}
//______________________________________________________________________________
SplitGaussian::SplitGaussian(double a, double m, double l, double r) :
  amplitude(a),
  mode(m),
  lsigma(l),
  rsigma(r),
  area(0),
  mean(0),
  sigma(1),
  skew(0)
{ 
  c1 = sqrt(2./M_PI);
  c2 = 1.-2./M_PI;
  c3 = 4./M_PI-1.;
  Derive();
}
//______________________________________________________________________________
SplitGaussian SplitGaussian::operator+(const SplitGaussian &g) const {
  SplitGaussian gg;
  
  // useful varibles
  double ss, ss1, ss2, x2l, x2r, x3l, x3r;
	ss1 = sigma*sigma;
	x2l = ss1*area + pow(mean, 2)*area; // sum (x^2)
	x3l = skew*area + 3*mean*x2l - 2*area*pow(mean,3); // sum (x^3)
	ss2 = g.sigma*g.sigma;
	x2r = ss2*g.area + pow(g.mean, 2)*g.area;
	x3r = g.skew*g.area + 3*g.mean*x2r - 2*g.area*pow(g.mean,3);
  
  // Solve for the area, mean, sigma, and skewness
  gg.area = area + g.area;
	gg.mean = (area*mean + g.area*g.mean)/gg.area;
	ss = (x2l+x2r - gg.area*pow(gg.mean,2))/gg.area;
	gg.sigma = sqrt(ss);
	gg.skew = ((x3l+x3r) - 3*gg.mean*(x2l+x2r) + 2*gg.area*pow(gg.mean,3))/gg.area;
	
  gg.InverseDerive();
  //gg.Derive();
  return gg;
}
//______________________________________________________________________________
SplitGaussian& SplitGaussian::operator+=(const SplitGaussian &rhs) {
  *this = *this + rhs;
  return *this;
}
//______________________________________________________________________________
double SplitGaussian::operator[](double x) const {
  double z;
  if (x<=mode) z = (x-mode)/lsigma;
  else         z = (x-mode)/rsigma;
  return amplitude * exp(-0.5*z*z);
}
//______________________________________________________________________________
double SplitGaussian::Area(double b1, double b2) const {
  // fix boundaries
  if (b2<b1) {
    double tmp = b1;
    b1 = b2;
    b2 = tmp;
  }
  return this->CF(b2)-this->CF(b1);
}
//______________________________________________________________________________
double SplitGaussian::BhattacharyyaDistance(SplitGaussian &g) const {
  // calculates the Bhattacharyya Distance
  double dm, ss1, ss2, part1, part2;
  dm = fabs(mode - g.mode);
  if (mode<g.mode) {
    ss1 = rsigma*rsigma;
    if (ss1 <= 0) ss1 = sigma*sigma;
    ss2 = g.lsigma*g.lsigma;
    if (ss2 <= 0) ss2 = g.sigma*g.sigma;
  }else {
    ss1 = lsigma*lsigma;
    if (ss1 <= 0) ss1 = sigma*sigma;
    ss2 = g.rsigma*g.rsigma;
    if (ss2 <= 0) ss2 = g.sigma*g.sigma;
  }
  part1 = 0.25*dm*dm / (ss1+ss2);
  part2 = 0.25*(ss1/ss2 + ss2/ss1 + 2.);
  return part1 + 0.25*log(part2);
}
//______________________________________________________________________________
double SplitGaussian::BhattacharyyaDistance2(SplitGaussian &g) const {
  // calculates the Bhattacharyya Distance: doesn't use the left or right sigma
  double dm, ss1, ss2, part1, part2;
  dm = fabs(mean - g.mean);
  ss1 = sigma*sigma;
  ss2 = g.sigma*g.sigma;
  part1 = 0.25*dm*dm / (ss1+ss2);
  part2 = 0.25*(ss1/ss2 + ss2/ss1 + 2.);
  return part1 + 0.25*log(part2);
}
//______________________________________________________________________________
double SplitGaussian::CF(double x) const {
  double c=sqrt(0.5);
  if (x < mode) return area*0.5*(1.+erf(c*(x-mode)/lsigma));
	return area*0.5*(1.+erf(c*(x-mode)/rsigma));
}
//______________________________________________________________________________
void SplitGaussian::Derive() {
  area  = amplitude*(lsigma+rsigma)/c1;
  mean  = c1*(rsigma-lsigma) + mode;
  sigma = sqrt(c2*pow(rsigma-lsigma,2) + rsigma*lsigma);
  skew  = c1*(rsigma-lsigma)*(c3*pow(rsigma-lsigma,2) + rsigma*lsigma);
}
//______________________________________________________________________________
double SplitGaussian::Distance(double x) const {
  if (x < mode) return fabs(x-mode)/lsigma;
  return fabs(x-mode)/rsigma;
}
//______________________________________________________________________________
double SplitGaussian::EarthsMoversDistance(SplitGaussian &g, double b1, double b2) const {
  if (b2<b1) {
    double tmp = b1;
    b1 = b2;
    b2 = tmp;
  }
  SplitGaussian gg = *this + g;
  gg = gg.Sub(b1, b2);
  double norm_p = gg.Area(b1, b2) + 1.e-4;
  double norm_q = g.Area(b1, b2) + 1.e-4;
  double P = (b2-gg.mode) * gg.CF(b2)+ gg[b2]
            -((b1-gg.mode) * gg.CF(b1) + gg[b1]);
  double Q = (b2-g.mode) * g.CF(b2) + g[b2]
            -((b1-g.mode) * g.CF(b1) + g[b1]);
  double edm = fabs(P/norm_p - Q/norm_q);
  //cout << P << " " << Q << " " << norm_p << " " << norm_q << " " << edm << endl;
  //cout << "\t";
  //gg.Print(4);
  return edm;
}
//______________________________________________________________________________
void SplitGaussian::InverseDerive() {
  
  double ss = sigma*sigma;
  double c, d, sl, sr, s, min=-1, delta, dd;
	complex<double> C, x1, x2, x3;
	c = ss/(3.*c1*c1-2.);
	d = -skew/c1/(3.*c1*c1-2.);
	complex<double> del0(-3.*c, 0);
	complex<double> del1(27.*d, 0);
	complex<double> u1(1, 0);
	complex<double> u2(-0.5, sqrt(3.)/2.);
	complex<double> u3(-0.5, -sqrt(3.)/2.);
	C = pow(0.5*(del1+sqrt(del1*del1-4.*del0*del0*del0)), 1./3.);
	x1 = -1./3.*(u1*C + del0/u1/C);
	x2 = -1./3.*(u2*C + del0/u2/C);
	x3 = -1./3.*(u3*C + del0/u3/C);
    
  dd = x1.real();
  sl = sqrt(fabs((1.-4*c2)*dd*dd+4*ss))-dd/2.;
  sr = sqrt(fabs((1.-4*c2)*dd*dd+4*ss))+dd/2.;
  s = c2*pow(sr-sl,2) + sr*sl;
  if (min > fabs(s-ss) || min == -1) {
    min = fabs(s-ss);
    delta = x1.real();
  }
  dd = x2.real();
  sl = fabs(sqrt(fabs((1.-4*c2)*dd*dd+4*ss))-dd)/2.;
  sr = fabs(sqrt(fabs((1.-4*c2)*dd*dd+4*ss))+dd)/2.;
  s = c2*pow(sr-sl,2) + sr*sl;
  if (min > fabs(s-ss) || min == -1) {
    min = fabs(s-ss);
    delta = x2.real();
  }
  dd = x3.real();
  sl = fabs(sqrt(fabs((1.-4*c2)*dd*dd+4*ss))-dd)/2.;
  sr = fabs(sqrt(fabs((1.-4*c2)*dd*dd+4*ss))+dd)/2.;
  s = c2*pow(sr-sl,2) + sr*sl;
  if (min > fabs(s-ss) || min == -1) {
    min = fabs(s-ss);
    delta = x3.real();
  }
	
	// Solve for the mode, l/r sigma, amplitude
	mode = mean - c1*delta;
	lsigma = fabs(sqrt(fabs((1.-4*c2)*delta*delta+4*ss))-delta)/2.;
	rsigma = fabs(sqrt(fabs((1.-4*c2)*delta*delta+4*ss))+delta)/2.;
	amplitude = area*c1/(rsigma+lsigma);
}
//______________________________________________________________________________
void SplitGaussian::Fit(vector<float>::iterator s, vector<float>::iterator e, double offset) {

  // useful varibles
  int l = int(e-s), max_idx=0;
  vector<double> xxcf(l+1,0), xcf(l+1,0), cf(l+1,0);
  double w, x=0, max=0;
  vector<float>::iterator it = s;
  for (; it != e; it++) {
    w = (*it > 0) ? *it : 0;
    xxcf[x+1] = x*x*w + xxcf[x];
    xcf[x+1]  = x*w + xcf[x];
    cf[x+1]   = w + cf[x];
    x += 1;
    if (max < *it) {
      max = *it;
      max_idx = x;
    }
  }
  
  if (cf[l] < 0.1) {
    area = cf[l];
    mode = max_idx + offset;
    mean = max_idx + offset;
    lsigma = rsigma = sigma = 0.1;
    skew = 0;
    amplitude = area*c1/(rsigma+lsigma);
    return;
  }
  
  // Find the maximum likelihood
  double L1, L2, L, minL;
  area = cf[l];
  mean = xcf[l]/area;
  it = s;
  sigma = skew = 0;
  minL = 0;
  for (int i=1; i<l-1; i++) {
    L1 = xxcf[i] - 2.*i*xcf[i] + i*i*cf[i];
    L2 = (xxcf[l]-xxcf[i+1]) - 2.*i*(xcf[l]-xcf[i+1]) + i*i*(cf[l]-cf[i+1]);
    L  = sign(L1)*pow(fabs(L1), 1./3.) + sign(L2)*pow(fabs(L2), 1./3.);
    if (L < minL || i==1) {
      minL = L;
      mode = i;
    }
    if (*it > 0) {
      sigma += pow(i-mean,2)*(*it)/area;
      skew  += pow(i-mean,3)*(*it)/area;
    }
    it++;
  }
  sigma = sqrt(sigma);
  
  // update to remove tail effect 
  bool redo = false;
  int start = mean - 2*sigma;
  int end   = mean + 2*sigma;
  if (start < 1) start = 1;
  if (end > l-1) end = l-1;
  if (start != 1 || end != l-1) redo = true;
  if (redo) {
    minL = 0;
    for (int i=start; i<end; i++) {
      L1 = (xxcf[i]-xxcf[start-1]) - 2.*i*(xcf[i]-xcf[start-1]) + i*i*(cf[i]-cf[start-1]);
      L2 = (xxcf[end+1]-xxcf[i+1]) - 2.*i*(xcf[end+1]-xcf[i+1]) + i*i*(cf[end+1]-cf[i+1]);
      L  = sign(L1)*pow(fabs(L1), 1./3.) + sign(L2)*pow(fabs(L2), 1./3.);
      if (L < minL || i==start) {
        minL = L;
        mode = i;
      }
    }
  }
  
  // Derive other quantities
  lsigma = xxcf[mode] - 2.*mode*xcf[mode] + mode*mode*cf[mode];
  lsigma = sqrt(fabs(minL/area)*pow(lsigma, 2./3.));
  rsigma = xxcf[l]-2.*mode*xcf[l]+mode*mode*cf[l];
  rsigma -= xxcf[mode+1]-2.*mode*xcf[mode+1]+mode*mode*cf[mode+1];
  rsigma = sqrt(fabs(minL/area)*pow(rsigma, 2./3.));
  amplitude = area*c1/(rsigma+lsigma);
  //Derive();
  
  // Adjust for offset
  mode += offset;
  mean += offset;
}
//______________________________________________________________________________
SplitGaussian SplitGaussian::Sub(double b1, double b2) const {
  if (b2<b1) {
    double tmp = b1;
    b1 = b2;
    b2 = tmp;
  }
  // useful varibles
  double c = sqrt(0.5);
  double kl = lsigma/sqrt(2*M_PI);
  double kr = rsigma/sqrt(2*M_PI);
  double g1l = exp(-0.5*(b1-mode)*(b1-mode)/lsigma/lsigma);
  double g1r = exp(-0.5*(b1-mode)*(b1-mode)/rsigma/rsigma);
  double g2l = exp(-0.5*(b2-mode)*(b2-mode)/lsigma/lsigma);
  double g2r = exp(-0.5*(b2-mode)*(b2-mode)/rsigma/rsigma);
  SplitGaussian g;
  // Calculate the area 
  g.area = Area(b1, b2)+1.e-6;
  // Calculate the mean
  if (b2<mode) g.mean = 0.5*mode*erf(c*(b2-mode)/lsigma)-kl*g2l;
  else         g.mean = 0.5*mode*erf(c*(b2-mode)/rsigma)-kr*g2r;
  if (b1<mode) g.mean -= 0.5*mode*erf(c*(b1-mode)/lsigma)-kl*g1l;
  else         g.mean -= 0.5*mode*erf(c*(b1-mode)/rsigma)-kr*g1r;
  g.mean *= area/g.area;
  // Calculate the sigma 
  if (b2<mode) g.sigma = 0.5*(mode*mode+lsigma*lsigma)*erf(c*(b2-mode)/lsigma) 
                         - kl*(mode+b2)*g2l;
  else         g.sigma = 0.5*(mode*mode+rsigma*rsigma)*erf(c*(b2-mode)/rsigma) 
                         - kr*(mode+b2)*g2r;
  if (b1<mode) g.sigma -= 0.5*(mode*mode+lsigma*lsigma)*erf(c*(b1-mode)/lsigma) 
                         - kl*(mode+b1)*g1l;
  else         g.sigma -= 0.5*(mode*mode+rsigma*rsigma)*erf(c*(b1-mode)/rsigma) 
                         - kr*(mode+b1)*g1r;
  g.sigma *= area;
  // Calculate the skewness 
  if (b2<mode) g.skew = 0.5*mode*(mode*mode+3*lsigma*lsigma)*erf(c*(b2-mode)/lsigma) 
                        - kl*(mode*mode+2*lsigma*lsigma+b2*b2+mode*b2)*g2l;
  else         g.skew = 0.5*mode*(mode*mode+3*rsigma*rsigma)*erf(c*(b2-mode)/rsigma) 
                        - kr*(mode*mode+2*rsigma*rsigma+b2*b2+mode*b2)*g2r;
  
  if (b1<mode) g.skew -= 0.5*mode*(mode*mode+3*lsigma*lsigma)*erf(c*(b1-mode)/lsigma) 
                        - kl*(mode*mode+2*lsigma*lsigma+b1*b1+mode*b1)*g1l;
  else         g.skew -= 0.5*mode*(mode*mode+3*rsigma*rsigma)*erf(c*(b1-mode)/rsigma) 
                        - kr*(mode*mode+2*rsigma*rsigma+b1*b1+mode*b1)*g1r;
  g.skew *= area;
  // make the central moments 
  g.skew = (g.skew - 3*g.mean*g.sigma + 2*g.area*pow(g.mean, 3))/g.area;
  g.sigma = sqrt(fabs(g.sigma - g.area*g.mean*g.mean)/g.area);
  //cout << g.area << " " << g.mean << " " << g.sigma << " " << g.skew << endl;
  g.InverseDerive();
  return g;
}
//______________________________________________________________________________
vector<float> SplitGaussian::Waveform(double b1, double b2, double spacing) const {
  if (spacing==0) spacing=1;
  if (spacing<0)  spacing*=-1;
  if (b2<b1) {
    double tmp = b1;
    b1 = b2;
    b2 = tmp;
  }
  // useful varibles
  double z;
  vector<float> wave;
  
  for (double x=b1; x<=b2; x+=spacing) {
    if (x<=mode) z = (x-mode)/lsigma;
    else         z = (x-mode)/rsigma;
    wave.push_back(amplitude * exp(-0.5*z*z));
  }
  return wave;
}


