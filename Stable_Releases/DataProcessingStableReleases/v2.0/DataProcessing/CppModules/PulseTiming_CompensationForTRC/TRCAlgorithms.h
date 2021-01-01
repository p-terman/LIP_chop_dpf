#ifndef __TRC_ALGORITHMS_H__
#define __TRC_ALGORITHMS_H__ 1
/*

Author: Sergey Uvarov
Version: 1.0 created at 3/14/2013

*/

#include <vector>
#include <cmath>

using std::sqrt;
using std::fabs;
using std::vector;
typedef std::vector<float> FVec;

//______________________________________________________________________________
double Interpolation(double x1, double y1, double x2, double y2, double y) {
  // Linear interpolation of the x value given for a given y value.
  if (y1==y2) return (x1+x2)/2.0;
  if (y1>y and y2>y) return (y1<y2) ? x1 : x2;
  return (x2-x1)/(y2-y1) * (y-y1) + x1;
}
//______________________________________________________________________________
FVec SimpleCdf(FVec *Y) {
	FVec cdf;
  double area = 0;
  int last = int(Y->size())-1;
  for (int i=0; i<=last; i++) {
    area += (*Y)[i];
    cdf.push_back(area);
  }
  return cdf;
}
//______________________________________________________________________________
FVec SmarterCdf(FVec *Y, double threshold) {
	FVec cdf;
  double area = 0;
  int last = int(Y->size())-1;
  bool above;
  for (int i=0; i<=last; i++) {
    above = (*Y)[i] > threshold;
    if (i>0) above = above || (*Y)[i-1] > threshold;
    if (i<last) above = above || (*Y)[i+1] > threshold;
    above = above && (*Y)[i] >= 0;
    
    if (above) area += (*Y)[i];
    cdf.push_back(area);
  }
  return cdf;
}
//______________________________________________________________________________
void AftTxTiming(FVec *Y, double threshold, double ratio, double &tlx, double &trx) {
  
  // clear values 
  tlx = trx = 0;
  int last = int(Y->size())-1;
  if (last < 1) return;
  
  // Declare the timing indicies
  int tlx_idx, trx_idx;
  
  // build cdf 
  FVec cdf = SmarterCdf(Y, threshold);
  double area = cdf[last];
  if (area <= 0) {
  	tlx = trx = -1;
  	return;
  } 
  
  // get %area timing 
  tlx_idx = 0;
  trx_idx = -1;
  for (size_t i=0; i<cdf.size(); i++) {
    if (cdf[i]/area <= ratio) tlx_idx = i;
    if (cdf[i]/area >= 1.0-ratio and trx_idx==-1) trx_idx = i;
  }
  if (tlx_idx == last) tlx_idx = last-1;
  if (trx_idx < tlx_idx) trx_idx = last;
  
  tlx = Interpolation(tlx_idx, cdf[tlx_idx]/area, tlx_idx+1., cdf[tlx_idx+1]/area, ratio);
  trx = Interpolation(trx_idx, cdf[trx_idx]/area, trx_idx-1., cdf[trx_idx-1]/area, 1.-ratio);
}
//______________________________________________________________________________
void SkinnyPulseTiming(FVec *Y, double threshold, int width, double &start, double &end) {
	// clear values 
  start = end = 0;
  int last = int(Y->size())-1;
  if (last < 1) return;
  
  // compensate for small pulses (length wise)
  if (last < width) {
  	start = 0;
  	end = last;
  	return;
  }  
  
  // search for maximum area 
  int forward = width - 1;
  FVec cdf = SmarterCdf(Y, threshold);
  double max = cdf[forward];
  start = 0;
  end = forward;
  for (int i=1; i<=last-forward; i++) {
  	if (max < cdf[i+forward]-cdf[i-1]) {
  		max = cdf[i+forward]-cdf[i-1];
  		start = i;
  		end   = i+forward;
  	}
  }
}

#endif
