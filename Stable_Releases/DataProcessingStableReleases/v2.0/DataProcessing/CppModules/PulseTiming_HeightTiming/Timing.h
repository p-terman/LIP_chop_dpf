#ifndef __TIMING_H__
#define __TIMING_H__ 1
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
int FirstLeftIndex(FVec *input, int start, int end, double above, double noise) {
  // This function returns the first index where next two values are above a 
  // given threshold or the next value is above a value plus the 5 times the 
  // baseline noise. Forward Iteration.
  int left = start;
  for (; left<end-1; left++) {
    if (input->at(left+1) > above) {
      if (input->at(left+2) > above) break;
    }else if (input->at(left+1) > above+3*noise) break;
  }
  return left;
}
//______________________________________________________________________________
int FirstRightIndex(FVec *input, int start, int end, double above, double noise) {
  // This function returns the first index where next two values are above a 
  // given threshold or the next value is above a value plus the 5 times the 
  // baseline noise. Reverse Iteration.
  int right = start;
  for (; right>end+1; right--) {
    if (input->at(right-1) > above) {
      if (input->at(right-2) > above) break;
    }else if (input->at(right-1) > above+3*noise) break;
  }
  return right;
}
//______________________________________________________________________________
void HeightTiming(FVec *Y, double threshold, double *timing) {
  // This function takes the pulse waveform and returns the start, first 10% 
  // first 50%, max, last 50%, last 10%, and end time of the pulse. 
  // INPUTS:
  //   Y -> pulse waveform in phe/sample 
  //   start -> start time of the pulse relative to the trigger
  //   threshold -> square sum of baseline variance of all contibuting PODs 
  // OUTPUTS:
  //   timing -> t0_idx, t10l_idx .... t10r_idx, t2_idx
  // ASSUMPTIONS:
  //   1) timing variable has double[8] allocated to it 
  //   2) threshold > 0
  
  // Declare the timing indicies
  int t0_idx, t10l_idx, t50l_idx, t1_idx, t50r_idx, t10r_idx, t2_idx;
  int last = int(Y->size())-1;
  
  // Find the maximum value and it's index
  double max = Y->at(0), raw_area = 0;
  t1_idx = 0;
  for (size_t i=1; i<Y->size(); i++) {
    raw_area += (*Y)[i];
    if (max < (*Y)[i]) {
      max = (*Y)[i];
      t1_idx = i;
    }
  }
  timing[3] = t1_idx;
  timing[7] = max;

  
  threshold = sqrt(threshold*threshold+0.05*0.05);
  if (threshold>0.15) threshold = 0.15;
  
  // skip noise pulses 
  if (max <= 3.*threshold) {
    if (t1_idx > 0) {
      timing[0] = t1_idx-1;
      timing[1] = t1_idx-1;
      timing[2] = t1_idx-1;
    }
    else {
      timing[0] = t1_idx;
      timing[1] = t1_idx;
      timing[2] = t1_idx;
    }
    if (t1_idx < last) {
      timing[4] = t1_idx+1;
      timing[5] = t1_idx+1;
      timing[6] = t1_idx+1;
    }
    else {
      timing[4] = t1_idx;
      timing[5] = t1_idx;
      timing[6] = t1_idx;
    }
    return;
  }
  
  // Find the start of the pulse 
  t0_idx = FirstLeftIndex(Y, 0, t1_idx, 1.28*threshold, (3-1.28)*threshold);
  if (t0_idx < 0) t0_idx = 0;
  timing[0] = t0_idx;
  
  // Find the 10% rise time of the pulse 
  t10l_idx = FirstLeftIndex(Y, t0_idx, t1_idx, 0.1*max, threshold);
  timing[1] = t10l_idx;// Interpolation(t10l_idx, Y->at(t10l_idx), t10l_idx+1.0, 
                       //     Y->at(t10l_idx+1), 0.1*max);
  
  // Find the 50% rise time of the pulse 
  t50l_idx = FirstLeftIndex(Y, t10l_idx, t1_idx, 0.5*max, threshold);
  timing[2] = t50l_idx;// Interpolation(t50l_idx, Y->at(t50l_idx), t50l_idx+1.0, 
                       //     Y->at(t50l_idx+1), 0.5*max);
  
  // Find the end of the pulse 
  t2_idx = FirstRightIndex(Y, last, t1_idx, 1.28*threshold, (3-1.28)*threshold);
  if (t2_idx > last) t2_idx = last;
  timing[6] = t2_idx;
  
  // Find the 10% fall time of the pulse 
  t10r_idx = FirstRightIndex(Y, t2_idx, t1_idx, 0.1*max, threshold);
  timing[5] = t10r_idx;// Interpolation(t10r_idx, Y->at(t10r_idx), t10r_idx-1.0, 
                       //     Y->at(t10r_idx-1), 0.1*max);
  
  // Find the 50% fall time of the pulse 
  t50r_idx = FirstRightIndex(Y, t10r_idx, t1_idx, 0.5*max, threshold);
  timing[4] = t50r_idx; // Interpolation(t50r_idx, Y->at(t50r_idx), t50r_idx-1.0, 
                        //    Y->at(t50r_idx-1), 0.5*max);
}
//______________________________________________________________________________
void AreaTiming(FVec *Y, double *timing) {
  
  // Declare the timing indicies
  int t01_idx, t05_idx, t10_idx, t25_idx, t50_idx, t75_idx, t90_idx, t95_idx, t99_idx;
  
  // build cdf 
  FVec cdf;
  double area = 0, threshold = 0.15, max=(*Y)[0];
  int max_idx=0, last = int(Y->size())-1;
  bool above;
  for (int i=0; i<=last; i++) {
    if (max < (*Y)[i]) { max = (*Y)[i]; max_idx = i; }
    above = (*Y)[i] > threshold;
    if (i>0) above = above || (*Y)[i-1] > threshold;
    if (i<last) above = above || (*Y)[i+1] > threshold;
    above = above && (*Y)[i] >= 0;
    
    if (above) area += (*Y)[i];
    cdf.push_back(area);
  }
  
  if (area <= 0) {
  	timing[4] = max_idx;
  	if (max_idx > 0) {
      timing[0] = max_idx-1;
      timing[1] = max_idx-1;
      timing[2] = max_idx-1;
      timing[3] = max_idx-1;
    }
    else {
      timing[0] = max_idx;
      timing[1] = max_idx;
      timing[2] = max_idx;
      timing[3] = max_idx;
    }
    if (max_idx < last) {
      timing[5] = max_idx+1;
      timing[6] = max_idx+1;
      timing[7] = max_idx+1;
      timing[8] = max_idx+1;
    }
    else {
      timing[5] = max_idx;
      timing[6] = max_idx;
      timing[7] = max_idx;
      timing[8] = max_idx;
    }
    timing[9] = cdf[int(timing[8])]-cdf[int(timing[0])];
    return;
  }
  
  // get %area timing 
  t01_idx = t05_idx =t10_idx = t25_idx = t50_idx = 0;
  t75_idx = t90_idx = t95_idx = t99_idx = -1;
  for (size_t i=0; i<cdf.size(); i++) {
    if (cdf[i]/area <= 0.01) t01_idx = i;
    if (cdf[i]/area <= 0.05) t05_idx = i;
    if (cdf[i]/area <= 0.10) t10_idx = i;
    if (cdf[i]/area <= 0.25) t25_idx = i;
    if (cdf[i]/area <= 0.50) t50_idx = i;
    if (cdf[i]/area >= 0.75 and t75_idx==-1) t75_idx = i;
    if (cdf[i]/area >= 0.90 and t90_idx==-1) t90_idx = i;
    if (cdf[i]/area >= 0.95 and t95_idx==-1) t95_idx = i;
    if (cdf[i]/area >= 0.99 and t99_idx==-1) t99_idx = i;
  }
  if (t50_idx == last) t50_idx = last-1;
  if (t75_idx == -1) t75_idx = last;
  if (t90_idx == -1) t90_idx = last;
  if (t95_idx == -1) t95_idx = last;
  if (t99_idx == -1) t99_idx = last;
  
  timing[0] = t01_idx;
  timing[1] = t05_idx;
  timing[2] = t10_idx;
  timing[3] = t25_idx;
  timing[4] = t50_idx;
  timing[5] = t75_idx;
  timing[6] = t90_idx;
  timing[7] = t95_idx;
  timing[8] = t99_idx;
  timing[9] = cdf[t99_idx]-cdf[t01_idx];
 


    
}
#endif
