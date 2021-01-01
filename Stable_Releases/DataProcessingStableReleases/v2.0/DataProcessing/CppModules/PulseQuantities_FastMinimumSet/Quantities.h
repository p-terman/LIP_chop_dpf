#ifndef __QUANTITIES_H__
#define __QUANTITIES_H__ 1

#include <vector>
#include <iostream>
#include <cmath>

using std::vector;
using std::fabs;
using std::floor;
using std::sqrt;

typedef std::vector<float> FVec;
typedef std::vector< vector<float> > D2Vec;

//______________________________________________________________________________
void S2Filter(FVec *wave, size_t S1_width, size_t S2_width, 
              double &s1_max, double &s2_max) {
  // return s1_max and s2_max
    
  vector<double> cdf;
  double sum = 0;
  for (size_t i=0; i<wave->size(); i++) {
    sum += wave->at(i);
    cdf.push_back(sum);
  }
  
  // setup the lower and upper limits
  int s1_lower = floor(S1_width/2.0);
  int s1_upper = floor((S1_width+1)/2.0);
  int s2_lower = floor(S2_width/2.0);
  int s2_upper = floor((S2_width+1)/2.0);  
  int last = int(cdf.size())-1;
  double s1_value1, s1_value2, s2_value1, s2_value2;
  
  for (int i=0; i<=last; i++) {
    // s1 section
    if (i-s1_lower-1 < 0) s1_value1 = 0;
    else s1_value1 = cdf[i-s1_lower-1];
    if (i+s1_upper >= last) s1_value2 = sum;
    else s1_value2 = cdf[i+s1_upper];
    if (i==0) s1_max = s1_value2-s1_value1;
    if (s1_max < s1_value2-s1_value1) s1_max = s1_value2-s1_value1;
    
    // s2 section
    if (i-s2_lower-1 < 0) s2_value1 = 0;
    else s2_value1 = cdf[i-s2_lower-1];
    if (i+s2_upper >= last) s2_value2 = sum;
    else s2_value2 = cdf[i+s2_upper];
    if (i==0) s2_max = s2_value2-s2_value1;
    if (s2_max < s2_value2-s2_value1) s2_max = s2_value2-s2_value1;
  }
}
//______________________________________________________________________________
double PromptFraction(FVec *wave, double pulse_area, size_t t10l, 
                      size_t prebins=5, size_t window=4) {
  double area = 0;
  size_t start = 0;
  if (int(t10l)-int(prebins) > 0) start = t10l-prebins;
  size_t end = t10l+window;
  if (end >= wave->size()) end = wave->size()-1;
  for (size_t i=start; i<=end; i++)
    area += wave->at(i);
  return area/pulse_area; 
}
//______________________________________________________________________________
double PseudoPromptFraction(FVec *wave, int start, int end, int center, 
                            int n=9, double threshold=0.1) {  
  double nsigma = 3.0;  
  double area = 0;
  for (int i=start; i<=end; i++) {
    if (wave->at(i) >= nsigma*threshold) area += wave->at(i);
  }
  if (area <= 0) return 1;
  
  double spf = 0;
  int counter = 0;
  for (int i=center; i>=start; i--) {
    if (wave->at(i) >= nsigma*threshold) {
      spf += wave->at(i);
      counter++;
    }
    if (counter > n) break;
  }
  counter = 0;
  for (int i=center+1; i<=end; i++) {
    if (wave->at(i) >= nsigma*threshold) {
      spf += wave->at(i);
      counter++;
    }
    if (counter > n) break;
  }
  return spf/area;
}
//______________________________________________________________________________
void PulseBasicQuantites(FVec *wave, size_t start, size_t end, double *results) {
  // results are
  // [0] pulse_area_phe
  // [1] pulse_height_phe_per_sample
  // [2] pulse_area_positive_phe
  // [3] pulse_area_negative_phe
  
  results[0] = results[1] = results[2] = results[3] = 0;
  results[1] = 0;
  for (size_t i=start; i<=end; i++) {
    results[0] += wave->at(i);
    if (wave->at(i) >= 0) results[2] += wave->at(i);
    else results[3] += wave->at(i);
    if (results[1] < wave->at(i)) results[1] = wave->at(i);
  }  
}
//______________________________________________________________________________
void PeakQuantities(FVec *wave, size_t start, size_t end, double *results) {
  // results are 
  // [0] peak_area_phe
  // [1] peak_height_phe_per_sample
  // [2] mean_first_last_pts_phe_per_sample
  
  results[0] = 0;
  results[1] = 0;
  for (size_t i=start; i<=end; i++) {
    results[0] += wave->at(i);
    if (results[1] < wave->at(i)) results[1] = wave->at(i);
  }
  if (wave->size() > 0)
    results[2] = (wave->at(0)+wave->at(wave->size()-1))/2.0;
  results[2] = 0;
}
//______________________________________________________________________________
void PODBaselineStats(FVec *wave, double *results) {
  // Results are 
  // [0] pod predetect mean 
  // [1] pod predetect std 
  // [2] pod postdetect mean 
  // [3] pod postdetect std 
  
  int last = int(wave->size())-1;
  results[0] = results[1] = results[2] = results[3] = 0;
  if (last < 31) return;
  
  // calculate the predetect 
  for (int i=0; i<24; i++) {
    results[0] += wave->at(i);
    results[1] += wave->at(i)*wave->at(i);
  }
  results[0] /= 24.0;
  results[1] = sqrt(fabs(results[1]-24.0*results[0]*results[0])/23);
  // calculate the postdetect 
  for (int i=0; i<31; i++) {
    results[2] += wave->at(last-i);
    results[3] += wave->at(last-i)*wave->at(last-i);
  }
  results[2] /= 31.0;
  results[3] = sqrt(fabs(results[3]-31.0*results[2]*results[2])/30);
}


#endif
