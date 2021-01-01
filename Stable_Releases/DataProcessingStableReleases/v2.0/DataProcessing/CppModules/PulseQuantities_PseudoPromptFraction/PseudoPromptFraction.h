#ifndef __PSEUDO_PROMPT_FRACTION_H__
#define __PSEUDO_PROMPT_FRACTION_H__ 1

#include <vector>
#include <iostream>
#include <cmath>

using std::vector;
using std::isfinite;
using std::cerr;
typedef std::vector<float> FVec;

//______________________________________________________________________________
double PseudoPromptFraction(FVec *wave, int start, int end, int center, 
                            int n=9, double threshold=0.1) {
  
  // error checking  
  if (!isfinite(start) or !isfinite(end) or !isfinite(center)) {
    cerr << "***Warning in PseudoPromptFraction: timing information is incorrect\n";
    cerr << "passing inf or nan\n";
    return 0;
  }
  if (end >= int(wave->size())) end = int(wave->size())-1;
  if (start < 0) start = 0;
  if (end - start <= 0) {
    cerr << "***Warning in PseudoPromptFraction: timing information is incorrect\n";
    cerr << "hft_t2_samples <= hft_t0_samples\n";
    return 0;
  }
  if (start > center or center > end) {
    cerr << "***Warning in PseudoPromptFraction: timing information is incorrect\n";
    cerr << "hft_t10l_samples < hft_t0_samples or hft_t10l_samples > hft_t2_samples\n";
    return 0;
  } 
  
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
#endif
