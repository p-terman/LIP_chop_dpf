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
void BasicQuantites(FVec *wave, size_t start, size_t end, double *results) {
  // results are
  // [0] area_phe
  // [1] height_phe_per_sample
  // [2] area_positive_phe
  // [3] area_negative_phe
  
  results[0] = results[1] = results[2] = results[3] = 0;
  results[1] = 0;
  for (size_t i=start; i<=end; i++) {
    results[0] += wave->at(i);
    if (wave->at(i) >= 0) results[2] += wave->at(i);
    else results[3] += wave->at(i);
    if (results[1] < wave->at(i)) results[1] = wave->at(i);
  }  
}


#endif
