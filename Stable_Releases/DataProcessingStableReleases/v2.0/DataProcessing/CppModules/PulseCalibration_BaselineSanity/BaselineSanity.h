#ifndef __BASELINE_SANITY_H__ 
#define __BASELINE_SANITY_H__ 1
#include <vector>

typedef std::vector<float> FVec;
typedef std::vector<double> DVec;

double BaselineSanity(FVec *wave, float max_shift, unsigned int min_width);
#endif

