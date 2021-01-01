#include <vector>
#include <cmath>
#include "BaselineSanity.h"

using namespace std;

double BaselineSanity(FVec *wave, float max_shift, unsigned int min_width) {
  
  // tracking variables
  vector<int> starts, ends;
  DVec means, cdf, sqcdf, dydx;
  double shift = 0;
  
  // compare variables
  double true_std, cal_mean, cal_std, cal_slope, std_err;
  int nsamples, current, min_step = int(min_width)-1;
  bool pass;
  
  // correcting variables
  double yscale, xscale, mid, tmp1, tmp2, tmp3, value, prev_value;
  
  /////////////////////////////
  // 1st find the true_std 
  /////////////////////////////
  cal_mean = cal_std = 0;
  vector<float>::iterator it = wave->begin();
  int counter = 0;
  for (++it; it != wave->begin()+22; it++) {
  	counter++;
  	cal_mean += *it;
    cal_std += (*it)*(*it);
  }
  cal_mean /= counter;
  cal_std = sqrt(fabs(cal_std-counter*cal_mean*cal_mean)/(counter-1.));
  true_std = cal_std;
  
  cal_mean = cal_std = 0;
  vector<float>::reverse_iterator rit = wave->rbegin();
  counter = 0;
  for (++rit; rit != wave->rbegin()+30; rit++) {
  	counter++;
  	cal_mean += *rit;
    cal_std += (*rit)*(*rit);
  }
  cal_mean /= counter;
  cal_std = sqrt(fabs(cal_std-counter*cal_mean*cal_mean)/(counter-1.));
  if (cal_std < true_std) true_std = cal_std;
  
  /////////////////////////////
  // 2nd, make the cdf and sqcdf. Allows for fast mean and std calculations
  /////////////////////////////
  tmp1 = tmp2 = tmp3 = value = prev_value = 0;
  for (it = wave->begin(); it != wave->end(); it++) {
    value = *it;
    tmp1 += value;
    tmp2 += value * value;
    tmp3 += value - prev_value;
		prev_value = value;
		
    cdf.push_back(tmp1);
    sqcdf.push_back(tmp2);
    dydx.push_back(tmp3);    
  }
  
  current = means.size();
  nsamples = min_step;
  
  
  /////////////////////////////
  // 3rd find flat regions within the waveform 
  /////////////////////////////
  std_err = pow(2.*pow(true_std, 4.)/(min_step-1.0), 0.25);
  for (size_t p=0/*ends[0]*/; p+min_width<wave->size(); p++) {
    // calulate the new mean, std and err 
    cal_mean = (cdf[p+min_step]-cdf[p])/min_step;
    cal_std = sqcdf[p+min_step]-sqcdf[p];
    cal_std = sqrt(fabs(cal_std-min_step*cal_mean*cal_mean)/(min_step-1.0));
    cal_slope = (dydx[p+min_step]-dydx[p])/2.0;
    
    // check if found 
    pass = cal_std < true_std+std_err and fabs(cal_mean) < max_shift 
           and min_width < wave->size() 
           and fabs(cal_slope) < true_std;
    if (!pass) continue;
    
    // record the new mean find and location
    current = means.size();
    nsamples = min_step;
    starts.push_back(p);
    ends.push_back(p+nsamples);
    means.push_back(cal_mean);
    
    // increase the size and update the end 
    do {
      // update values 
      means[current] = cal_mean;
      ends[current] = p+nsamples;
      
      // calulate the new mean and std
      nsamples++;
      cal_mean = (cdf[p+nsamples]-cdf[p])/nsamples;
      cal_std = sqcdf[p+nsamples]-sqcdf[p];
      cal_std = sqrt(fabs(cal_std-nsamples*cal_mean*cal_mean)/(nsamples-1.0));
      cal_slope = (dydx[p+nsamples]-dydx[p])/2.0;
      std_err = pow(2.*pow(true_std, 4.)/(nsamples-1.0), 0.25);
      
      pass = cal_std < true_std+std_err and fabs(cal_mean) < max_shift  
             and p+nsamples+1<wave->size() 
             and fabs(cal_slope) < true_std;
    }while(pass);
    p = ends[current];
    std_err = pow(2.*pow(true_std, 4.)/(min_step-1.0), 0.25);
  }
  
  /////////////////////////////
  // 4th guarantee last points to be part the means -> DO NOT guarantee 
  /////////////////////////////
  current = means.size()-1;
  if (current < 0) return 0;
  
  if (ends[current] > int(wave->size()-31)) ends[current] = wave->size()-1;
  /*
  else {
    size_t p = wave->size()-min_width-1;
    cal_mean = (cdf[p+min_step]-cdf[p])/min_step;
    starts.push_back(p);
    ends.push_back(p+min_step);
    means.push_back(cal_mean);
  }*/
  
  
  // 5.a step: adjust all the points before the first found mean 
  cal_mean = means[0];
  for (it = wave->begin(); it != wave->begin()+starts[0]; it++) {
		*it = *it - cal_mean;
  }
  
  /////////////////////////////
  // 5th subtract the means from the waveform with smooth transitions
  /////////////////////////////
  for (size_t m=0; m+1<means.size(); m++) {
    yscale = (means[m+1]-means[m])/2.0;
    xscale = (starts[m+1]-ends[m])/5.0+0.05;
    mid = (starts[m+1]-ends[m])/2.0+ends[m];
    
    if (shift < yscale*2.0) shift = yscale*2.0;
    counter = starts[m];
    for (it = wave->begin()+starts[m]; it != wave->begin()+starts[m+1]; it++) {
      *it -= yscale*(1.0+erf((counter-mid)/xscale)) + means[m];
      counter++;
    }
  }
  
  /////////////////////////////
  // 6th subtract the last mean 
  /////////////////////////////
  current = means.size()-1;
  cal_mean = means[current];
  for (it = wave->begin()+ends[current]; it != wave->end(); it++) {
    *it = *it - cal_mean;
  }
  
  return shift;
}
