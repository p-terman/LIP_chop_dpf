
#include <cmath>
#include <algorithm>
#include "PulseFinder.h"

#define DEBUG 0
#define SMALLNUMBER 0.000001

using namespace std;

////////////////////////////////////////////////////////////////////////////////
// Independent Functions
////////////////////////////////////////////////////////////////////////////////
double PulseFinder::BDistance(int l, int r) {
  double ratio, mean, sigma, asymmetry, alpha, w1, w2, min, max;
  double ampl1, ampl2, s2like, s2like1, s2like2;
  double sq3 = sqrt(3.0);
  double sq2 = sqrt(2.0);
  if (pars[r].area < pars[l].area) {
    min = pars[r].area;
    max = pars[l].area;
    mean = pars[l].mean;
    sigma = pars[l].sigma;
    ampl1 = pars[l].amplitude;
    ampl2 = pars[r].amplitude;
    s2like1 = pars[l].sigma*2.35/se_width_mean_samples;
    s2like2 = pars[r].sigma*2.35/se_width_mean_samples;
  }else {
    min = pars[l].area;
    max = pars[r].area;
    mean = pars[r].mean;
    sigma = pars[r].sigma;
    ampl1 = pars[r].amplitude;
    ampl2 = pars[l].amplitude;
    s2like1 = pars[r].sigma*2.35/se_width_mean_samples;
    s2like2 = pars[l].sigma*2.35/se_width_mean_samples;
  }
  ratio = min/max;
  alpha = pow(2*ratio/(ratio*ratio+1.), 0.25);
  s2like = pow(2*s2like1*s2like2/(s2like1+s2like2),4);
  w1 = pow((ampl1-ampl2)/(ampl1+ampl2),2);
  w2 = pow(sigma/sqrt(pow(sigma,2)+pow(se_width_mean_samples/2.35,2)),2);
  asymmetry = (pars[l].lsigma/(pars[l].rsigma)-1.)/pars[l].amplitude
            + (pars[r].rsigma/(pars[r].lsigma)-1.)/pars[r].amplitude;
  asymmetry *= ampl1<ampl2 ? ampl1 : ampl2;
  if (asymmetry<0) asymmetry = 0;
  SplitGaussian g = pars[l] + pars[r];
  
  // calculate distances 
  double d1 = pars[l].BhattacharyyaDistance(pars[r])*alpha;// + dis_offset;
  double d2 = pars[l].BhattacharyyaDistance2(pars[r])*alpha;
  //double d3 = (pars[r].mean-sq3*pars[r].sigma-pars[l].mean-sq3*pars[l].sigma)/(sq2*se_gap); //d2;
  double d3 = ((pars[r].mode-sq3*pars[r].lsigma)-(pars[l].mode+sq3*pars[l].rsigma))/(sq2*se_gap);
  double d4 = fabs(g.mean-mean)/(sigma/sqrt(max/se_area_mean_phe));
  
  // modify distances 
  //d1 -= asymmetry*s2like1*s2like1 + dis_scale/(s2like2+1./s2like1);
  //d1 -= (asymmetry + dis_scale)/(s2like2+1./s2like1);
  //d1 -= (asymmetry+dis_scale)*s2like1*s2like2/(s2like1+s2like2);
  d1 = (d1/s2like + d2*s2like)/(s2like+1/s2like) + dis_offset;
  d2 += dis_offset;
  d3 -= dis_scale*(2*ampl1/(ampl1+ampl2));
  d4 -= asymmetry*s2like1*se_area_mean_phe/min + dis_scale/(s2like2+1./s2like1);
  if (d1<0) d1 = 0;
  if (d2<0) d2 = 0;
  if (d3<0) d3 = 0;
  if (d4<0) d4 = 0;
  
  // overall distance with weight factors included 
  double distance = d1*(1-w1)*(1-w2) + d2*w1*(1-w2) + d3*(1-w1)*w2 + d4*w1*w2;
  
#if DEBUG>1
cout << "\t\tpar left " << pars[l].amplitude << " " << pars[l].mode << " "
     << pars[l].lsigma << " " << pars[l].rsigma << " "
     << pars[l].area << " " << pars[l].mean << " " << pars[l].sigma << endl;
cout << "\t\tpar right " << pars[r].amplitude << " " << pars[r].mode << " "
     << pars[r].lsigma << " " << pars[r].rsigma << " "
     << pars[r].area << " " << pars[r].mean << " " << pars[r].sigma << endl;
cout << "\t\tDistance " << d1 << " " << d2 << " " << d3 << " " << d4
     << " " << w1 << " " << w2 << "\n"
     << "\t\t\t" << d1*(1-w1)*(1-w2) << " + " << d2*w1*(1-w2) 
     << " + " << d3*(1-w1)*w2 << " + " << d4*w1*w2 << "\n"
     << "\t\t\t (" << asymmetry << "," << s2like1 << "," << s2like2 << "," << s2like <<")\n";
cout.flush();
#endif  
  return distance;
}

////////////////////////////////////////////////////////////////////////////////
// Pulse Finder : Main Operations Functions
////////////////////////////////////////////////////////////////////////////////
void PulseFinder::Initilize() {
	
	// FIND THE SMOOTH MAX THRESHOLD
	vector<float> tmp_smooth(kz_samples*kz_iterations, 0);
	int index = (kz_samples*kz_iterations)/2;
	tmp_smooth[index] = max_threshold;
	
	int lower = floor(-0.5*kz_samples);
  int upper = floor(0.5*(kz_samples+1));
  int length = tmp_smooth.size();
  double weight, total, ty;
	vector<float> tmp;
	
	for (int i=0; i<kz_iterations; i++) {    
    tmp = tmp_smooth;
    for (int j=0; j<length; j++) { 
      total = 0;
      ty = 0;
      for (int k=lower; k<=upper; k++) {
        weight = 0.5*(1.0-cos((2.0*M_PI*(k-lower))/(kz_samples)));
        total += weight;
        if (j+k >= 0 && j+k < length) {
          ty += tmp[j+k] * weight;
        }
      }
      tmp_smooth[j] = ty/total;
    }
  }	
	smooth_max_threshold = tmp_smooth[index];
	
	// other 
	se_gap = se_width_mean_samples/2.35*4.29/3.464;
	se_amp = 2*se_area_mean_phe/(se_width_mean_samples/2.35*4.29); // =0.304
  se_err = se_amp*sqrt(pow(se_area_sigma_phe/se_area_mean_phe,2)+
                       pow(se_width_sigma_samples/se_width_mean_samples,2)); // = 0.108
}
//______________________________________________________________________________
void PulseFinder::KZFilter() {
	// KZ filter is a repeatition of a normalized box filter 
  int lower = floor(-0.5*kz_samples);
  int upper = floor(0.5*(kz_samples+1));
  int length = wave.size();
  double weight, total, ty, threshold = max_threshold;
  FVec tmp;
  // remove majority of the baseline 
  for (int i=0; i<length; i++) {
  	if (i==0 || i+1==length) smooth.push_back(0);
  	else if (wave[i]>threshold || wave[i+1]>threshold || wave[i-1]>threshold) {
  		if (wave[i]>0) smooth.push_back(wave[i]);
  		else smooth.push_back(0);
  	} 
  	else smooth.push_back(0);
  }
  
  // smooth out to minimize other baseline 
  for (int i=0; i<kz_iterations; i++) {    
    tmp = smooth;
    for (int j=0; j<length; j++) { 
      total = 0;
      ty = 0;
      for (int k=lower; k<=upper; k++) {
        weight = 0.5*(1.0-cos((2.0*M_PI*(k-lower))/(kz_samples)));
        total += weight;
        if (j+k >= 0 && j+k < length) {
          ty += tmp[j+k] * weight;
        }
      }
      smooth[j] = ty/total;
    }
  }
}
//______________________________________________________________________________
void PulseFinder::GetMaximums() {
	
  int width=floor(0.5*max_seperation), max_index=0, length=smooth.size();
  bool max_found, min_found;
  double min_value, min, max_value=0, threshold = smooth_max_threshold;
  
  for (int i=width; i<length-width; i++) {
    max_found = true;
    // 1st condition: above the max threshold
    if (smooth[i] < threshold) continue; 
    
    // 2nd condition: two maximums are seperated by a specified distance
    if (maxs.size()) {
      if (i-max_index < max_seperation) {
        if (smooth[i] > smooth[max_index]) {
          max_index = i;
          maxs[int(maxs.size())-1] = max_index;
        }
        continue;
      } 
    }
    
    // 3rd condition: highest value in the range 
    for (int j=1; j<width; j++) {
      max_found = max_found && smooth[i] >= smooth[i+j] 
                            && smooth[i] >= smooth[i-j];
    	if (!max_found) break;
    }
    if (!max_found) continue;
    
    //cout << "max potential " << i << " " << smooth[i] << endl;
    
    // 4th condition: a valid minimum exist between two maximums 
    if (maxs.size()) {      
      // Find the minimum between the current candidate and the previous max.
      // Using the square of the value to avoid sudden charge dips due to 
      // large pulses. 
      min_value = smooth[i];
      min = smooth[i]*smooth[i];
      for (int j=i-1; j>max_index; j--)
        if (min > smooth[j]*smooth[j]) {
          min = smooth[j]*smooth[j];
          min_value = smooth[j];
        }
      // a valid minimum is when the min lower by a specified percent than the 
      // max and the difference between the max and min is above the noise.
      // This applies to both maximums. 
      min_found = min_max_ratio > fabs(min_value/smooth[i])
                  && min_max_ratio > fabs(min_value/max_value)
                  && min_value < smooth[i] - threshold
                  && min_value < max_value - threshold;
            
      if (!min_found) {
        // If no minimum is found, then one of the maximums is invalid.
        // So the larger max is kept. 
        if (max_value > smooth[i]) max_found = false; 
        else {
          max_found = false;
          max_index = i;
          max_value = smooth[max_index];
          maxs[int(maxs.size())-1] = max_index;
        }
      }
    } // end of looking for minimums
    
    if (max_found) {
      max_index = i;
      max_value = smooth[max_index];
      maxs.push_back(max_index);
      i += max_seperation-1;
    }
  }
#if DEBUG>0
	cout << "Maximum:: ";
	for (size_t i=0; i<maxs.size(); i++) 
		cout << maxs[i] << " ";
	cout << endl;
	cout.flush();
#endif
}
//______________________________________________________________________________
void PulseFinder::GetMinimums() {
	// The minimums include the first and last points 
	mins.push_back(0);
  double min=1e6, value;
  int idx_big, idx_small, min_idx, nmaxs=maxs.size(), step;
  int width = floor(0.5*max_seperation);
  
  for (int i=1; i<nmaxs; i++) {
    // want to search from smaller maximum towards the bigger maximum 
    if (smooth[maxs[i]] > smooth[maxs[i-1]]) {
    	idx_big = maxs[i]-width;
    	idx_small = maxs[i-1]+width;
    	step = 1;
    }else {
    	idx_big = maxs[i-1]+width;
    	idx_small = maxs[i]-width;
    	step = -1;
    }
    // find the average minimum of over the sample of points. 
    min_idx = idx_small;
    for (int j=idx_small; j!=idx_big; j+=step) {
    	value = smooth[j]+wave[j]; 
    	for (int k=1; k<width; k++) {
    		value += smooth[j+k]+wave[j+k] + smooth[j-k]+wave[j-k];
    	}
    	
    	if (min>value || j==idx_small) {
    		min = value;
    		min_idx = j;
    	}    	
    }
    mins.push_back(min_idx);
  }
  mins.push_back(int(smooth.size()));
  
  // assume the minimums give correct borders 
	borders = mins;
#if DEBUG>0
	cout << "Minimum:: ";
	for (size_t i=0; i<mins.size(); i++) 
		cout << mins[i] << " ";
	cout << endl;
	cout.flush();
#endif
}
//______________________________________________________________________________
void PulseFinder::BuildSplitGaussians() {
  //if (borders.size() < 2) return;
  pars.resize(borders.size()-1);
  changed.resize(borders.size()-1);

  for (size_t i=0; i<pars.size(); i++) {
    pars[i].Fit(smooth.begin()+borders[i], smooth.begin()+borders[i+1], borders[i]);
    //pars[i].Derive();
    changed[i] = true;
#if DEBUG>0
	  cout << "Build SP: " 
	       << i << " [" << borders[i] << "," << borders[i+1] << "]\t"
			   << "par(a,m,ls,rs,mm,s,sk)=(" << pars[i].area << "," << pars[i].mode 
			   << "," << pars[i].lsigma << "," << pars[i].rsigma << ","
			   << pars[i].mean << "," << pars[i].sigma << ","
			   << pars[i].skew << ")\n";
		cout.flush();
#endif
  }
}
//______________________________________________________________________________
void PulseFinder::SEFinder(bool includeGrouping) {
  
  SplitGaussian g, prev_g;
  size_t startIdx, endIdx;
  double width, distance, d1, d2, d3, d4, sep, gap;
  double prev_distance, cur_distance, mod_nsigma;
  
  if (includeGrouping) mod_nsigma = nsigma;
  else mod_nsigma = nsigma+0.25;
  
  // amplitude parameters
  double upper = se_amp + mod_nsigma*se_err;
  
  // loop through the SplitGaussian
  for (size_t i=0; i+1<pars.size(); i++) {
    // Set initial conditions
    startIdx = endIdx = i;
    prev_g = g = pars[startIdx];
    d3 = pars[startIdx].amplitude;
    d4 = 0;
    prev_distance = cur_distance = 300;
    
    // loop through possibilities 
    for (size_t j=startIdx+1; j<pars.size(); j++) {
      g = g + pars[j];
      width = g.sigma*2.35;
      gap = ((pars[j].mean-pars[j].sigma) - (pars[j-1].mean+pars[j-1].sigma));
      
      d1 = (g.area-se_area_mean_phe)/se_area_sigma_phe;
      d2 = (width-se_width_mean_samples)/se_width_sigma_samples;
      d3 = d3<pars[j].amplitude ? pars[j].amplitude : d3;
      d4 = gap>d4 ? gap : d4;
      distance = sqrt(pow(d1, 2) + pow(d2, 2));
      
      sep = 1;
      if (j<pars.size() && includeGrouping) {
        sep = (pars[j+1].mean-pars[j].mean)/pars[j+1].sigma - (pars[j].mean-prev_g.mean)/prev_g.sigma;
      }
      
      if (d3<upper && distance<mod_nsigma && sep >= 1 && d4<mod_nsigma*se_gap) {
        prev_distance = cur_distance;
        cur_distance = distance;
        endIdx = j;
      }else if (endIdx != startIdx) { // check if the previous SE was found 
        prev_distance = cur_distance;
        cur_distance = 300;
        endIdx = j;
      }
#if DEBUG>2
	  cout << "  SE Loop: "
	       << i << " " << j << " "
			   << "par(a,mm,s)=(" << g.area << "," << g.mean << "," << g.sigma << ") "
			   << "d(1,2,3)=(" << d1 << "," << d2 << "," << d3 << "," << d4 << ")->" 
			   << distance << "\n" 
			   << "\t\t\tpasses " << bool(d3<upper) << " " << bool(sep>=1) << " " 
			   << bool(d4<mod_nsigma*se_gap) << "\n";
		cout.flush();
#endif 
      
      if (prev_distance < mod_nsigma) {
        // found an SE 
        if (prev_distance < cur_distance) break;
      } 
      else if (d1 > mod_nsigma || d2 > mod_nsigma) break; // cannot possibly find any more SE
      prev_g = g;
    }
    
    // select the minimum distance and combine 
    if (prev_distance<mod_nsigma || cur_distance<mod_nsigma) {
      if (prev_distance<cur_distance) {
        endIdx--;
      }
      if (startIdx == endIdx) endIdx++;
      
      // combine 
      borders.erase(borders.begin()+startIdx+1, borders.begin()+endIdx+1);
      pars.erase(pars.begin()+startIdx+1, pars.begin()+endIdx+1);
      changed.erase(changed.begin()+startIdx+1, changed.begin()+endIdx+1);
      BuildPars(borders[startIdx], borders[startIdx+1], pars[startIdx]);
      
#if DEBUG>1
	  cout << "SEFinder: "
	       << "[" << borders[startIdx] << "," << borders[startIdx+1] << "] "
	       << startIdx << " " << endIdx << "\t"
	       << "distance(p,c,)=(" << prev_distance << ","
	       << cur_distance << ")\n";
	  cout.flush();
#endif      
      i--;
    }
  }
  
#if DEBUG>0
	cout << "SEFinder Results\n";
	for (size_t i=0; i+1<borders.size(); i++) 
	cout << i << " [" << borders[i] << "," << borders[i+1] << "]\t"
			   << "par(a,m,ls,rs)=(" << pars[i].area << "," << pars[i].mode 
			   << "," << pars[i].lsigma << "," << pars[i].rsigma << ")\n";
  cout.flush();
#endif
}
//______________________________________________________________________________
void PulseFinder::IteritiveCluster() {
	
	// helpful variables 
	size_t idx;
	double llDis, rrDis, lD, rD;
	double lSeperated, rSeperated;
	bool lCombine, rCombine, lExist, llExist, rExist, rrExist;
	
	// Make list of smallest to largest pulse candidates 
	MakeOrderList();
	list<AreaOrder>::iterator it = order_list.begin();
	list<AreaOrder>::iterator it_last = order_list.end();
	
	// Loop through regions 
	for (; it != it_last && pars.size()>1;) {
		
		// ignore pulse candiates which had nothing happen around them 
		if (!changed[it->index]) continue;
		
		// setup  
		idx = it->index;
		lExist = rExist = llExist = rrExist = false;
		if (idx>0) if (pars[idx].area>pars[idx-1].area) lExist = true;
		if (idx>1) llExist = true; // if (pars[idx-1].area<pars[idx-2].area) llExist = true;
		if (idx+1<pars.size()) if (pars[idx].area>pars[idx+1].area) rExist = true;
		if (idx+2<pars.size()) rrExist = true; //if (pars[idx+1].area<pars[idx+2].area) rrExist = true;
		lD = rD = llDis = rrDis = 100*nsigma;
		lSeperated = rSeperated = 100*nsigma;
		lCombine = rCombine = false;
		
#if DEBUG>1
	  cout << "IteritiveCluster: "
	       << idx << " [" << borders[idx] << "," << borders[idx+1] << "]\t"
			   << "par(a,m,ls,rs)=(" << pars[idx].area << "," << pars[idx].mode 
			   << "," << pars[idx].lsigma << "," << pars[idx].rsigma << ")\n";
		cout.flush();
#endif
		
		// get some quick stats 
		if (lExist) {
			lD = BDistance(idx-1, idx);
			lSeperated = pars[idx].mode - pars[idx-1].mode;
			lSeperated /= pars[idx].mean - pars[idx-1].mean;
			
			if (llExist) {
			  llDis = BDistance(idx-2, idx-1);
			}
			
			lCombine = (lD<nsigma && llDis>lD) || lSeperated < 0.35;
#if DEBUG>1
	  cout << "\tleft "
	       << lCombine << " " << lD << " " << llDis << " " << lSeperated << "\n";
	  cout.flush();
#endif
		}		
		
		if (rExist) {
			rD = BDistance(idx, idx+1);
			rSeperated = pars[idx+1].mode - pars[idx].mode;
			rSeperated /= pars[idx+1].mean - pars[idx].mean;
			
			if (rrExist) {
			  rrDis = BDistance(idx+1, idx+2);
			}
			
			rCombine = (rD<nsigma && rrDis>rD) || rSeperated < 0.35;
#if DEBUG>1
	  cout << "\tright "
	       << rCombine << " " << rD << " " << rrDis << " " << rSeperated << "\n";
#endif
		}
    
    // combine pulse candidates
		if (lCombine || rCombine) {
			RemoveAndUpdate(idx, lCombine, rCombine);
			it = order_list.begin();
			it_last = order_list.end();
			continue;
		}
		
		// set that this pulse candidate cannot combine with its neighbors
		changed[idx] = false;
		it++;
	}
	
	
#if DEBUG>0
	cout << "IteritiveCluster Results\n";
	for (size_t i=0; i+1<borders.size(); i++) 
	cout << i << " [" << borders[i] << "," << borders[i+1] << "]\t"
			   << "par(a,m,ls,rs)=(" << pars[i].area << "," << pars[i].mode 
			   << "," << pars[i].lsigma << "," << pars[i].rsigma << ")\n";
	cout.flush();
#endif
}
//______________________________________________________________________________
void PulseFinder::BaselineTrim() {
	
	// find the standard deviation of the baseline
	double threshold = 0.03;
	int border, length = wave.size(), last=int(borders.size()-1);
	
	// the new maximums between the borders 
	vector<int> newmaxs;
	for (size_t j=1; j<borders.size(); j++) {
		double max=0;
		int max_location = -1; 
		for (size_t k=0; k<maxs.size(); k++) {
			if (maxs[k]>borders[j]) break;
			if (maxs[k]>borders[j-1] && max<smooth[maxs[k]]) {
				max = smooth[maxs[k]];
				max_location = maxs[k];
			}
		}
		if (max_location == -1) { // if no maximum is found 
			for (int k=borders[j-1]; k<borders[j]; k++) {
				if (max<wave[k]) {
					max = wave[k];
					max_location = k;
				}
			}
		}
		newmaxs.push_back(max_location);
	}	
	
	threshold = max_threshold;
	IVec newborders;
	
	// remove prebuffer baseline  
	border = borders[0];
	for (int j=borders[0]+1; j<newmaxs[0]; j++) {
		border = j-1;
		if (wave[j] > threshold) break;
	}
	newborders.push_back(border);
	
	// any borders in between the original waveform
	for (int i=1; i<last; i++) {
		border = borders[i];
		for (int j=borders[i]-2; j>newmaxs[i-1]; j--) {
			border = j+2;
			if (wave[j] > threshold) break;
		}
		// append 20 extra samples to border 
		border = (border+30<borders[i]) ? border+30 : borders[i];
		newborders.push_back(border);
		
		border = borders[i];
		for (int j=borders[i]+1; j<newmaxs[i]; j++) {
			border = j-1;
			if (wave[j] > threshold) break;
		}
		newborders.push_back(border);
	}
	// remove postbuffer baseline
	border = borders[last];
	for (int j=borders[last]-2; j>newmaxs[last-1]; j--) {
		border = j+2;
		if (wave[j] > threshold) break;
	}
	border = (border+30<borders[last]) ? border+30 : borders[last];
	newborders.push_back(border);
	
	borders = newborders;
	
	// add any extra space if possible 
	last = borders.size()-1;
	if (borders[0] > 0) borders[0]--;
	for (size_t i=1; i+1<borders.size(); i+=2) {
		if (borders[i+1]-borders[i]>2) {
			borders[i]+=2;
			borders[i+1]--;
		}
	}
	if (borders[last]+2 < length) borders[last]+=2;
	
#if DEBUG>0
	cout << "BaselineTrim: ";
	for (size_t i=0; i<borders.size(); i++) 
		cout << borders[i] << " ";
	cout << endl;
	cout.flush();
#endif
}

////////////////////////////////////////////////////////////////////////////////
// Pulse Finder : Supporting Functions
////////////////////////////////////////////////////////////////////////////////
bool SortCompareFunc(AreaOrder a, AreaOrder b) {
	return a.area > b.area;
}
void PulseFinder::MakeOrderList() {
	order_list.clear();
	AreaOrder a;
	for (unsigned int j=0; j<pars.size(); j++) {
		if (!changed[j]) continue;
		a.area = pars[j].area; a.index = j;
		order_list.push_back(a);
	}
	order_list.sort(SortCompareFunc);
}
//______________________________________________________________________________
void PulseFinder::BuildPars(int start, int end, SplitGaussian &par) {
	par.Fit(smooth.begin()+start, smooth.begin()+end, start);
}
//______________________________________________________________________________
void PulseFinder::RemoveAndUpdate(int index, bool left, bool right) {
	
	
	if (left && right) {
		// erase pulse candidate
		borders.erase(borders.begin()+index+1);
		pars.erase(pars.begin()+index+1);
		changed.erase(changed.begin()+index+1);
		
		borders.erase(borders.begin()+index);
		pars.erase(pars.begin()+index-1);
		changed.erase(changed.begin()+index-1);
    
		index--;
		
	}else if (left) {
    // erase pulse candidate
		borders.erase(borders.begin()+index);
		pars.erase(pars.begin()+index-1);
		changed.erase(changed.begin()+index-1);
    // update changed array 
		index--;
	
	}else if (right) {
    // erase pulse candidate
		borders.erase(borders.begin()+index+1);
		pars.erase(pars.begin()+index+1);
    changed.erase(changed.begin()+index+1);
	}
	
	// update changed array 
	changed[index] = true;
	BuildPars(borders[index], borders[index+1], pars[index]);
	if (index>0) changed[index-1] = true;
	if (index+1<int(changed.size())) changed[index+1] = true;
	MakeOrderList();
}

