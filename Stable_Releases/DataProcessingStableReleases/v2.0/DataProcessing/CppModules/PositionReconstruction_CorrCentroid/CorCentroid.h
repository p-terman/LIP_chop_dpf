#include "../Utilities/PMTInfo.h"
#include <vector>
#include <cmath>


using std::vector;
using std::sqrt;

//take in a vector of doubles and return two doubles
void CorrectedCentroidPosRecon(vector<double> area, double &corr_x ,double &corr_y) {
  // This is the algorithm:
  //        \Sum ( x_i * A_i * F_i * (1+(R_i/24.765)) )
  //  x = ------------------------------------------------
  //           \Sum ( A_i * F_i * (1+(R_i/24.765)) )
  // 
  
  corr_x = 0;
  corr_y = 0;
   
  // There is no S2 pulse, don't try to do pos recon.
  double x_numerator = 0;
  double x_denominator = 0;
  double y_numerator = 0;
  double y_denominator = 0;
  double total_peak_light = 0;
  double dist_to_ptfe = 24.765; // cm
  bool top;
  
  // Loop over all channels to calculate total light
  for (size_t i=0; i < area.size(); i++) {    
    top = i+1 == 121 || i+1 < 61;
    if (!top) continue;     // Continue if not a top pmt.
    
    total_peak_light += area[i];
  }
  
  // Loop over all channels.
  double x, y, r, f;
  for (size_t i=0; i<area.size() && total_peak_light>0; i++) {    
    // Grab the peak for current channel.
    top = i+1 == 121 || i+1 < 61;
    if (!top) continue;     // Continue if not a top pmt.
    
    x = PMT[i][0]; // PMT position
    y = PMT[i][1]; // PMT position
    f = area[i] / total_peak_light;
    r = sqrt(x*x + y*y);
    x_numerator   += x * area[i] * f * (1+ r/dist_to_ptfe);
    x_denominator +=     area[i]  * f * (1+ r/dist_to_ptfe);
    y_numerator   += y * area[i] * f * (1+ r/dist_to_ptfe);
    y_denominator +=     area[i] * f * (1+ r/dist_to_ptfe);
  }
  
  if (total_peak_light > 0){
		if (x_denominator > 0) corr_x = x_numerator / x_denominator;
  	if (y_denominator > 0) corr_y = y_numerator / y_denominator;
	}
  // Add a scaling factor to push the events to the edge of the detector:
  corr_x *= 1.22;
  corr_y *= 1.22;
}
