#ifndef __PULSE_FINDER_H__
#define __PULSE_FINDER_H__ 1

#include <vector>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <list>
#include "SplitGaussian.h"

typedef std::vector<float> FVec;
typedef std::vector<double> DVec;
typedef std::vector<int> IVec;

typedef struct {
	double area;
	int index;
}AreaOrder;

double Overlap(SplitGaussian &lhs, SplitGaussian &rhs, double b1, double b2);

class PulseFinder {
	protected:
		// Main operation functions
		double BDistance(int left, int right);
		
		// Supporting functions
		void MakeOrderList();
		void BuildPars(int start, int end, SplitGaussian &g);
		void RemoveAndUpdate(int index, bool left, bool right);
		
		double se_gap;
		double se_amp;
		double se_err;
		double smooth_max_threshold;
		std::list<AreaOrder> order_list;
		std::vector<bool> changed;
		
	public:
		// Parameters
		int kz_samples;
		int kz_iterations;
		double min_max_ratio;
		double max_threshold;
		int max_seperation;
		double nsigma;
		double se_area_mean_phe;
		double se_area_sigma_phe;
		double se_width_mean_samples;
		double se_width_sigma_samples;
		double dis_inv_ratio;
		double dis_scale;
		double dis_offset;
		
		// Containers
		FVec wave;
		FVec smooth;
		IVec maxs;
		IVec mins;
		IVec borders;
		std::vector<SplitGaussian> pars;
		
		void KZFilter();
		void GetMaximums();
		void GetMinimums();
		void BuildSplitGaussians();
		void SEFinder(bool includeGrouping);
		void IteritiveCluster();
		void BaselineTrim();
		
		PulseFinder() { 
			kz_samples=10; kz_iterations=2; min_max_ratio=0.6; max_threshold=0.15;
			max_seperation=8; nsigma=3.0; se_area_mean_phe=24.55;
			se_area_sigma_phe=5.799; se_width_mean_samples=44.23;
			se_width_sigma_samples=11.81;dis_inv_ratio=0.01;dis_scale=1.0;
			dis_offset=1.2;
		}
		void Initilize();
		void Execute() {
			Clear();
			KZFilter();
			GetMaximums();
			GetMinimums();
			BuildSplitGaussians();
			SEFinder(true);
			SEFinder(false);
			IteritiveCluster();
			BaselineTrim();
		}
		void Clear() {
			smooth.clear();
			maxs.clear();
			mins.clear();
			borders.clear();
			pars.clear();
			order_list.clear();
			changed.clear();
		}
};

#endif

