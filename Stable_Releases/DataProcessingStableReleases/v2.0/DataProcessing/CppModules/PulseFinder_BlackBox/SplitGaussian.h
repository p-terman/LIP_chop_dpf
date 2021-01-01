#ifndef __SPLIT_GAUSSIAN_H__
#define __SPLIT_GAUSSIAN_H__
#include <cmath>
#include <vector>
#include <iostream>

class SplitGaussian {
  public: 
  	double amplitude;
    double mode;
    double lsigma;
    double rsigma;
    double area;
    double mean;
    double sigma;
    double skew;
    
    SplitGaussian();
    SplitGaussian(double a, double m, double l, double r);
    SplitGaussian operator+(const SplitGaussian &rhs) const;
    SplitGaussian& operator+=(const SplitGaussian &rhs);
    double operator[](double x) const;
    double Area(double b1, double b2) const;
    double CF(double x) const;
    void Derive();
    void InverseDerive();
    double Distance(double x) const;
    double BhattacharyyaDistance(SplitGaussian &g) const;
    double BhattacharyyaDistance2(SplitGaussian &g) const;
    double EarthsMoversDistance(SplitGaussian &g, double b1, double b2) const;
    void Fit(std::vector<float>::iterator s, std::vector<float>::iterator e, double offset);
    void Print(size_t verbosity=0) const {
      std::cout << amplitude << " " << mode << " " << lsigma << " " << rsigma;
      if (verbosity>0) std::cout << " " << area;
      if (verbosity>1) std::cout << " " << mean;
      if (verbosity>2) std::cout << " " << sigma;
      if (verbosity>3) std::cout << " " << skew;
      std::cout << "\n";
    }
    SplitGaussian Sub(double b1, double b2) const;
    std::vector<float> Waveform(double b1, double b2, double spacing) const;
    
    
  private:
    double c1, c2, c3;
};

#endif
