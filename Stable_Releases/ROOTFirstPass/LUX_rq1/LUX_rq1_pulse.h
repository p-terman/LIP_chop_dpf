#ifndef LUX_RQ1_PULSE_H
#define LUX_RQ1_PULSE_H
#include "TObject.h"

#include <iostream>
using namespace std;

class LUX_rq1_event;
class LUX_rq1_header;
class LUX_rq1_channel;

class LUX_rq1_pulse : public TObject{
 public:
  LUX_rq1_pulse();
  LUX_rq1_pulse(Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t);
  LUX_rq1_pulse(LUX_rq1_pulse *lrq1p);
  ~LUX_rq1_pulse();

  void Set_t0(Double_t x){t0 = x;};
  void Set_t10l(Double_t x){t10l = x;};
  void Set_t50l(Double_t x){t50l = x;};
  void Set_t1(Double_t x){t1 = x;};
  void Set_t50r(Double_t x){t50r = x;};
  void Set_t10r(Double_t x){t10r = x;};
  void Set_t2(Double_t x){t2 = x;};
  void Set_ft10l(Double_t x){ft10l = x;};
  void Set_ft50l(Double_t x){ft50l = x;};
  void Set_ft1(Double_t x){ft1 = x;};
  void Set_ft50r(Double_t x){ft50r = x;};
  void Set_ft10r(Double_t x){ft10r = x;};
  void Set_peak_height(Double_t x){peak_height = x;};
  void Set_whichpeak(Double_t x){whichpeak = x;};
  void Set_whichpulse(Double_t x){whichpulse = x;};
  void Set_s2candidate_areadiffs(Double_t x){s2candidate_areadiffs = x;};
  void Set_evt_sum_per_ch(Double_t x){evt_sum_per_ch = x;};
  void Set_evt_mean(Double_t x){evt_mean = x;};
  void Set_evt_std(Double_t x){evt_std = x;};
  void Set_npts_above_thresh(Double_t x){npts_above_thresh = x;};
  void Set_t_mean(Double_t x){t_mean = x;};
  void Set_t_std(Double_t x){t_std = x;};
  void Set_prompt_fraction(Double_t x){prompt_fraction = x;};
  void Set_preS1(Double_t x){preS1 = x;};
  void Set_sat(Double_t x){sat = x;};
  void Set_evt_list(Double_t x){evt_list = x;};
  void Set_timestamp(Double_t x){timestamp = x;};

  Double_t Get_t0(){return t0;};
  Double_t Get_t10l(){return t10l;};
  Double_t Get_t50l(){return t50l;};
  Double_t Get_t1(){return t1;};
  Double_t Get_t50r(){return t50r;};
  Double_t Get_t10r(){return t10r;};
  Double_t Get_t2(){return t2;};
  Double_t Get_ft10l(){return ft10l;};
  Double_t Get_ft50l(){return ft50l;};
  Double_t Get_ft1(){return ft1;};
  Double_t Get_ft50r(){return ft50r;};
  Double_t Get_ft10r(){return ft10r;};
  Double_t Get_peak_height(){return peak_height;};
  Double_t Get_whichpeak(){return whichpeak;};
  Double_t Get_whichpulse(){return whichpulse;};
  Double_t Get_s2candidate_areadiffs(){return s2candidate_areadiffs;};
  Double_t Get_evt_sum_per_ch(){return evt_sum_per_ch;};
  Double_t Get_evt_mean(){return evt_mean;};
  Double_t Get_evt_std(){return evt_std;};
  Double_t Get_npts_above_thresh(){return npts_above_thresh;};
  Double_t Get_t_mean(){return t_mean;};
  Double_t Get_t_std(){return t_std;};
  Double_t Get_prompt_fraction(){return prompt_fraction;};
  Double_t Get_preS1(){return preS1;};
  Double_t Get_sat(){return sat;};
  Double_t Get_evt_list(){return evt_list;};
  Double_t Get_timestamp(){return timestamp;};

  void Fill(Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t);
  void Fill(Char_t** a, UInt_t* b, Double_t* c);
  void Clear(const Option_t* option="");

private : 

  Double_t t0;
  Double_t t10l;
  Double_t t50l;
  Double_t t1;
  Double_t t50r;
  Double_t t10r;
  Double_t t2;
  Double_t ft10l;
  Double_t ft50l;
  Double_t ft1;
  Double_t ft50r;
  Double_t ft10r;
  Double_t peak_height;
  Double_t whichpeak;
  Double_t whichpulse;
  Double_t s2candidate_areadiffs;
  Double_t evt_sum_per_ch;
  Double_t evt_mean;
  Double_t evt_std;
  Double_t npts_above_thresh;
  Double_t t_mean;
  Double_t t_std;
  Double_t prompt_fraction;
  Double_t preS1;
  Double_t sat;
  Double_t evt_list;
  Double_t timestamp;
  ClassDef(LUX_rq1_pulse,2);
};

#endif
