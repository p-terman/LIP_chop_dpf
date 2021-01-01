#ifndef LUX_RQ1_CHANNEL_H
#define LUX_RQ1_CHANNEL_H
#include "TObject.h"

#include <iostream>
using namespace std;

class LUX_rq1_event;
class LUX_rq1_header;
class LUX_rq1_pulse;
class LUX_rq1_channel : public TObject{
 public:
  LUX_rq1_channel();
  LUX_rq1_channel(LUX_rq1_channel *lrq1c);
  ~LUX_rq1_channel();
  LUX_rq1_channel(Double_t, Double_t, Double_t, Double_t, Double_t, Double_t);

  void Set_head(Double_t x){head=x;};
  void Set_tail(Double_t x){tail=x;};
  void Set_peak_area_per_ch(Double_t x){peak_area_per_ch=x;};
  void Set_peak_height_per_ch(Double_t x){peak_height_per_ch=x;};
  void Set_peak_height_mV_per_ch(Double_t x){peak_height_mV_per_ch=x;};
  void Set_baseline_mV(Double_t x){baseline_mV=x;};

  Double_t Get_head(){return head;};
  Double_t Get_tail(){return tail;};
  Double_t Get_peak_area_per_ch(){return peak_area_per_ch;};
  Double_t Get_peak_height_per_ch(){return peak_height_per_ch;};
  Double_t Get_peak_height_mV_per_ch(){return peak_height_mV_per_ch;};
  Double_t Get_baseline_mV(){return baseline_mV;};

  void Fill(Double_t, Double_t, Double_t, Double_t, Double_t, Double_t);
  void Fill(Char_t** a, UInt_t* b, Double_t* c);
  void Clear(const Option_t* option="");

 private:
  Double_t head;
  Double_t tail;
  Double_t peak_area_per_ch;
  Double_t peak_height_per_ch;
  Double_t peak_height_mV_per_ch;
  Double_t baseline_mV;
  ClassDef(LUX_rq1_channel,1);
};
#endif
