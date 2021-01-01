#include <iostream>
#include "LUX_rq1_event.h"
#include "LUX_rq1_pulse.h"
#include "LUX_rq1_header.h"
#include "LUX_rq1_channel.h"
using namespace std;

ClassImp(LUX_rq1_pulse);

LUX_rq1_pulse::LUX_rq1_pulse(){
  t0 = 0;
  t10l = 0;
  t50l = 0;
  t1 = 0;
  t50r = 0;
  t10r = 0;
  t2 = 0;
  ft10l = 0;
  ft50l = 0;
  ft1 = 0;
  ft50r = 0;
  ft10r = 0;
  peak_height = 0;
  whichpeak = 0;
  whichpulse = 0;
  s2candidate_areadiffs = 0;
  evt_sum_per_ch = 0;
  evt_mean = 0;
  evt_std = 0;
  npts_above_thresh = 0;
  t_mean = 0;
  t_std = 0;
  prompt_fraction = 0;
  preS1 = 0;
  sat = 0;
  evt_list = 0;
  timestamp = 0;
}

void LUX_rq1_pulse::Clear(const Option_t* option){
  t0 = 0;
  t10l = 0;
  t50l = 0;
  t1 = 0;
  t50r = 0;
  t10r = 0;
  t2 = 0;
  ft10l = 0;
  ft50l = 0;
  ft1 = 0;
  ft50r = 0;
  ft10r = 0;
  peak_height = 0;
  whichpeak = 0;
  whichpulse = 0;
  s2candidate_areadiffs = 0;
  evt_sum_per_ch = 0;
  evt_mean = 0;
  evt_std = 0;
  npts_above_thresh = 0;
  t_mean = 0;
  t_std = 0;
  prompt_fraction = 0;
  preS1 = 0;
  sat = 0;
  evt_list = 0;
  timestamp = 0;
}
LUX_rq1_pulse::LUX_rq1_pulse(Double_t a0, Double_t a1, Double_t a2, Double_t a3, Double_t a4, Double_t a5, Double_t a6, Double_t a7, Double_t a8, Double_t a9, Double_t a10, Double_t a11, Double_t a12, Double_t a13, Double_t a14, Double_t a15, Double_t a16, Double_t a17, Double_t a18, Double_t a19, Double_t a20, Double_t a21, Double_t a22, Double_t a23, Double_t a24, Double_t a25, Double_t a26){
  t0 = a0;
  t10l = a1;
  t50l = a2;
  t1 = a3;
  t50r = a4;
  t10r = a5;
  t2 = a6;
  ft10l = a7;
  ft50l = a8;
  ft1 = a9;
  ft50r = a10;
  ft10r = a11;
  peak_height = a12;
  whichpeak = a13;
  whichpulse = a14;
  s2candidate_areadiffs = a15;
  evt_sum_per_ch = a16;
  evt_mean = a17;
  evt_std = a18;
  npts_above_thresh = a19;
  t_mean = a20;
  t_std = a21;
  prompt_fraction = a22;
  preS1 = a23;
  sat = a24;
  evt_list = a25;
  timestamp = a26;
}

LUX_rq1_pulse::LUX_rq1_pulse(LUX_rq1_pulse *lrq1p){
  t0 = lrq1p->Get_t0();
  t10l = lrq1p->Get_t10l();
  t50l = lrq1p->Get_t50l();
  t1 = lrq1p->Get_t1();
  t50r = lrq1p->Get_t50r();
  t10r = lrq1p->Get_t10r();
  t2 = lrq1p->Get_t2();
  ft10l = lrq1p->Get_ft10l();
  ft50l = lrq1p->Get_ft50l();
  ft1 = lrq1p->Get_ft1();
  ft50r = lrq1p->Get_ft50r();
  ft10r = lrq1p->Get_ft10r();
  peak_height = lrq1p->Get_peak_height();
  whichpeak = lrq1p->Get_whichpeak();
  whichpulse = lrq1p->Get_whichpulse();
  s2candidate_areadiffs = lrq1p->Get_s2candidate_areadiffs();
  evt_sum_per_ch = lrq1p->Get_evt_sum_per_ch();
  evt_mean = lrq1p->Get_evt_mean();
  evt_std = lrq1p->Get_evt_std();
  npts_above_thresh = lrq1p->Get_npts_above_thresh();
  t_mean = lrq1p->Get_t_mean();
  t_std = lrq1p->Get_t_std();
  prompt_fraction = lrq1p->Get_prompt_fraction();
  preS1 = lrq1p->Get_preS1();
  sat = lrq1p->Get_sat();
  evt_list = lrq1p->Get_evt_list();
  timestamp = lrq1p->Get_timestamp();
}

LUX_rq1_pulse::~LUX_rq1_pulse(){}

void LUX_rq1_pulse::Fill(Double_t a0, Double_t a1, Double_t a2, Double_t a3, Double_t a4, Double_t a5, Double_t a6, Double_t a7, Double_t a8, Double_t a9, Double_t a10, Double_t a11, Double_t a12, Double_t a13, Double_t a14, Double_t a15, Double_t a16, Double_t a17, Double_t a18, Double_t a19, Double_t a20, Double_t a21, Double_t a22, Double_t a23, Double_t a24, Double_t a25, Double_t a26){
  t0 = a0;
  t10l = a1;
  t50l = a2;
  t1 = a3;
  t50r = a4;
  t10r = a5;
  t2 = a6;
  ft10l = a7;
  ft50l = a8;
  ft1 = a9;
  ft50r = a10;
  ft10r = a11;
  peak_height = a12;
  whichpeak = a13;
  whichpulse = a14;
  s2candidate_areadiffs = a15;
  evt_sum_per_ch = a16;
  evt_mean = a17;
  evt_std = a18;
  npts_above_thresh = a19;
  t_mean = a20;
  t_std = a21;
  prompt_fraction = a22;
  preS1 = a23;
  sat = a24;
  evt_list = a25;
  timestamp = a26;
}

void LUX_rq1_pulse::Fill(Char_t**a, UInt_t* b, Double_t* c){
  t0 = c[0];
  t10l = c[1];
  t50l = c[2];
  t1 = c[3];
  t50r = c[4];
  t10r = c[5];
  t2 = c[6];
  ft10l = c[7];
  ft50l = c[8];
  ft1 = c[9];
  ft50r = c[10];
  ft10r = c[11];
  peak_height = c[12];
  whichpeak = c[13];
  whichpulse = c[14];
  s2candidate_areadiffs = c[15];
  evt_sum_per_ch = c[16];
  evt_mean = c[17];
  evt_std = c[18];
  npts_above_thresh = c[19];
  t_mean = c[20];
  t_std = c[21];
  prompt_fraction = c[22];
  preS1 = c[23];
  sat = c[24];
  evt_list = c[25];
  timestamp = c[26];
}