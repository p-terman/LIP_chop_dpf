#include <iostream>
#include "LUX_rq1_event.h"
#include "LUX_rq1_pulse.h"
#include "LUX_rq1_header.h"
#include "LUX_rq1_channel.h"
using namespace std;

ClassImp(LUX_rq1_channel);

LUX_rq1_channel::LUX_rq1_channel(){
  head = 0;
  tail = 0;
  peak_area_per_ch = 0;
  peak_height_per_ch = 0;
  peak_height_mV_per_ch = 0;
  baseline_mV = 0;
}

LUX_rq1_channel::LUX_rq1_channel(Double_t a, Double_t b, Double_t c, Double_t d, Double_t e, Double_t f){
  head = a;
  tail = b;
  peak_area_per_ch = c;
  peak_height_per_ch = d;
  peak_height_mV_per_ch = e;
  baseline_mV = f;
}

LUX_rq1_channel::LUX_rq1_channel(LUX_rq1_channel *lrq1c){
  head = lrq1c->Get_head();
  tail = lrq1c->Get_tail();
  peak_area_per_ch = lrq1c->Get_peak_area_per_ch();
  peak_height_per_ch = lrq1c->Get_peak_height_per_ch();
  peak_height_mV_per_ch = lrq1c->Get_peak_height_mV_per_ch();
  baseline_mV = lrq1c->Get_baseline_mV();
}

LUX_rq1_channel::~LUX_rq1_channel(){}

void LUX_rq1_channel::Fill(Double_t a, Double_t b, Double_t c, Double_t d, Double_t e, Double_t f){
  head = a;
  tail = b;
  peak_area_per_ch = c;
  peak_height_per_ch = d;
  peak_height_mV_per_ch = e;
  baseline_mV = f;
}

void LUX_rq1_channel::Fill(Char_t** a, UInt_t* b, Double_t* c){
  head = c[0];
  tail = c[1];
  peak_area_per_ch = c[2];
  peak_height_per_ch = c[3];
  peak_height_mV_per_ch = c[4];
  baseline_mV = c[5];
}

void LUX_rq1_channel::Clear(const Option_t* option){
  head = 0;
  tail = 0;
  peak_area_per_ch = 0;
  peak_height_per_ch = 0;
  peak_height_mV_per_ch = 0;
  baseline_mV = 0;
}

