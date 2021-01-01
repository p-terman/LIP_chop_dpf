#include <iostream>
#include <string>
#include "LUX_rq1_pulse.h"
#include "LUX_rq1_header.h"
#include "LUX_rq1_channel.h"
#include "LUX_rq1_event.h"

using namespace std;

ClassImp(LUX_rq1_header)

LUX_rq1_header::LUX_rq1_header(){
  dataset_name = "";
  first_evt_in_file = 0;
  nb_evts_in_file = 0;
}

LUX_rq1_header::LUX_rq1_header(string a, UInt_t b, UInt_t c){
  dataset_name = a;
  first_evt_in_file = b;
  nb_evts_in_file = c;
}

LUX_rq1_header::LUX_rq1_header(LUX_rq1_header *lrq1h){
  dataset_name = lrq1h->Get_dataset_name();
  first_evt_in_file = lrq1h->Get_first_evt_in_file();
  nb_evts_in_file = lrq1h->Get_nb_evts_in_file();
}

LUX_rq1_header::~LUX_rq1_header(){}

void LUX_rq1_header::Fill(string a, UInt_t b, UInt_t c){
  dataset_name = a;
  first_evt_in_file = b;
  nb_evts_in_file = c;
}

void LUX_rq1_header::Fill(string str){
  dataset_name = str.substr(0,19).c_str();
  first_evt_in_file = atoi(str.substr(19,22).c_str());
  nb_evts_in_file = atoi(str.substr(23,26).c_str());
}

void LUX_rq1_header::Fill(Char_t** a, UInt_t* b, Double_t* c){
  dataset_name = a[0];
  first_evt_in_file = b[0];
  nb_evts_in_file = b[1];
}

void LUX_rq1_header::Clear(const Option_t* option){
  dataset_name = "";
  first_evt_in_file = 0;
  nb_evts_in_file = 0;
}
