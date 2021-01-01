#ifndef LUX_RQ1_HEADER_H
#define LUX_RQ1_HEADER_H

#include "TObject.h"

#include <iostream>
#include <string>
using namespace std;

class LUX_rq1_event;
class LUX_rq1_pulse;
class LUX_rq1_channel;
class LUX_rq1_header : public TObject{

public:
  LUX_rq1_header();
  LUX_rq1_header(LUX_rq1_header *lrq1h);
  ~LUX_rq1_header();
  LUX_rq1_header(string,UInt_t,UInt_t);

  void Set_dataset_name(string x){dataset_name = x;};
  void Set_first_evt_in_file(UInt_t x){first_evt_in_file = x;};
  void Set_nb_evts_in_file(UInt_t x){nb_evts_in_file = x;};

  string Get_dataset_name(){return dataset_name;};
  UInt_t Get_first_evt_in_file(){return first_evt_in_file;};
  UInt_t Get_nb_evts_in_file(){return nb_evts_in_file;};

  void Fill(string,UInt_t,UInt_t);
  void Fill(string str);
  void Fill(Char_t** a, UInt_t* b, Double_t* c);
  void Clear(const Option_t* option="");

private:
  string dataset_name;
  UInt_t first_evt_in_file;
  UInt_t nb_evts_in_file;

ClassDef(LUX_rq1_header,1);
};
#endif
