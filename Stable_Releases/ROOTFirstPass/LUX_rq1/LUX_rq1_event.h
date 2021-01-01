// Analsysis version this class was created with
// 4.00000
// Don't change the above lines!!
#ifndef LUX_RQ1_EVENT_H
#define LUX_RQ1_EVENT_H
#include "TObject.h"
#include "TClonesArray.h"
#include <vector>
class LUX_rq1_pulse;
class LUX_rq1_header;
class LUX_rq1_channel;
class LUX_rq1_event : public TObject{

public :

  LUX_rq1_event();
  LUX_rq1_event(LUX_rq1_event *lrq1e);
  LUX_rq1_event(UInt_t a, UInt_t b);
  ~LUX_rq1_event();

  void Set_num_pulses(UInt_t x){num_pulses=x;};

  UInt_t Get_num_pulses(){return num_pulses;};
  LUX_rq1_pulse* Get_pulse(UInt_t num);

  void AddPulse(LUX_rq1_pulse *lrq1p, Int_t num);
  void Clear(const Option_t* option="");

  void Set_num_channels(UInt_t x){num_channels=x;};

  UInt_t Get_num_channels(){return num_channels;};
  LUX_rq1_channel* Get_channel(UInt_t n_pulse, UInt_t n_channel);

  void AddChannel(LUX_rq1_channel *lrq1c, Int_t num);
private :

  TClonesArray *pulse;
  UInt_t num_pulses;

  TClonesArray *channel;
  UInt_t num_channels;

  ClassDef(LUX_rq1_event,1);
};

#endif
