#ifndef LUX_EVT_PULSE_H
#define LUX_EVT_PULSE_H
#include "TObject.h"
#include "TClonesArray.h"
#include "TH1I.h"

#include <iostream>
using namespace std;

class Lux_EVT_Event;
class Lux_EVT_Channel;
class Lux_EVT_Header;

class Lux_EVT_Pulse : public TObject{  // the class inherits from TObject for ROOT compatability
 public:
  Lux_EVT_Pulse();               				// empty constructor
  Lux_EVT_Pulse(ULong_t, Long_t, ULong_t, ULong_t, ULong_t, TH1I*);	//constructor that fills all class members 
  Lux_EVT_Pulse(Lux_EVT_Pulse *p);				//constructor that fills all class members from another Lux_EVT_Pulse instance p
  ~Lux_EVT_Pulse();						// destructor

  void Set_pulseNumber(ULong_t x){pulseNumber=x;};		// Sets a value x to pulseNumber
  void Set_pulseStarts(Long_t x){pulseStarts=x;};		// Sets a value x to pulseStarts
  void Set_pulseLength(ULong_t x){pulseLength=x;};		// Sets a value x to pulseLength
  void Set_pulseBaseline(ULong_t x){pulseBaseline=x;};		// Sets a value x to pulseBaseline
  void Set_channelNumber(ULong_t x){channelNumber=x;};		// Sets a value x to channelNumber
  void Set_pulseData(TH1I* x){pulseData=x;}; 			// Sets a value x to pulseData  

  ULong_t Get_pulseNumber(){return pulseNumber;};		// Get the value of pulseNumber
  Long_t Get_pulseStarts(){return pulseStarts;};		// Get the value of pulseStarts
  ULong_t Get_pulseLength(){return pulseLength;};		// Get the value of pulseLength
  ULong_t Get_pulseBaseline(){return pulseBaseline;};		// Get the value of pulseBaseline
  ULong_t Get_channelNumber(){return channelNumber;};		// Get the value of channelNumber
  TH1I* Get_pulseData(){return pulseData;}; 			// Get the value of pulseData

  void Fill(ULong_t, Long_t, ULong_t, ULong_t, ULong_t, TH1I*); // fills all class members
  void Clear(Option_t* /*option*/ =""); 				// clears all class members
  
private:
  ULong_t pulseNumber;						// pulse number
  Long_t pulseStarts;						//number of samples relative to trigger that pulse begins
  ULong_t pulseLength;						// number of samples in the pulse (includes pulse detect pretrigger and pulse end posttrigger
  ULong_t pulseBaseline;					// baseline value of pulse as measured by the ADC (in ADC counts)
  ULong_t channelNumber;					// channel number of this pulse
  TH1I* pulseData;						// pulse data

  ClassDef(Lux_EVT_Pulse,1); 					// pulse information
};
#endif

