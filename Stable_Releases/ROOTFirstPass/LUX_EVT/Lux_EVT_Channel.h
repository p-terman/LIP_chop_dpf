#ifndef __LUX_EVT_CHANNEl__
#define __LUX_EVT_CHANNEl__
#include "TObject.h"
#include "TClonesArray.h"

#include <iostream>
using namespace std;

class Lux_EVT_Event;

class Lux_EVT_Channel : public TObject{  			// the class inherits from TObject for ROOT compatability
 public : 
  Lux_EVT_Channel();                                        	// empty constructor
  Lux_EVT_Channel(Char_t bdt, ULong_t cn, Double_t vr, Double_t vo, Double_t tr, Long_t pt, ULong_t es, ULong_t pd, ULong_t pe, ULong_t np);                                                        	// constructor that fills all class members
  Lux_EVT_Channel(Lux_EVT_Channel *c);             	   	// constructor that fills all class members from another Lux_EVT_Channel instance c
  ~Lux_EVT_Channel(){}                                      	// destructor
  
  void Set_binDataType(Char_t x){binDataType = x;};        	// Sets a value x to binDataType  
  void Set_channelNumber(ULong_t x){channelNumber = x;};	// Sets a value x to channelNumber
  void Set_voltageRes(Double_t x){voltageRes = x;};		// Sets a value x to voltageRes
  void Set_voltageOff(Double_t x){voltageOff = x;};		// Sets a value x to voltageOff 
  void Set_timeRes(Double_t x){timeRes = x;};			// Sets a value x to timeRes 
  void Set_preTrigger(Long_t x){preTrigger = x;};		// Sets a value x to preTrigger 
  void Set_eventSize(ULong_t x){eventSize = x;};		// Sets a value x to eventSize 
  void Set_pulseDetect(ULong_t x){pulseDetect = x;};		// Sets a value x to pulseDetect 
  void Set_pulseEnd(ULong_t x){pulseEnd = x;};			// Sets a value x to pulseEnd 
  void Set_numPulses(ULong_t x){numPulses = x;};		// Sets a value x to numPulses
  
  Char_t Get_binDataType(){return binDataType;};		// Get the value of binDataType
  ULong_t Get_channelNumber(){return channelNumber;};		// Get the value of channelNumber
  Double_t Get_voltageRes(){return voltageRes;};		// Get the value of voltageRes
  Double_t Get_voltageOff(){return voltageOff;};		// Get the value of voltageOff
  Double_t Get_timeRes(){return timeRes;};			// Get the value of timeRes
  Long_t Get_preTrigger(){return preTrigger;};			// Get the value of preTrigger
  ULong_t Get_eventSize(){return eventSize;};			// Get the value of eventSize
  ULong_t Get_pulseDetect(){return pulseDetect;};		// Get the value of pulseDetect
  ULong_t Get_pulseEnd(){return pulseEnd;};			// Get the value of pulseEnd
  ULong_t Get_numPulses(){return numPulses;};			// Get the value of numPulses
  
  void Fill(Char_t bdt, ULong_t cn, Double_t vr, Double_t vo, Double_t tr, Long_t pt, ULong_t es, ULong_t pd, ULong_t pe, ULong_t np);  // fills all class members
  void Clear(Option_t* /*option*/ ="");					// clears all members of this Lux_EVT_Channel instance
  
 private : 

  Char_t binDataType;						// binary data type
  ULong_t channelNumber;					// channel number
  Double_t voltageRes;						// voltage resolution
  Double_t voltageOff;						// voltage offset
  Double_t timeRes;						// time resolution
  Long_t preTrigger;						// pretrigger
  ULong_t eventSize;						// event size in bytes
  ULong_t pulseDetect;						// pulse detect pretrigger
  ULong_t pulseEnd;						// pulse end posttrigger
  ULong_t numPulses;						// number of pulses in this channel
 
  ClassDef(Lux_EVT_Channel,1);					// channel information
};
#endif
