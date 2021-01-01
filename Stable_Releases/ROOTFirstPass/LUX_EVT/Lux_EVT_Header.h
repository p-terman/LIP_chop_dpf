#ifndef LUX_EVT_HEADER_H
#define LUX_EVT_HEADER_H

#include "TObject.h"
#include "TClonesArray.h"
#include "TArrayL.h"

#include <iostream>
using namespace std;

class Lux_EVT_Event;
class Lux_EVT_Pulse;
class Lux_EVT_Channel;

class Lux_EVT_Header : public TObject{  // the class inherits from TObject for ROOT compatability
public:
  Lux_EVT_Header();                                   		// empty constructor
  Lux_EVT_Header(ULong_t, ULong_t, ULong_t, ULong_t, ULong_t); 	// constructor that fills all class members
  Lux_EVT_Header(Lux_EVT_Header *h);                  		// constructor that fills all class members from another Lux_Header instance h
  ~Lux_EVT_Header();                                  		// destructor

  void Set_dateTime(ULong_t x){dateTime = x;};        		// Sets a value x to dateTime
  void Set_location(ULong_t x){location = x;};        		// Sets a value x to location
  void Set_firstEvent(ULong_t x){firstEvent = x;};    		// Sets a value x to firstEvent
  void Set_numEvents(ULong_t x){numEvents = x;};      		// Sets a value x to numEvents
  void Set_numSeq(ULong_t x){numSeq = x;};                      // Sets a value x to numSeq
  void Set_timeStampLatch(TArrayL* x){timeStampLatch = x;};     // Sets the array timeStampLatch equal to the array x  
  void Set_timeStampLatch(ULong_t x, Int_t index){if(index>0){timeStampLatch->Set(index+1);timeStampLatch->AddAt(x,index);}};  // sets a value x to timeStampLatch for sequence <index>
  void Set_timeStampEnd(TArrayL* x){timeStampEnd = x;};         // Sets the array timeStampEnd equal to the array x  
  void Set_timeStampEnd(ULong_t x, Int_t index){if(index>0){timeStampEnd->Set(index+1);timeStampEnd->AddAt(x,index);}};  // sets a value x to timeStampEnd for sequence <index>
  
  ULong_t Get_dateTime(){return dateTime;};           		// Get the value of dateTime
  ULong_t Get_location(){return location;};           		// Get the value of location
  ULong_t Get_firstEvent(){return firstEvent;};       		// Get the value of firstEvent
  ULong_t Get_numEvents(){return numEvents;};         		// Get the value of numEvents
  UShort_t Get_numSeq(){return numSeq;};			// get the value of numSeq
  TArrayL* Get_timeStampLatch(){return timeStampLatch;};        // Get the TArrayL containing timeStampLatch information
  ULong_t Get_timeStampLatch(Int_t index){return timeStampLatch->At(index);};  // Get the timeStampLatch value for sequence <index>
  TArrayL* Get_timeStampEnd(){return timeStampEnd;};            // Get the TArrayL containing timeStampEnd information
  ULong_t Get_timeStampEnd(Int_t index){return timeStampEnd->At(index);};  // Get the timeStampEnd value for sequence <index>

  void Fill(ULong_t, ULong_t, ULong_t, ULong_t, ULong_t);      	// Fills all class members
  void Clear(Option_t* /*option*/ ="");                    		// Clears all members of this Lux_EVT_Header instance

private:
  ULong_t dateTime;                                   		// Date and time of the run
  ULong_t location;                                   		// location, experiment, run version
  ULong_t firstEvent;                                 		// first event in the file
  ULong_t numEvents;                                  		// number of events in the file  
  UShort_t numSeq;						// number of sequences in the header
  TArrayL *timeStampLatch;				// array of each sequences timestamp latch  
  TArrayL *timeStampEnd;				// array of each sequences timestamp end
  
  ClassDef(Lux_EVT_Header,1);                         		// Run Header
};
#endif
