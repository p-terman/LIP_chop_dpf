#ifndef __LUX_EVT_EVENT__
#define __LUX_EVT_EVENT__
#include "TObject.h"
#include "TClonesArray.h"

class Lux_EVT_Channel;
class Lux_EVT_Pulse;

class Lux_EVT_Event : public TObject{  // the class inherits from TObject for ROOT compatability

 public :

  Lux_EVT_Event();                  				//empty constructor
  Lux_EVT_Event(ULong_t, ULong_t, ULong_t, ULong_t, ULong_t, ULong_t, ULong_t, ULong_t, ULong64_t); // constructor that fills all class members excluding channel and pulse (these are added seperately)
  Lux_EVT_Event(Lux_EVT_Event *e);  				//constructor that fills all class members (including channel and pulse) from another Lux_EVT_Event instance e
  ~Lux_EVT_Event();   						//destructor

  void Set_eDateTime(ULong_t x){eDateTime = x;};       		// Sets a value x to eDateTime
  void Set_eLocation(ULong_t x){eLocation = x;};       		// Sets a value x to eLocation
  void Set_eventNum(ULong_t x){eventNum = x;};         		// Sets a value x to eventNum
  void Set_numChannels(ULong_t x){numChannels = x;};   		// Sets a value x to numChannels
  void Set_numPulses(ULong_t x){numPulses = x;};       		// Sets a value x to numPulses
  void Set_eventSize(ULong_t x){eventSize = x;};       		// Sets a value x to eventSize
  void Set_recordForm(ULong_t x){recordForm = x;};     		// Sets a value x to recordForm
  void Set_recordSize(ULong_t x){recordSize = x;};     		// Sets a value x to recordSize
  void Set_timeStamp(ULong64_t x){timeStamp = x;};     		// Sets a value x to timeStamp

  ULong_t Get_eDateTime(){return eDateTime;};          		// Get the value of eDateTime
  ULong_t Get_eLocation(){return eLocation ;};         		// Get the value of eLocation
  ULong_t Get_eventNum(){return eventNum;};            		// Get the value of eventNum
  ULong_t Get_numChannels(){return numChannels;};      		// Get the value of numChannels
  ULong_t Get_numPulses(){return numPulses;};          		// Get the value of numPulses
  ULong_t Get_eventSize(){return eventSize;};          		// Get the value of eventSize
  ULong_t Get_recordForm(){return recordForm;};        		// Get the value of recordForm
  ULong_t Get_recordSize(){return recordSize;};        		// Get the value of recordSize
  ULong64_t Get_timeStamp(){return timeStamp;};        		// Get the value of timeStamp
  TClonesArray* Get_channel(){return channel;};        		// Get the TClonesArray containig all channels for this event
  TClonesArray* Get_pulse(){return pulse;};            		// Get the TClonesArray containig all pulses for this event
  Lux_EVT_Channel* Get_channel(ULong_t index);         		// Get the Lux_EVT_Channel instance at <index> in the channel TClonesArray
  Lux_EVT_Pulse* Get_pulse(ULong_t index);             		// Get the Lux_EVT_Pulse instance at <index> in the pulse TClonesArray
  
  void Fill(ULong_t edt, ULong_t el,ULong_t en, ULong_t nc, ULong_t np, ULong_t es, ULong_t rf, ULong_t rs, ULong64_t ts);  // fills all class members excluding channel and pulse (these are added seperately)
  void AddChannel(Lux_EVT_Channel *t, Int_t index);    		// adds a Lux_EVT_Channel to the channel TClonesArray at <index>
  void AddPulse(Lux_EVT_Pulse *p, Int_t index);        		// adds a Lux_EVT_Pulse to the pulse TClonesArray at <index>
  void Clear(Option_t* /*option*/ ="");                		// Clears all members of this Lux_EVT_Event instance
 

 private :

  ULong_t eDateTime;                                   		// date and time of the run
  ULong_t eLocation;                                   		// location, experiment, run version
  ULong_t eventNum;                                    		// event number
  ULong_t numChannels;                                 		// number of channels in the event
  ULong_t numPulses;                                   		// number of pulses in the event (over all channels)
  ULong_t eventSize;                                   		// event size in bytes
  ULong_t recordForm;                                  		// 
  ULong_t recordSize;                                  		// 
  ULong64_t timeStamp;                                 		// 
  TClonesArray *channel;                               		// array of all the channels in this event
  TClonesArray *pulse;                                 		// array of all pulses in this event

  ClassDef(Lux_EVT_Event,1);                           		// Event Information
};
#endif
