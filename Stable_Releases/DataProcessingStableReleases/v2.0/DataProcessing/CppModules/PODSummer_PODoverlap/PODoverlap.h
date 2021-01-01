#ifndef __POD_OVERLAP_H__
#define __POD_OVERLAP_H__ 1
/*

Author: Sergey Uvarov
Version: 1.0 created at 2/27/2013

*/
#include "../Utilities/LCvtEvent.h"
#include "../Utilities/LCvtChannel.h"
#include <vector>
#include <algorithm>
#include <iostream>

using std::vector;
using std::sort;
using std::cout;
using std::endl;
typedef std::vector<int> IVec;

//______________________________________________________________________________
void PODoverlapXe(LCvtEvent *event, IVec &tstart, IVec &tend) {
  // This function finds the start and end times of pulse by checking which 
  // PODs times overlap 
  // INPUTS:
  //   event -> cvt file event 
  // OUTPUTS:
  //   start -> pulses start time relative to the trigger 
  //   end -> pulses end time relative to the trigger   
  
  // clear return values
  tstart.clear();
  tend.clear();
  
  // first build up the POD windows in this event
  IVec start, end;
  int pod_start, pod_end, sumpod_start, sumpod_end;
  bool new_sumpod;
  
  // Loop through the pods looking for overlap 
  for (size_t c=0; c<122; c++) {
    if (!event->channels[c]) continue;
    
    for (unsigned int j=0; j<event->channels[c]->number_of_pods; j++) {
      pod_start = event->channels[c]->pod_starts[j];
      pod_end = int(event->channels[c]->pod_lengths[j])+pod_start-1;
      if (pod_start >= pod_end) continue;
      new_sumpod = true;
      
      // loop throught the sumpods  
      for (size_t k=0; k<start.size(); k++) {
        // check if pod does not overlap  
        if (pod_start > end[k] or pod_end < start[k]) continue;
        
        // adjust the sumpods limits
        if (pod_start < start[k]) start[k] = pod_start;
        if (pod_end > end[k]) end[k] = pod_end;
        new_sumpod = false;
      }
      if (new_sumpod) {
        start.push_back(pod_start);
        end.push_back(pod_end);
      }
    }
  }
  
  // time order the sumpods
  sort(start.begin(), start.end());
  sort(end.begin(), end.end());
  
  // compare the sumpods limits to themselfs to look for overlaps
  for (size_t i=0; i<start.size(); i++) {
    sumpod_start = start[i];
    sumpod_end = end[i];
    new_sumpod = true;
    
    // loop through a 
    for (size_t k=0; k<tstart.size(); k++) {
      // check if sumpod does not overlaps with another sumpod 
      if (sumpod_start > tend[k] or sumpod_end < tstart[k]) continue;
      
      // keep only the early sumpods 
      if (tstart[k] > sumpod_start) tstart[k] = sumpod_start;
      if (tend[k] < sumpod_end) tend[k] = sumpod_end;
      new_sumpod = false;
      break;
    }
    
    if (new_sumpod) {
      tstart.push_back(sumpod_start);
      tend.push_back(sumpod_end);
    }
  }
}
//______________________________________________________________________________
void PODoverlapWa(LCvtEvent *event, IVec &tstart, IVec &tend) {
  // This function finds the start and end times of pulse by checking which 
  // PODs times overlap 
  // INPUTS:
  //   event -> cvt file event 
  // OUTPUTS:
  //   start -> pulses start time relative to the trigger 
  //   end -> pulses end time relative to the trigger   
  
  // clear return values
  tstart.clear();
  tend.clear();
  
  // first build up the POD windows in this event
  IVec start, end;
  int pod_start, pod_end, sumpod_start, sumpod_end;
  bool new_sumpod;
  
  // Loop through the pods looking for overlap 
  for (size_t c=128; c<136; c++) {    
    if (!event->channels[c]) continue;
    
    for (unsigned int j=0; j<event->channels[c]->number_of_pods; j++) {
      pod_start = event->channels[c]->pod_starts[j];
      pod_end = int(event->channels[c]->pod_lengths[j])+pod_start-1;
      if (pod_start >= pod_end) continue;
      new_sumpod = true;
      
      // loop throught the sumpods  
      for (size_t k=0; k<start.size(); k++) {
        // check if pod does not overlap  
        if (pod_start > end[k] or pod_end < start[k]) continue;
        
        // adjust the sumpods limits
        if (pod_start < start[k]) start[k] = pod_start;
        if (pod_end > end[k]) end[k] = pod_end;
        new_sumpod = false;
      }
      if (new_sumpod) {
        start.push_back(pod_start);
        end.push_back(pod_end);
      }
    }
  }
  
  // time order the sumpods
  sort(start.begin(), start.end());
  sort(end.begin(), end.end());
  
  // compare the sumpods limits to themselfs to look for overlaps
  for (size_t i=0; i<start.size(); i++) {
    sumpod_start = start[i];
    sumpod_end = end[i];
    new_sumpod = true;
    
    // loop through a 
    for (size_t k=0; k<tstart.size(); k++) {
      // check if sumpod does not overlaps with another sumpod 
      if (sumpod_start > tend[k] or sumpod_end < tstart[k]) continue;
      
      // keep only the early sumpods 
      if (tstart[k] > sumpod_start) tstart[k] = sumpod_start;
      if (tend[k] < sumpod_end) tend[k] = sumpod_end;
      new_sumpod = false;
      break;
    }
    
    if (new_sumpod) {
      tstart.push_back(sumpod_start);
      tend.push_back(sumpod_end);
    }
  }
}
#endif
