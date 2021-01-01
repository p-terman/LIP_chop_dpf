#ifndef __POD_OVERLAP_H__
#define __POD_OVERLAP_H__ 1
/*

Author: Sergey Uvarov
Version: 1.0 created at 2/27/2013

*/
#include "../Utilities/LEvtEvent.h"
#include <vector>
#include <algorithm>
#include <map>

using std::vector;
using std::map;
using std::sort;
typedef std::vector<int> IVec;
typedef std::vector< vector<int> > I2Vec;

//______________________________________________________________________________
void WhichPODs(LEvtEvent *event, IVec *start, IVec *end, 
               I2Vec &first_pod, I2Vec &last_pod) {
  // This function finds which PODs are part on the pulse. It does this by 
  // identifing the first and last POD in a channel contributes to the pulse.
  // Note. indexing starts at 1 not 0.
  // INPUTS:
  //   event -> evt file event 
  //   start -> pulses start time relative to the trigger 
  //   end -> pulses end time relative to the trigger  
  // OUTPUTS:
  //   first_pod -> index of the first contributing POD to the pulse 
  //   last_pod -> index of the last contributing POD to the pulse 
  
  // reset the POD indicies 
  first_pod.clear();
  last_pod.clear();
  first_pod.resize(start->size(), IVec(122, 0));
  last_pod.resize(start->size(), IVec(122, 0));
  
  // find the first_pod and last_pod 
  int pod_start;
  for (size_t i=0; i<start->size(); i++) {
    for (size_t c=0; c<122; c++) {   
      for (unsigned int j=0; j<event->channels[c]->number_of_pods; j++) {
        pod_start = event->channels[c]->pod_starts[j];
        if (pod_start >= start->at(i) and !first_pod[i][c]) {
          first_pod[i][c] = j+1;
          last_pod[i][c] = j+1;
        }
        if (pod_start < end->at(i) and first_pod[i][c]) 
          last_pod[i][c] = j+1;
      }
    }
  }  
}
//______________________________________________________________________________
void PODoverlap(LEvtEvent *event, size_t max_pulses, IVec &start, IVec &end) {
  // This function finds the start and end times of pulse by checking which 
  // PODs times overlap 
  // INPUTS:
  //   event -> evt file event 
  //   max_pulses -> the maximum number of pulses to keep 
  // OUTPUTS:
  //   start -> pulses start time relative to the trigger 
  //   end -> pulses end time relative to the trigger   
  
  // clear return values
  start.clear();
  end.clear();
  
  // first build up the POD windows in this event
  IVec tstart, tend;
  int pod_start, pod_end, sumpod_start, sumpod_end;
  bool new_sumpod;
  
  // Loop through the pods looking for overlap 
  for (size_t c=0; c<122; c++) {    
    for (unsigned int j=0; j<event->channels[c]->number_of_pods; j++) {
      pod_start = event->channels[c]->pod_starts[j];
      pod_end = int(event->channels[c]->pod_lengths[j])+pod_start-1;
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
  
  // skip finding n max pulses
  if (tstart.size() <= max_pulses) {
    start = tstart;
    end = tend;
    return;
  }
  
  // find the areas of all sumpods
  map<double, int> area_sumpod;
  double area;
  vector<double> *pod_waveform;
  for (size_t i=0; i<tstart.size(); i++) {
    area = 0;
    for (size_t c=0; c<122; c++) {    
      for (unsigned int j=0; j<event->channels[c]->number_of_pods; j++) {
        pod_start = event->channels[c]->pod_starts[j];
        pod_end = int(event->channels[c]->pod_lengths[j])+pod_start-1;
        if (pod_start > tend[i] or pod_end < tstart[i]) continue;
        pod_waveform = &event->channels[c]->pod_data_phe[j];
        for (size_t k=0; k<pod_waveform->size(); k++)
          area += pod_waveform->at(k);
      }
    }
    area_sumpod[area] = i;
  }
  
  // grab "max_pulses" number of the biggest sumpods
  map<double, int>::reverse_iterator rit = area_sumpod.rbegin();
  start.clear(); end.clear();
  for (size_t i=0; rit != area_sumpod.rend() and i<max_pulses; i++) {
    start.push_back(tstart[rit->second]);
    end.push_back(tend[rit->second]);
    rit++;
  }
  
  // time order the sumpods again
  sort(start.begin(), start.end());
  sort(end.begin(), end.end());
}
#endif
