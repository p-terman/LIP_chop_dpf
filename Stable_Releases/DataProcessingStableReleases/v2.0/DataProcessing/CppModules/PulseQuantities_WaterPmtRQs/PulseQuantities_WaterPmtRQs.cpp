#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <map>

#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>

#include "../Utilities/LEvtFile.h"
#include "../Utilities/LCvtChannel.h"
#include "../Utilities/LCvtEvent.h"
#include "../Utilities/LCvtFile.h"
#include "../Utilities/RQFile_IO.hh"

#include "PODoverlap.h"
#include "Timing.h"
#include "Quantities.h"

using namespace std;
using boost::property_tree::ptree;

typedef vector<float> FVec;


//______________________________________________________________________________
double RawArea(FVec* v) {
  double area = 0;
  for (size_t i=0; i<v->size(); i++) area += v->at(i);
  return area;
}
//______________________________________________________________________________
vector<size_t> NBiggest(LCvtChannel *ch, size_t keep) {
  vector<size_t> index;
  map<double, size_t> area_ordered;
  for (size_t i=0; i<ch->pod_data_phe.size(); i++) {
    area_ordered[ RawArea(&ch->pod_data_phe[i]) ] = i;
  }
  map<double, size_t>::reverse_iterator rit = area_ordered.rbegin();
  size_t counter=0;
  for (;rit != area_ordered.rend(); rit++) {
    index.push_back(rit->second);
    counter++;
    if (counter >= keep) break;
  }  
  return index;
}
//______________________________________________________________________________
double FindPeakStd(LCvtChannel *ch, int pstart, int pend) {
  int pod_start, pod_end, npods = 0;
  double pod_std, avg_pod_std = 0;
  for (size_t pp=0; pp<ch->pod_data_phe.size(); pp++) {
    pod_start = ch->pod_starts[pp];
    pod_end = pod_start+int(ch->pod_lengths[pp])-1;
    if (pod_start > pend or pod_end < pstart) continue;
    // find the standard deviation of a pod 
    pod_std = 0;
    for (size_t i=0; i<24; i++) 
      pod_std += ch->pod_data_phe[pp][i]*ch->pod_data_phe[pp][i];
    avg_pod_std += pod_std/24;
    npods++;
  }
  if (npods) avg_pod_std = sqrt(avg_pod_std/npods);
  return avg_pod_std;
}
//______________________________________________________________________________
void FindModule(ptree &full_tree, ptree& module, string module_name) {  
  string dp = "data_processing_settings";
  string tag = "module_name";
  BOOST_FOREACH(ptree::value_type &v, full_tree.get_child(dp)) {
    if (v.second.count(tag)==1) {
      if (v.second.get<string>(tag)==module_name) {
        module = v.second;
      }
    }
  }  
}
//______________________________________________________________________________
int main(int argc, char* argv[]) {
  
  //parse inputs_______________________________________________________________
  if( argc !=8  ) {
      cerr << "Please specify:\n-evt_filename, evt_dir, rq_filename, rq_dir, module no., dp_xml_full_path, iq_xml_full_path" << endl;
      return -1;
  }
  
  string evt_dir = string(argv[2]);
  string evt_filename = string(argv[1]);
  string rq_dir = string(argv[4]);
  string rq_filename =string(argv[3]);
  string dp_xml_full_path = string(argv[6]);
  string iq_xml_full_path = string(argv[7]);
  
  // Read input xml files into BOOST ptrees____________________________________
  ptree dp_settings;
  
  ifstream dp_settings_ifs(dp_xml_full_path.c_str());
  read_xml(dp_settings_ifs,dp_settings);
  dp_settings_ifs.close();
  
  // Fill variables from xml __________________________________________________
  ptree module_settings;
  FindModule(dp_settings, module_settings, "PulseQuantities_WaterPmtRQs");
  size_t wsumpod_dim=module_settings.get<size_t>("parameters.max_num_sumpods");
  string wsumpod_dim_str=module_settings.get<string>("parameters.max_num_sumpods");
  double daq_sat=module_settings.get<double>("parameters.daq_saturation_mV");
  
  //Load data into LCvtFile object_____________________________________________
  string cvt_filename(evt_filename);
  LEvtFile *evtfile = new LEvtFile(evt_dir, evt_filename);
  LCvtFile* cvt = new LCvtFile(evtfile);
  delete evtfile;
  
  vector<double> mVns_per_phe(136, 69.4/25.);  
  cvt->Convert(mVns_per_phe);
  
  /*
  cvt_filename.replace(cvt_filename.rfind(".evt"), 4, ".cvt");
  {
    // This code quickly checks if a cvt file exist 
    string name = evt_dir + "/" + cvt_filename;
    ifstream tmp(name.c_str());
    if (!tmp.good()) {
      cerr << "Error: Cannot find cvt file.\n";
      tmp.close();
      return 1;
    }
    tmp.close();
  }  
  LCvtFile* cvt = new LCvtFile(evt_dir, cvt_filename);
  */
  //Load data into RQFileIO object_____________________________________________  
  RQFileIO *rqio = new RQFileIO;
  rqio->ReadFile(rq_dir+"/"+rq_filename);
  
  // Add the RQs that this module is going to generate__________________________
  if (!rqio->events.Get("wnsumpods"))
    rqio->events.AddRQ("wnsumpods", "uint32", "1");
  
  if (!rqio->events.Get("wsumpod_start_samples"))
    rqio->events.AddRQ("wsumpod_start_samples", "int32", wsumpod_dim_str);
  if (!rqio->events.Get("wsumpod_end_samples"))
    rqio->events.AddRQ("wsumpod_end_samples", "int32", wsumpod_dim_str);
  if (!rqio->events.Get("wsumpod_t0_samples"))
    rqio->events.AddRQ("wsumpod_t0_samples", "double", wsumpod_dim_str);
  if (!rqio->events.Get("wsumpod_t10l_samples"))
    rqio->events.AddRQ("wsumpod_t10l_samples", "double", wsumpod_dim_str);
  if (!rqio->events.Get("wsumpod_t50l_samples"))
    rqio->events.AddRQ("wsumpod_t50l_samples", "double", wsumpod_dim_str);
  if (!rqio->events.Get("wsumpod_t1_samples"))
    rqio->events.AddRQ("wsumpod_t1_samples", "double", wsumpod_dim_str);
  if (!rqio->events.Get("wsumpod_t50r_samples"))
    rqio->events.AddRQ("wsumpod_t50r_samples", "double", wsumpod_dim_str);
  if (!rqio->events.Get("wsumpod_t10r_samples"))
    rqio->events.AddRQ("wsumpod_t10r_samples", "double", wsumpod_dim_str);
  if (!rqio->events.Get("wsumpod_t2_samples"))
    rqio->events.AddRQ("wsumpod_t2_samples", "double", wsumpod_dim_str);
  if (!rqio->events.Get("wsumpod_std_phe_per_sample"))
    rqio->events.AddRQ("wsumpod_std_phe_per_sample", "double", wsumpod_dim_str);
  if (!rqio->events.Get("wsumpod_area_phe"))
    rqio->events.AddRQ("wsumpod_area_phe", "double", wsumpod_dim_str);
  if (!rqio->events.Get("wsumpod_height_phe_per_sample"))
    rqio->events.AddRQ("wsumpod_height_phe_per_sample", "double", wsumpod_dim_str);
  if (!rqio->events.Get("wsumpod_positive_area_phe"))
    rqio->events.AddRQ("wsumpod_positive_area_phe", "double", wsumpod_dim_str);
  if (!rqio->events.Get("wsumpod_negative_area_phe"))
    rqio->events.AddRQ("wsumpod_negative_area_phe", "double", wsumpod_dim_str);
  
  if (!rqio->events.Get("wpeak_t0_samples"))
    rqio->events.AddRQ("wpeak_t0_samples", "double", wsumpod_dim_str+",8");
  if (!rqio->events.Get("wpeak_t10l_samples"))
    rqio->events.AddRQ("wpeak_t10l_samples", "double", wsumpod_dim_str+",8");
  if (!rqio->events.Get("wpeak_t50l_samples"))
    rqio->events.AddRQ("wpeak_t50l_samples", "double", wsumpod_dim_str+",8");
  if (!rqio->events.Get("wpeak_t1_samples"))
    rqio->events.AddRQ("wpeak_t1_samples", "double", wsumpod_dim_str+",8");
  if (!rqio->events.Get("wpeak_t50r_samples"))
    rqio->events.AddRQ("wpeak_t50r_samples", "double", wsumpod_dim_str+",8");
  if (!rqio->events.Get("wpeak_t10r_samples"))
    rqio->events.AddRQ("wpeak_t10r_samples", "double", wsumpod_dim_str+",8");
  if (!rqio->events.Get("wpeak_t2_samples"))
    rqio->events.AddRQ("wpeak_t2_samples", "double", wsumpod_dim_str+",8");
  if (!rqio->events.Get("wpeak_std_phe_per_sample"))
    rqio->events.AddRQ("wpeak_std_phe_per_sample", "double", wsumpod_dim_str+",8");
  if (!rqio->events.Get("wpeak_area_phe"))
    rqio->events.AddRQ("wpeak_area_phe", "double", wsumpod_dim_str+",8");
  if (!rqio->events.Get("wpeak_height_phe_per_sample"))
    rqio->events.AddRQ("wpeak_height_phe_per_sample", "double", wsumpod_dim_str+",8");
  
  if (!rqio->events.Get("wdaq_saturation_flag"))
    rqio->events.AddRQ("wdaq_saturation_flag", "bool", "8");
  
  // grap rqs 
  Rq* nsumpod = rqio->events.Get("wnsumpods");
  Rq* wsstart = rqio->events.Get("wsumpod_start_samples");
  Rq* wsend = rqio->events.Get("wsumpod_end_samples");
  Rq* wst0 = rqio->events.Get("wsumpod_t0_samples");
  Rq* wst10l = rqio->events.Get("wsumpod_t10l_samples");
  Rq* wst50l = rqio->events.Get("wsumpod_t50l_samples");
  Rq* wst1 = rqio->events.Get("wsumpod_t1_samples");
  Rq* wst50r = rqio->events.Get("wsumpod_t50r_samples");
  Rq* wst10r = rqio->events.Get("wsumpod_t10r_samples");
  Rq* wst2 = rqio->events.Get("wsumpod_t2_samples");
  Rq* wsstd = rqio->events.Get("wsumpod_std_phe_per_sample");
  Rq* wsarea = rqio->events.Get("wsumpod_area_phe");
  Rq* wsheight = rqio->events.Get("wsumpod_height_phe_per_sample");
  Rq* wsparea = rqio->events.Get("wsumpod_positive_area_phe");
  Rq* wsnarea = rqio->events.Get("wsumpod_negative_area_phe");
  Rq* wpt0 = rqio->events.Get("wpeak_t0_samples");
  Rq* wpt10l = rqio->events.Get("wpeak_t10l_samples");
  Rq* wpt50l = rqio->events.Get("wpeak_t50l_samples");
  Rq* wpt1 = rqio->events.Get("wpeak_t1_samples");
  Rq* wpt50r = rqio->events.Get("wpeak_t50r_samples");
  Rq* wpt10r = rqio->events.Get("wpeak_t10r_samples");
  Rq* wpt2 = rqio->events.Get("wpeak_t2_samples");
  Rq* wpstd = rqio->events.Get("wpeak_std_phe_per_sample");
  Rq* wparea = rqio->events.Get("wpeak_area_phe");
  Rq* wpheight = rqio->events.Get("wpeak_height_phe_per_sample");
  Rq* wdaq = rqio->events.Get("wdaq_saturation_flag");
  
  //-----EVENT PROCESSING LOOP BEGINS-----
  double results[8], threshold, tmp;
  int npods, start, end;
  LCvtChannel *ch, *wsumpod;
  vector<float> waveform;
  vector<size_t> index;
  vector<int> pstart, pend;
  for (size_t e=0; e<cvt->cvt_events.size(); e++) {
    rqio->GetEvent(e);     
    
    if (cvt->cvt_events[e]->channels.size() < 138) {
    	PODoverlapWa(cvt->cvt_events[e], pstart, pend);
    	cvt->cvt_events[e]->AddWaSumPodCh(pstart, pend);
    }
    wsumpod = cvt->cvt_events[e]->channels.at(137);
    
    // do the event level rqs 
    npods = wsumpod->pod_data_phe.size();
    nsumpod->SetInt(npods);
    
    // get the only the nbiggest pulses 
    index = NBiggest(wsumpod, wsumpod_dim);
    
    for (size_t p=0; p<index.size(); p++) {
      
      // do the sumpod level rqs 
      start = wsumpod->pod_starts[index[p]];
      end = start+int(wsumpod->pod_lengths[index[p]])-1;
      wsstart->SetInt(start, p);
      wsend->SetInt(end, p);
      // get the waveform 
      waveform = wsumpod->GetPeak(start, end);  
      
      // do the timing
      threshold = 0;
      for (size_t c=0; c<8; c++) {
        ch = cvt->cvt_events[e]->channels.at(c+128);
        tmp = FindPeakStd(ch, start, end);
        threshold += tmp*tmp;
      }
      threshold = sqrt(threshold);  
      HeightTiming(&waveform, threshold, results);
      wst0->SetDouble(results[0], p);
      wst10l->SetDouble(results[1], p);
      wst50l->SetDouble(results[2], p);
      wst1->SetDouble(results[3], p);
      wst50r->SetDouble(results[4], p);
      wst10r->SetDouble(results[5], p);
      wst2->SetDouble(results[6], p);
      wsstd->SetDouble(threshold, p);
      
      // do the basic quantities 
      BasicQuantites(&waveform, results[0], results[6], results);
      wsarea->SetDouble(results[0], p);
      wsheight->SetDouble(results[1], p);
      wsparea->SetDouble(results[2], p);
      wsnarea->SetDouble(results[3], p);
      
      // calculate the peak level rqs 
      for (int pp=0; pp<8; pp++) {
        // grab the water pmt channel 
        ch = cvt->cvt_events[e]->channels.at(pp+128);
        waveform = ch->GetPeak(start, end);
        
        // do the peak timing information
        threshold = FindPeakStd(ch, start, end);        
        HeightTiming(&waveform, threshold, results);
        wpt0->SetDouble(results[0], p, pp);
        wpt10l->SetDouble(results[1], p, pp);
        wpt50l->SetDouble(results[2], p, pp);
        wpt1->SetDouble(results[3], p, pp);
        wpt50r->SetDouble(results[4], p, pp);
        wpt10r->SetDouble(results[5], p, pp);
        wpt2->SetDouble(results[6], p, pp);
        wpstd->SetDouble(threshold, p, pp);
        
        // do the area and height 
        BasicQuantites(&waveform, results[0], results[6], results);
        wparea->SetDouble(results[0], p, pp);
        wpheight->SetDouble(results[1], p, pp);
        
        // do the daq saturation flag 
        if (results[1] >= daq_sat) {
          wdaq->SetInt(1, pp);
        }
      }      
    }
  }
  //-----EVENT PROCESSING LOOP ENDS-----
  
  // Write the output rqs
  rqio->WriteFile(rq_dir+"/"+rq_filename);

  //clean up___________________________________________________________________
  delete cvt;
  delete rqio;

  return 0;
}
