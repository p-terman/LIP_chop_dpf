#include <iostream>
#include <string>
#include <exception>

#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>

#include "../Utilities/LCvtChannel.h"
#include "../Utilities/LCvtEvent.h"
#include "../Utilities/LCvtFile.h"
#include "../Utilities/RQFile_IO.hh"

#include "PHET.h"

using namespace std;
using boost::property_tree::ptree;

//______________________________________________________________________________
void FindModule(ptree &full_tree, ptree& module, string module_name) {  
  string dp = "data_processing_settings";
  string tag = "module_name";
  BOOST_FOREACH(ptree::value_type &v, full_tree.get_child(dp)) {
    if (v.second.count(tag)==1) {
      if (v.second.get<string>(tag)==module_name){
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
  
  //string evt_dir = string(argv[2]);
  //string evt_filename = string(argv[1]);
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
  FindModule(dp_settings, module_settings, "PulseClassifier_PHET");
  double parameters[8];
  try { parameters[5] = module_settings.get<float>("parameters.ppwid_cut");}
  catch (exception const& ex) { parameters[5] = 200; }
  try { parameters[1] = module_settings.get<float>("parameters.pmt_coincidence_cut");}
  catch (exception const& ex) { parameters[1] = 2; }
  try { parameters[0] = module_settings.get<float>("parameters.s1_threshold_phe_cut");}
  catch (exception const& ex) { parameters[0] = 1.0; }
  try { parameters[2] = module_settings.get<float>("parameters.s2_threshold_phe_cut");}
  catch (exception const& ex) { parameters[2] = 50; }
  try { parameters[3] = module_settings.get<float>("parameters.sp_threshold_phe_cut");}
  catch (exception const& ex) { parameters[3] = 2.0; }
  try { parameters[4] = module_settings.get<float>("parameters.s1_width_cut");}
  catch (exception const& ex) { parameters[4] = 10.0; }
  try { parameters[6] = module_settings.get<float>("parameters.peak_area_coincidence_cut");}
  catch (exception const& ex) { parameters[6] = 0.25; }
  try { parameters[7] = module_settings.get<float>("parameters.peak_height_coincidence_cut");}
  catch (exception const& ex) { parameters[7] = 0.1; }
  
  string pulse_dim_str=dp_settings.get<string>("data_processing_settings.global.max_num_pulses");
  int pulse_dim = dp_settings.get<int>("data_processing_settings.global.max_num_pulses");
  
  //Load data into RQFileIO object ____________________________________________ 
  RQFileIO *rqio = new RQFileIO;
  rqio->ReadFile(rq_dir+"/"+rq_filename);
  
  // Add the RQs that this module is going to generate ________________________
  rqio->events.AddRQ("pulse_classification", "uint32", pulse_dim_str); 
  rqio->events.AddRQ("npeaks", "uint32", pulse_dim_str); 
  rqio->events.AddRQ("pmt_coincidence", "uint32", pulse_dim_str); 
  rqio->events.AddRQ("max_peak_area_phe", "float", pulse_dim_str); 
  
  // grap rqs _________________________________________________________________
  bool cal_npulses = false;
  Rq* npulses     = rqio->events.Get("npulses");
  if (!npulses) {
    rqio->events.AddRQ("npulses", "uint32", "1"); 
    cal_npulses=true;
    npulses = rqio->events.Get("npulses");
  }
  Rq* kept_idx    = rqio->events.Get("index_kept_sumpods");
  Rq* ptype       = rqio->events.Get("pulse_classification");
  Rq* parea       = rqio->events.Get("pulse_area_phe");
  Rq* aft05       = rqio->events.Get("aft_t05_samples");
  Rq* aft25       = rqio->events.Get("aft_t25_samples");
  Rq* aft75       = rqio->events.Get("aft_t75_samples");
  Rq* aft95       = rqio->events.Get("aft_t95_samples");
  Rq* peak_area   = rqio->events.Get("peak_area_phe");
  Rq* peak_height = rqio->events.Get("peak_height_phe_per_sample");
  Rq* npeaks      = rqio->events.Get("npeaks");
  Rq* coin        = rqio->events.Get("pmt_coincidence");
  Rq* max_peak    = rqio->events.Get("max_peak_area_phe");
  Rq* tbAsym      = rqio->events.Get("top_bottom_asymmetry");
  Rq* filter      = rqio->events.Get("s2filter_max_area_diff");
  
  // Intermediate variables 
  double values[6], area, height, max_value;
  int type, coincidence, count;
  
  //-----EVENT PROCESSING LOOP BEGINS-----
  for (size_t e=0; e<rqio->events.GetNSequences(); e++) {
    rqio->GetEvent(e);
    
    if (cal_npulses) {
      count=0;
      for (int p=0; p<pulse_dim; p++) {
        if (kept_idx->GetInt(p) == 1) count++;
      }
      npulses->SetInt(count);
    }
    
  	for (int p=0; p<npulses->GetInt(); p++) {
      // Calculate the pmt coincidence and max_peak_area_phe 
      coincidence = 0;
      max_value = -1;
      for (int pp=0; pp<122; pp++) {
      	area   = peak_area->GetDouble(p, pp);
      	height = peak_height->GetDouble(p, pp);
      	if (area == 0 && height == 0) continue;
      	if (area >= parameters[6] && height >= parameters[7]) coincidence++;
      	if (max_value < area) max_value = area;
      }
      
      // Do the pulse classification 
      values[0] = parea->GetDouble(p);
      values[1] = 0.5*((aft95->GetDouble(p)-aft05->GetDouble(p))/3.28972 
                     + (aft75->GetDouble(p)-aft25->GetDouble(p))/1.34898);
      values[2] = max_value/parea->GetDouble(p);
      values[3] = coincidence;
      values[4] = tbAsym->GetDouble(p);
      values[5] = filter->GetDouble(p);

      type = PulseClassify(values, parameters);
      
      // Store the results 
      ptype->SetInt(type, p);  
      npeaks->SetInt(coincidence, p);
      coin->SetInt(coincidence, p);
      max_peak->SetDouble(max_value, p);
  	}
  }
  //-----EVENT PROCESSING LOOP ENDS-----
  
  // Write the output rqs
  rqio->WriteFile(rq_dir+"/"+rq_filename);

  //clean up___________________________________________________________________
  delete rqio;

  return 0;
}
