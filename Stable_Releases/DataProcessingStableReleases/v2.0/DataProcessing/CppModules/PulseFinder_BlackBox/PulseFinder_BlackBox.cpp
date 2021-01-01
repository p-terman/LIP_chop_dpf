#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <algorithm>
#include <vector>
#include <exception>

#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>

#include "../Utilities/LCvtChannel.h"
#include "../Utilities/LCvtEvent.h"
#include "../Utilities/LCvtFile.h"
#include "../Utilities/RQFile_IO.hh"

#include "PulseFinder.h"
#include "SplitGaussian.h"

using namespace std;
using boost::property_tree::ptree;

////////////////////////////////////////////////////////////////////////////////
// Supporting Functions
////////////////////////////////////////////////////////////////////////////////
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
void QuickStat(vector<float> *wave, float &area, float &height) {
  area = height = 0;
  for (size_t i=0; i<wave->size(); i++) {
    area += wave->at(i);
    if (height < wave->at(i)) height = wave->at(i);
  }
}
//______________________________________________________________________________
int FindCoincidence(LCvtEvent *ev, int start, int end, float area, float height) {
  LCvtChannel *ch;
  int counter = 0, pstart, pend, cstart, cend;
  float parea, pheight;
  vector<float> wave;
  for (size_t c=0; c<122; c++) {
    ch = ev->channels[c];
    for (size_t pp=0; pp<ch->pod_starts.size(); pp++) {
      // check if peak in the pulse timing region 
      pstart = ch->pod_starts[pp];
      pend = pstart + int(ch->pod_lengths[pp])-1;
      if (pstart > end or pend < start) continue;
      
      // find the height and area of the peak 
      cstart = pstart>start ? pstart : start;
      cend   = pend<end ? pend : end;
      wave = ch->GetPeak(cstart, cend);
      QuickStat(&wave, parea, pheight);
      
      // coincidence area and height requirements 
      if (parea < area or pheight < height) continue;
      
      counter++;
      break;
    }
  }
  return counter;
}
////////////////////////////////////////////////////////////////////////////////
// Main Function
////////////////////////////////////////////////////////////////////////////////
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
  FindModule(dp_settings, module_settings, "PulseFinder_BlackBox");
  cout<<module_settings.get<string>("relative_path")<<endl;  
  
  // Get parameters from module xml ___________________________________________
  PulseFinder pf;
  try {	pf.kz_samples = module_settings.get<int>("parameters.kz_width");}
  catch (exception const& ex) { pf.kz_samples = 10; }
  try {	pf.kz_iterations = module_settings.get<int>("parameters.kz_iterations");}
  catch (exception const& ex) { pf.kz_iterations = 2; }
  try {	pf.max_seperation = module_settings.get<int>("parameters.max_seperation");}
  catch (exception const& ex) { pf.max_seperation = 8; }
  try {	pf.max_threshold = module_settings.get<float>("parameters.max_thres");}
  catch (exception const& ex) { pf.max_threshold = 0.15; }
  try {	pf.min_max_ratio = module_settings.get<float>("parameters.min_max_ratio");}
  catch (exception const& ex) { pf.min_max_ratio = 0.75; }
  try {	pf.nsigma = module_settings.get<float>("parameters.nsigmas");}
  catch (exception const& ex) { pf.nsigma = 3; }
  try {	pf.se_area_mean_phe = module_settings.get<float>("parameters.se_area_mean_phe");}
  catch (exception const& ex) { pf.se_area_mean_phe = 24.47; }
  try {	pf.se_area_sigma_phe = module_settings.get<float>("parameters.se_area_sigma_phe");}
  catch (exception const& ex) { pf.se_area_sigma_phe = 5.799; }
  try {	pf.se_width_mean_samples = module_settings.get<float>("parameters.se_width_mean_samples");}
  catch (exception const& ex) { pf.se_width_mean_samples = 44.23; }
  try {	pf.se_width_sigma_samples = module_settings.get<float>("parameters.se_width_sigma_samples");}
  catch (exception const& ex) { pf.se_width_sigma_samples = 11.81; }
  try {	pf.dis_inv_ratio = module_settings.get<float>("parameters.dis_inv_ratio");}
  catch (exception const& ex) { pf.dis_inv_ratio = 0.01; }
  try {	pf.dis_scale = module_settings.get<float>("parameters.dis_scale");}
  catch (exception const& ex) { pf.dis_scale = 0.05; }
  try {	pf.dis_offset = module_settings.get<float>("parameters.dis_offset");}
  catch (exception const& ex) { pf.dis_offset = 0.5; }
  pf.Initilize();
  
  int s1coincidence, evt_coin, look_foward, keep_wide;
  float coin_area, coin_height, evt_thresh, width_cut;
  try {	s1coincidence = module_settings.get<int>("parameters.coincidence");} 
  catch (exception const& ex) { s1coincidence = 2; }
  try {	coin_area = module_settings.get<float>("parameters.coin_area_thres");}
  catch (exception const& ex) { coin_area = 0.25; }
  try {	coin_height = module_settings.get<float>("parameters.coin_height_thres");} 
  catch (exception const& ex) { coin_height = 0.1; }
  try {	evt_thresh = module_settings.get<float>("parameters.event_rq.threshold");}
  catch (exception const& ex) { evt_thresh = 15; }
  try {	evt_coin = module_settings.get<int>("parameters.event_rq.coincidence");}
  catch (exception const& ex) { evt_coin = 2; }
  try {	look_foward = module_settings.get<int>("parameters.group_adj_sumpod_samples");}
  catch (exception const& ex) { look_foward = 50; }
  try {	keep_wide = module_settings.get<int>("parameters.keep_wide_pulses");}
  catch (exception const& ex) { keep_wide = 5; }
  try {	width_cut = module_settings.get<float>("parameters.width_cut_samples");}
  catch (exception const& ex) { width_cut = 40; }
  
  int pulse_dim = dp_settings.get<int>("data_processing_settings.global.max_num_pulses");
  string pulse_dim_str=dp_settings.get<string>("data_processing_settings.global.max_num_pulses");
  
  //Load data into LCvtFile object_____________________________________________
  string cvt_filename(evt_filename);
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
  
  //Load data into RQFileIO object ____________________________________________ 
  RQFileIO *rqio = new RQFileIO;
  rqio->ReadFile(rq_dir+"/"+rq_filename);
  
  // Add the RQs that this module is going to generate ________________________
  rqio->events.AddRQ("full_evt_area_phe", "float", "1"); 
  rqio->events.AddRQ("n_samples_in_evt", "uint32", "1"); 
  rqio->events.AddRQ("total_npulses", "uint32", "1"); 
  rqio->events.AddRQ("total_nsumpods", "uint32", "1"); 
  rqio->events.AddRQ("total_npulses_thresh", "uint32", "1"); 
  rqio->events.AddRQ("total_npulses_coincidence", "uint32", "1"); 
  rqio->events.AddRQ("npulses", "uint32", "1");
  rqio->events.AddRQ("npeaks", "uint32", pulse_dim_str);
  rqio->events.AddRQ("pulse_start_samples", "int32", pulse_dim_str, -999999);
  rqio->events.AddRQ("pulse_end_samples", "int32", pulse_dim_str, -999999);
  rqio->events.AddRQ("pulse_length_samples", "uint32", pulse_dim_str);
  rqio->events.AddRQ("index_kept_sumpods", "uint32", pulse_dim_str); 
  rqio->events.AddRQ("sg_amplitude_phe_per_sample", "float", pulse_dim_str);
  rqio->events.AddRQ("sg_mode_samples", "float", pulse_dim_str);
  rqio->events.AddRQ("sg_lsigma_samples", "float", pulse_dim_str);
  rqio->events.AddRQ("sg_rsigma_samples", "float", pulse_dim_str);
  rqio->events.AddRQ("sg_mean_samples", "float", pulse_dim_str);
  rqio->events.AddRQ("sg_sigma_samples", "float", pulse_dim_str);
  rqio->events.AddRQ("sg_skew3_samples", "float", pulse_dim_str);
  
  // grap rqs 
  Rq* npulses = rqio->events.Get("npulses");
  Rq* npeaks = rqio->events.Get("npeaks");
  Rq* pulse_start_samples = rqio->events.Get("pulse_start_samples");
  Rq* pulse_end_samples = rqio->events.Get("pulse_end_samples");
  Rq* pulse_length_samples = rqio->events.Get("pulse_length_samples");
  Rq* index_kept_sumpods = rqio->events.Get("index_kept_sumpods");
  Rq* event_area_phe = rqio->events.Get("full_evt_area_phe");
  Rq* total_npulses = rqio->events.Get("total_npulses");
  Rq* total_nsumpods = rqio->events.Get("total_nsumpods");
  Rq* total_npulses_thresh = rqio->events.Get("total_npulses_thresh");
  Rq* total_npulses_coincidence = rqio->events.Get("total_npulses_coincidence");
  Rq* nsamples = rqio->events.Get("n_samples_in_evt");
  
  Rq* sg_ampl = rqio->events.Get("sg_amplitude_phe_per_sample");
  Rq* sg_mode = rqio->events.Get("sg_mode_samples");
  Rq* sg_lsigma = rqio->events.Get("sg_lsigma_samples");
  Rq* sg_rsigma = rqio->events.Get("sg_rsigma_samples");
  Rq* sg_mean = rqio->events.Get("sg_mean_samples");
  Rq* sg_sigma = rqio->events.Get("sg_sigma_samples");
  Rq* sg_skew3 = rqio->events.Get("sg_skew3_samples");
  
  // Intermediate variables 
  vector<SplitGaussian> sgaus;
  SplitGaussian g;
  vector<int> pstart, pend, coincidence;
  map<int, int> time_ordered;
  map<int, int>::iterator it;
  map<double, int> area_ordered, height_ordered;
  map<double, int>::reverse_iterator rit;
  double total_area, width;
  int pulses, counter, coin, pulses_above_thresh, pulses_with_coin, keep_skiny;
  LCvtChannel *sumpod;
  int cur_spod_start, cur_spod_end, next_spod_start, next_spod_end;
  int last_wide=0, total_length;
  size_t cur_spod_idx;
  
  vector<double> width_map;
  
  //-----EVENT PROCESSING LOOP BEGINS-----
  for (size_t e=0; e<cvt->cvt_events.size(); e++) {
    rqio->GetEvent(e);
    
    //if (cvt->cvt_events[e]->gid_event_num != 614314) continue;
    //cout << "\nevent = " << e+1 << "/" << cvt->cvt_events.size() << endl;
    // clear tracking variables
    area_ordered.clear();
    height_ordered.clear();
    time_ordered.clear();
    total_area = 0;
    total_length = 0;
    pulses = 0;
    pulses_above_thresh = 0;
    pulses_with_coin = 0;
    sumpod = cvt->cvt_events[e]->channels.at(136);
    pstart.clear();
    pend.clear();
    coincidence.clear();
    
    // Find The Pulses ________________________________________________________
    cur_spod_idx = 0;
	  while (cur_spod_idx < sumpod->pod_data_phe.size()) {

	  	// grap starting sumpod and any nearby sumpods
	  	cur_spod_start = sumpod->pod_starts[cur_spod_idx];
	  	cur_spod_end = cur_spod_start+sumpod->pod_lengths[cur_spod_idx]-1;
	  	while (cur_spod_idx+1 < sumpod->pod_data_phe.size()) {
			  next_spod_start = sumpod->pod_starts[cur_spod_idx+1];
			  next_spod_end = next_spod_start+sumpod->pod_lengths[cur_spod_idx+1]-1;
			  if (next_spod_start - cur_spod_end <= look_foward) {
			    cur_spod_end = next_spod_end;
			    cur_spod_idx++;
			  }else {
			    break;
			  }
	  	}
	  	//cout << "sumpod start & end :: " << cur_spod_start << " " << cur_spod_end << endl;
	  	
    	// do the pulse finding 
    	pf.wave = sumpod->GetPeak(cur_spod_start, cur_spod_end);
    	pf.Execute();
      
      total_length += cur_spod_end-cur_spod_start+1;
      
			// loop through the found pulses 
    	counter = 0;
    	for (size_t p=0; p<pf.borders.size(); p+=2) {
        //cout << "p" << p << " " <<pf.borders[p] << " " << pf.borders[p+1]  << " "
        //     << cur_spod_start+pf.borders[p] << " " << cur_spod_start+pf.borders[p+1] << endl;
        
        // get the pulse start and end times  
        pstart.push_back(cur_spod_start+pf.borders[p]);
        pend.push_back(cur_spod_start+pf.borders[p+1]);
        g = pf.pars[p/2];
        g.Derive();
        g.mode += cur_spod_start;
        g.mean += cur_spod_start;
        width = 1.1774*(g.lsigma+g.rsigma); // FWHM conversion
        sgaus.push_back(g);
        
        // get the pulse coincidence
        coin   = FindCoincidence(cvt->cvt_events[e], pstart[pulses], 
                  pend[pulses], coin_area, coin_height);
        coincidence.push_back(coin);
        
        // event level rq 
        total_area += g.area;
        if (g.area > evt_thresh) pulses_above_thresh++;
        if (coin >= evt_coin) pulses_with_coin++;
        
        // map pulse wide & skiny pulses 
        if (width < width_cut) area_ordered[g.area] = pulses;
        else height_ordered[g.amplitude] = pulses;
        
        counter++;
        pulses++;
    	}
    	cur_spod_idx++;
    }
    
    // Store event level rq information 
    event_area_phe->SetDouble(total_area);
    nsamples->SetInt(total_length);
    total_nsumpods->SetInt(sumpod->pod_starts.size());
    total_npulses->SetInt(pulses);
    total_npulses_thresh->SetInt(pulses_above_thresh);
    total_npulses_coincidence->SetInt(pulses_with_coin);
    
    // Choose Which Pulses Are Recorded _______________________________________
    // select the wide pulses to keep 
    //cout << "wide pulse selection\n";
    counter = 0;
    last_wide = pulses;
    for (rit=height_ordered.rbegin(); rit!=height_ordered.rend(); rit++) {
      if (counter >= keep_wide) break;
      time_ordered[rit->second] = rit->second;
      if (counter==0 || last_wide<rit->second) last_wide = rit->second;
      counter++;
      //cout << rit->second << " last=" << last_wide << endl;
    }
    keep_skiny = pulse_dim - counter; // fill avalible slots with skinny pulses
    
    // select the skiny coincidence pulses to keep 
    //cout << "skiny coincidence pulse selection\n";
    counter = 0;
    for (rit=area_ordered.rbegin(); rit!=area_ordered.rend(); rit++) {
      if (counter >= keep_skiny) break;
      
      if (coincidence[rit->second] < s1coincidence) continue; // s1 priority 
      if (last_wide < rit->second) continue; // look backwards 
      time_ordered[rit->second] = rit->second;           
      counter++;
      //cout << rit->second << " last=" << last_wide << endl;
    }
    
    // store pulse the rest that can fit 
    //cout << "the rest selection\n";
    for (int i=last_wide; i--;) {
      if (counter >= keep_skiny) break;
      
      // check if already added 
      if (time_ordered.count(i)) continue;
      time_ordered[i] = i; 
      counter++;
      //cout << i << " last=" << last_wide << endl;
    }
    //cout << "# of skiny pulses keep = " << time_ordered.size() << endl;
    
    // store the results
    counter = 0;
    for (it=time_ordered.begin(); it != time_ordered.end(); it++) {
      index_kept_sumpods->SetInt(1, counter);
      pulse_start_samples->SetInt(pstart[it->first], counter);
      pulse_end_samples->SetInt(pend[it->first], counter);
      pulse_length_samples->SetInt(pend[it->first]-pstart[it->first]+1, counter);
      sg_ampl->SetDouble(sgaus[it->first].amplitude, counter);
      sg_mode->SetDouble(sgaus[it->first].mode, counter);
      sg_lsigma->SetDouble(sgaus[it->first].lsigma, counter);
      sg_rsigma->SetDouble(sgaus[it->first].rsigma, counter);
      sg_mean->SetDouble(sgaus[it->first].mean, counter);
      sg_sigma->SetDouble(sgaus[it->first].sigma, counter);
      sg_skew3->SetDouble(sgaus[it->first].skew, counter);
      npeaks->SetInt(coincidence[it->first], counter);
      counter++;
    }
    npulses->SetInt(counter);    
  }
  //-----EVENT PROCESSING LOOP ENDS-----
  
  // Write the output rqs
  rqio->WriteFile(rq_dir+"/"+rq_filename);

  //clean up___________________________________________________________________
  delete cvt;
  delete rqio;

  return 0;
}
