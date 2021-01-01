#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>

#include "../Utilities/LEvtChannel.h"
#include "../Utilities/LEvtEvent.h"
#include "../Utilities/LEvtFile.h"
#include "../Utilities/RQFile_IO.hh"

#include "MonteCarlo_AnswerKey.h"

#define P(x) cout << #x":\t" << x << endl

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
double FindPeakStd(LEvtChannel *ch, int pstart, int pend) {
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

  // Look to see if we have a cvt file instead of a evt file. We can't count
  // on a cvt file still containing the mc_answer_key stuff.
  if(false) return -1;
  
  // Read input xml files into BOOST ptrees____________________________________
  ptree dp_settings;  
  
  ifstream dp_settings_ifs(dp_xml_full_path.c_str());
  read_xml(dp_settings_ifs,dp_settings);
  dp_settings_ifs.close();
  
  // Fill variables from xml __________________________________________________   
  ptree module_settings;
  FindModule(dp_settings, module_settings, "MonteCarlo_AnswerKey");
  
  //Load data into LEvtFile object_____________________________________________
  LEvtFile* evt = new LEvtFile(evt_dir, evt_filename);
  
  //Load data into RQFileIO object_____________________________________________  
  RQFileIO *rqio = new RQFileIO;
  rqio->ReadFile(rq_dir+"/"+rq_filename);
  
  // Add the RQs that this module is going to generate__________________________
  //size_t pulse_dim=dp_settings.get<size_t>("data_processing_settings.global.max_num_pulses");
  //string pulse_dim_str=dp_settings.get<string>("data_processing_settings.global.max_num_pulses");
 

  // Get some useful RQs
  //Rq* event_number = rqio->events.Get("event_number");

  // Add MC Truth RQs
  rqio->events.AddRQ("mc_key_index", "int", "1");
  rqio->events.AddRQ("mc_luxsim_random_seed", "int", "1");
  rqio->events.AddRQ("mc_event_number_luxsim", "int", "1");
  rqio->events.AddRQ("mc_prim_par_type", "char", "20");
  rqio->events.AddRQ("mc_prim_par_E_keV", "double", "1");
  rqio->events.AddRQ("mc_prim_par_pos_x_cm", "double", "1");
  rqio->events.AddRQ("mc_prim_par_pos_y_cm", "double", "1");
  rqio->events.AddRQ("mc_prim_par_pos_z_cm", "double", "1");
  rqio->events.AddRQ("mc_prim_par_dir_x", "double", "1");
  rqio->events.AddRQ("mc_prim_par_dir_y", "double", "1");
  rqio->events.AddRQ("mc_prim_par_dir_z", "double", "1");
  rqio->events.AddRQ("mc_par_type", "char", "20");
  rqio->events.AddRQ("mc_NEST_num_gammas", "int",  "1");
  rqio->events.AddRQ("mc_NEST_num_electrons", "int",  "1");
  rqio->events.AddRQ("mc_E_keV", "double", "1");
  rqio->events.AddRQ("mc_x_cm", "double", "1");
  rqio->events.AddRQ("mc_y_cm", "double", "1");
  rqio->events.AddRQ("mc_z_cm", "double", "1");
  rqio->events.AddRQ("mc_missing_E_keV", "double", "1");
  rqio->events.AddRQ("mc_num_scats", "int", "1");
  rqio->events.AddRQ("mc_msd_cm", "double", "1");
  rqio->events.AddRQ("mc_timestamp", "unsigned long long", "1");
  rqio->events.AddRQ("mc_luxsim2evt_random_seed", "long long", "1");
  rqio->events.AddRQ("mc_photon_id", "int", "10");

  // Get MC Truth RQs
  Rq* mc_key_index = rqio->events.Get("mc_key_index");
  Rq* mc_luxsim_random_seed = rqio->events.Get("mc_luxsim_random_seed");
  Rq* mc_event_number_luxsim = rqio->events.Get("mc_event_number_luxsim");
  Rq* mc_prim_par_type = rqio->events.Get("mc_prim_par_type");
  Rq* mc_prim_par_energy_keV = rqio->events.Get("mc_prim_par_E_keV");
  Rq* mc_prim_par_pos_x_cm = rqio->events.Get("mc_prim_par_pos_x_cm");
  Rq* mc_prim_par_pos_y_cm = rqio->events.Get("mc_prim_par_pos_y_cm");
  Rq* mc_prim_par_pos_z_cm = rqio->events.Get("mc_prim_par_pos_z_cm");
  Rq* mc_prim_par_dir_x = rqio->events.Get("mc_prim_par_dir_x");
  Rq* mc_prim_par_dir_y = rqio->events.Get("mc_prim_par_dir_y");
  Rq* mc_prim_par_dir_z = rqio->events.Get("mc_prim_par_dir_z");
  Rq* mc_par_type = rqio->events.Get("mc_par_type");
  Rq* mc_NEST_num_gammas = rqio->events.Get("mc_NEST_num_gammas");
  Rq* mc_NEST_num_electrons = rqio->events.Get("mc_NEST_num_electrons");
  Rq* mc_E_keV = rqio->events.Get("mc_E_keV");
  Rq* mc_x_cm = rqio->events.Get("mc_x_cm");
  Rq* mc_y_cm = rqio->events.Get("mc_y_cm");
  Rq* mc_z_cm = rqio->events.Get("mc_z_cm");
  Rq* mc_missing_E_keV = rqio->events.Get("mc_missing_E_keV");
  Rq* mc_num_scats = rqio->events.Get("mc_num_scats");
  Rq* mc_msd_cm = rqio->events.Get("mc_msd_cm");
  Rq* mc_timestamp = rqio->events.Get("mc_timestamp");
  Rq* mc_luxsim2evt_random_seed = rqio->events.Get("mc_luxsim2evt_random_seed");
  Rq* mc_photon_id = rqio->events.Get("mc_photon_id");

  // Answer key variables.

  // Prepare the input evt file.
  string fullpath=evt_dir+'/'+evt_filename;

  evt->input.open(fullpath.data(), std::ios_base::binary);
  GoToEndOfLastWaveform(evt);
  vector<answer_key> answer_keys;
  ReadAnswerKeys(evt, answer_keys);
  evt->input.close();

  /*
  for(int i=0; i<answer_keys.size(); i++) {
    answer_key key = answer_keys[i];
    cout << "Answer key number:    \t" << key.key_number << endl;
    cout << "Sim Random Seed:      \t" << key.luxsim_random_seed << endl;
    cout << "Sim Event Number:     \t" << key.event_number_luxsim << endl;
    cout << "Passed trigger cond:  \t" << key.caused_trigger << endl;
    cout << "Primary Part Str Size:\t" << key.len_of_prim_par_type_string << endl;
    cout << "Primary Particle Type:\t" << key.prim_par_type_string << endl;
    cout << "Prim Par Energy:      \t" << key.prim_par_energy_keV << endl;
    cout << "Prim Par X,Y,Z:       \t" << key.prim_par_pos_cm[0] << ", " << key.prim_par_pos_cm[1]<< ", " << key.prim_par_pos_cm[2] << endl;
    cout << "Prim Par Dir X,Y,Z:   \t" << key.prim_par_dir[0] << ", " << key.prim_par_dir[1]<< ", " << key.prim_par_dir[2] << endl;
    cout << "Particle Str Size:    \t" << key.len_of_par_type_string << endl;
    cout << "Particle Type:        \t" << key.par_type_string << endl;
    cout << "Event energy:         \t" << key.event_energy_keV << endl;
    cout << "X,Y,Z:                \t" << key.event_pos_cm[0] << ", " << key.event_pos_cm[1] << ", " << key.event_pos_cm[2] << endl;
    cout << "Max scatter distance: \t" << key.msd_cm << endl;
    cout << "Missing energy:       \t" << key.missing_energy_keV << endl;
    cout << "Number of scatters:   \t" << key.number_of_scatters << endl;
    for(unsigned int scat_ind = 0; scat_ind <key.number_of_scatters; scat_ind++) {
      cout << "\t" << scat_ind+1 << ") " << "\t";
      cout << key.energy_deps[scat_ind] << " keV @ ";
      cout << key.x_scats[scat_ind] << ", ";
      cout << key.y_scats[scat_ind] << ", ";
      cout << key.z_scats[scat_ind] << " cm";
      cout << endl;
    }
    cout << "Photons per channel:  \t" << endl;
    for(unsigned int pmt_ind = 0; pmt_ind < 122; pmt_ind++) {
      if(pmt_ind%5 ==0) cout << "\t" << pmt_ind+1 << ") ";
      cout << "\t" << key.photons_per_pmt[pmt_ind] << "\t";
      if(pmt_ind%5 ==4) cout << endl;
    }
    cout << endl;
    cout << "Timestamp (nonsense): \t" << key.timestamp << endl;
    cout << "LUXSim2evt Rndm Seed: \t" << key.luxsim2evt_random_seed << endl;
    cout << "Photon type counts:   \t" << key.photonID[0] << " " << key.photonID[1] << " " << key.photonID[2] << " " << key.photonID[3] << " " << key.photonID[4] << " " << key.photonID[5] << " " << key.photonID[6] << " " << key.photonID[7] << " " << key.photonID[8] << " " << key.photonID[9] << endl;
    cout << endl;
    cout << endl;
  }
  */


  // This is a sanity check. If the number of answer keys do not match the
  // number of events in the file, something has gone wrong with book keeping.
  // This needs to be realized.
  if(answer_keys.size() != evt->evt_events.size()) {
    cout << "The number of answer keys does not match the number of" << endl;
    cout << "events in this file. Please contact Michael Woods." << endl;
    cout << endl << "ABORTING" << endl;
    return -1;
  }

  //-----EVENT PROCESSING LOOP BEGINS-----
  for (size_t e=0; e<evt->evt_events.size(); e++) {
    rqio->GetEvent(e);
    answer_key key = answer_keys[e];
    
    mc_key_index->SetInt(key.key_number);
    mc_luxsim_random_seed->SetInt(key.luxsim_random_seed);
    mc_event_number_luxsim->SetInt(key.event_number_luxsim);
    mc_prim_par_type->SetString(key.prim_par_type_string);
    mc_prim_par_energy_keV->SetDouble(key.prim_par_energy_keV);
    mc_prim_par_pos_x_cm->SetDouble(key.prim_par_pos_cm[0]);
    mc_prim_par_pos_y_cm->SetDouble(key.prim_par_pos_cm[1]);
    mc_prim_par_pos_z_cm->SetDouble(key.prim_par_pos_cm[2]);
    mc_prim_par_dir_x->SetDouble(key.prim_par_dir[0]);
    mc_prim_par_dir_y->SetDouble(key.prim_par_dir[1]);
    mc_prim_par_dir_z->SetDouble(key.prim_par_dir[2]);
    mc_par_type->SetString(key.par_type_string);
    mc_NEST_num_gammas->SetInt(key.NEST_num_gammas);
    mc_NEST_num_electrons->SetInt(key.NEST_num_electrons);
    mc_E_keV->SetDouble(key.event_energy_keV);
    mc_x_cm->SetDouble(key.event_pos_cm[0]);
    mc_y_cm->SetDouble(key.event_pos_cm[1]);
    mc_z_cm->SetDouble(key.event_pos_cm[2]);
    mc_msd_cm->SetDouble(key.msd_cm);
    mc_missing_E_keV->SetDouble(key.missing_energy_keV);
    mc_num_scats->SetInt(key.number_of_scatters);
    mc_timestamp->SetInt(key.timestamp);
    mc_luxsim2evt_random_seed->SetInt(key.luxsim2evt_random_seed);
    for(int i=0; i<10; i++)  {
      mc_photon_id->SetInt(key.photonID[i], i);
    }

  }
  //-----EVENT PROCESSING LOOP ENDS-----
  
  // Write the output rqs
  rqio->WriteFile(rq_dir+"/"+rq_filename);

  //clean up___________________________________________________________________
  delete evt;
  delete rqio;

  /*
  */
  return 0;
}
