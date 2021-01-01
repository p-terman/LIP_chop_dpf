#ifndef __MONTECARLO_ANSWERKEY_H__
#define __MONTECARLO_ANSWERKEY_H__ 1

#include <vector>
#include <iostream>
#include <cmath>
#include "../Utilities/LEvtFile.h"

#define binaryread(x,b) read((char*)&x,sizeof(x)); if (b) EndianByteSwap(x)
#define binaryreadl(x,l,b) read((char*)x,l); if(b) ByteSwap((unsigned char*)x,l)

using namespace std;

using std::vector;

// This struct is used to hold an "answer key". This is all information that
// LUXSim2evt writes at the end of the .evt file that contains the MC trutch
// information. A vector of these are generated later.
struct answer_key {
  long long key_number;
  long long luxsim_random_seed;
  long long event_number_luxsim;
  short caused_trigger;
  int len_of_prim_par_type_string;
  char* prim_par_type_string;
  double prim_par_energy_keV;
  double prim_par_pos_cm[3];
  double prim_par_dir[3];
  int len_of_par_type_string;
  char* par_type_string;
  int NEST_num_gammas;
  int NEST_num_electrons;
  double event_energy_keV;
  double event_pos_cm[3];
  double missing_energy_keV;
  int number_of_scatters;
  vector<int> NEST_num_gammas_scatters;
  vector<int> NEST_num_electrons_scatters;
  vector<double> energy_deps, x_scats, y_scats, z_scats;
  double msd_cm;
  int photons_per_pmt[122];
  unsigned long long timestamp;
  long long luxsim2evt_random_seed;
  int photonID[10];
};


void GoToEndOfLastWaveform(LEvtFile* evt) {
  ifstream &input = evt->input;
  input.seekg(evt->fh_gid_byte_loc[evt->fh_gid_byte_loc.size()-1]);

  char tempchar;
  unsigned int tempint;
  unsigned short tempshort;
  unsigned long long templonglong;
  double tempdouble;
  bool byteswap = evt->fh_endianness != 0x01020304;
  input.binaryread(tempint, byteswap);     // datetime
  input.binaryread(tempint, byteswap);     // location
  input.binaryread(tempint, byteswap);     // event_number
  unsigned int num_chans;
  input.binaryread(num_chans, byteswap);   // number of channels
  input.binaryread(tempint, byteswap);     // event size in bytes;
  input.binaryread(templonglong, byteswap);// ddc trig timestamp
  input.binaryread(tempint, byteswap);     // ddc trig seq num
  input.binaryread(tempint, byteswap);     // max filter resp
  input.binaryread(tempchar, byteswap);    // max chan id
  input.binaryread(tempshort, byteswap);   // num ddc boards
  for(int board=0; board<tempshort; board++) {
    input.binaryread(tempchar, byteswap);  // s1 hit vec
    input.binaryread(tempchar, byteswap);  // s2 hit vec
  }
  input.binaryread(tempchar, byteswap);    // ddc check byte
  input.binaryread(tempint, byteswap);     // record format
  input.binaryread(tempint, byteswap);     // record size
  input.binaryread(templonglong, byteswap);// trigger timestamp

  for(unsigned int chan=0; chan<num_chans; chan++) {
    input.binaryread(tempchar, byteswap);    // binary datatype
    input.binaryread(tempdouble, byteswap);  // voltage res
    input.binaryread(tempdouble, byteswap);  // voltage offset
    input.binaryread(tempdouble, byteswap);  // time res
    input.binaryread(tempint, byteswap);     // pretrigger
    input.binaryread(tempint, byteswap);     // eventsize
    input.binaryread(tempint, byteswap);     // pulse detect pretrig
    input.binaryread(tempint, byteswap);     // pulse end posttrig
    unsigned int num_pulses;
    input.binaryread(num_pulses, byteswap);  // number of pulses
    for(unsigned int p=0; p<num_pulses; p++) {
      input.binaryread(tempint, byteswap);
    }
    vector<unsigned int> lengths;
    for(unsigned int p=0; p<num_pulses; p++) {
      unsigned int length;
      input.binaryread(length, byteswap);
      lengths.push_back(length);
    }
    for(unsigned int p=0; p<num_pulses; p++) {
      input.binaryread(tempint, byteswap);
    }
    for(unsigned int p=0; p<num_pulses; p++) {
      input.ignore(sizeof(tempshort) * lengths[p]);
    }
  }
  return;
}



void ReadAnswerKeys(LEvtFile* evt, vector<answer_key> &answer_keys) {

  ifstream &input = evt->input;
  bool byteswap = evt->fh_endianness != 0x01020304;

  unsigned int number_answer_keys;
  input.binaryread(number_answer_keys, byteswap);

  vector<long long> key_locations;
  for(unsigned int key_ind=0; key_ind < number_answer_keys+1; key_ind++) {
    long long key_location;
    input.binaryread(key_location, byteswap);
    key_locations.push_back(key_location);
  }

  for(unsigned int key_ind=0; key_ind < number_answer_keys; key_ind++) {
    answer_key key;

    input.binaryread(key.key_number, byteswap);
    input.binaryread(key.luxsim_random_seed, byteswap);
    input.binaryread(key.event_number_luxsim, byteswap);
    input.binaryread(key.caused_trigger , byteswap);
    input.binaryread(key.len_of_prim_par_type_string , byteswap);
    key.prim_par_type_string = new char[key.len_of_prim_par_type_string+1];
    input.binaryreadl(key.prim_par_type_string, key.len_of_prim_par_type_string , byteswap);
    key.prim_par_type_string[key.len_of_prim_par_type_string] = '\0';
    input.binaryread(key.prim_par_energy_keV, byteswap);
    input.binaryread(key.prim_par_pos_cm, byteswap);
    input.binaryread(key.prim_par_dir, byteswap);
    input.binaryread(key.len_of_par_type_string, byteswap);
    key.par_type_string = new char[key.len_of_par_type_string+1];
    input.binaryreadl(key.par_type_string, key.len_of_par_type_string , byteswap);
    key.par_type_string[key.len_of_par_type_string] = '\0';
    input.binaryread(key.NEST_num_gammas, byteswap);
    input.binaryread(key.NEST_num_electrons, byteswap);
    input.binaryread(key.event_energy_keV, byteswap);
    input.binaryread(key.event_pos_cm, byteswap);
    input.binaryread(key.missing_energy_keV, byteswap);
    input.binaryread(key.number_of_scatters, byteswap);

    for(unsigned int scat_ind = 0; scat_ind<key.number_of_scatters; scat_ind++) {
      int temp_int;
      double temp_dbl;
      input.binaryread(temp_int, byteswap);   // NEST Numbers
      key.NEST_num_gammas_scatters.push_back(temp_int);
      input.binaryread(temp_int, byteswap);   // NEST Numbers
      key.NEST_num_electrons_scatters.push_back(temp_int);
      input.binaryread(temp_dbl, byteswap);
      key.energy_deps.push_back(temp_dbl);
      input.binaryread(temp_dbl, byteswap);
      key.x_scats.push_back(temp_dbl);
      input.binaryread(temp_dbl, byteswap);
      key.y_scats.push_back(temp_dbl);
      input.binaryread(temp_dbl, byteswap);
      key.z_scats.push_back(temp_dbl);
      // Read in the empty reserved space.
      double empty_field;
      for(int num_empty_fields=0; num_empty_fields<10; num_empty_fields++)
        input.binaryread(empty_field, byteswap);
    }
    input.binaryread(key.msd_cm, byteswap);
    input.binaryread(key.photons_per_pmt, byteswap);
    input.binaryread(key.timestamp, byteswap);
    input.binaryread(key.luxsim2evt_random_seed, byteswap);
    input.binaryreadl(key.photonID, 10*4, byteswap);
    // Read in the empty reserved space. 122*3*8 is was the original size of
    // the reserved empty space.
    //double empty_field;
    int empty_bytes = 122*3*8 - 8 - 4*10;
    input.ignore(empty_bytes);
    //for(int num_empty_fields=0; num_empty_fields<(122*3-1); num_empty_fields++) {
      //input.binaryread(empty_field, byteswap);
    //}
    answer_keys.push_back(key);

    /*
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
    cout << endl;
    cout << endl;
    */
  }

  return;
}

#endif
