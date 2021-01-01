#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <ctime>
#include <vector>

#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>

#include "../Utilities/LEvtEvent.h"
#include "../Utilities/LEvtFile.h"
#include "../Utilities/RQFile_IO.hh"


using std::vector;
using std::cout;
using std::endl;
using std::string;
using std::atoi;
using boost::property_tree::ptree;

typedef unsigned long long uint64;

//------++++++------++++++------++++++------++++++------++++++------++++++------
//  Functions Declarations
//------++++++------++++++------++++++------++++++------++++++------++++++------

// ptree functions
ptree LoadXMLFile(string filename);
ptree FindModule(ptree &fullPT, string module_name);

uint64 TimeSinceRunStart(string dataset);
string GetDatasetName(string evtFileName);
string GenerateRQSeedFilename(string evtFullName);

//------++++++------++++++------++++++------++++++------++++++------++++++------
//  Main Function
//------++++++------++++++------++++++------++++++------++++++------++++++------
int main(int argc, char* argv[]) {
  
  //parse inputs
  if( argc !=8  ) {
      cerr << "Please specify:\n-evt_filename, evt_dir, rq_filename, " 
           << "rq_dir, module no., dp_xml_full_path, iq_xml_full_path" << endl;
      return -1;
  }
  
  string evt_dir = argv[2];
  string evt_filename = argv[1];
  string rq_dir = argv[4];
  string rq_filename = argv[3];
  ptree dp_settings = LoadXMLFile(argv[6]);
  ptree module_settings = FindModule(dp_settings, "InitializeRQFile_Initialize");
  
  // check if rq seed file exist
  string seedRQfilename = GenerateRQSeedFilename(evt_dir+'/'+evt_filename);
  ifstream file(seedRQfilename.c_str());
  if (file.is_open()) {
  	// seed file exist, just copy it and move on.
  	file.close();
  	RQFileIO rqio;
  	rqio.ReadFile(seedRQfilename);
  	rqio.WriteFile(rq_dir+'/'+rq_filename);
  	return 0;
  }
  file.close();
  
  //Load data into LEvtFile object
  LEvtFile* ef = new LEvtFile(evt_dir,evt_filename);
  
  //Load data into RQFileIO object
  RQFileIO *rqio = new RQFileIO;
  
  // transfer the xml information
  rqio->xml = ef->settings_string;
  
  // transfer the evt header to the rq block 1 
  rqio->header.AddRQ("dataset_name", "char", "19");
  rqio->header.AddRQ("first_evt_in_file", "uint32", "1");
  rqio->header.AddRQ("nb_evts_in_file", "uint32", "1");
  rqio->header.SetNSequences(1);  
  rqio->GetHeader(0);
  
  // Remove the '.' in the dataset_name rq 20130717 - SU
  if (evt_filename[0] == '.') {
    evt_filename.erase(0, 1);
  }

  rqio->header.Get("dataset_name")->SetString(evt_filename.substr(0,19));
  if (ef->fh_gid_event_num.size())
    rqio->header.Get("first_evt_in_file")->SetInt(ef->fh_gid_event_num[0]);
  else
    rqio->header.Get("first_evt_in_file")->SetInt(0);
  rqio->header.Get("nb_evts_in_file")->SetInt(ef->fh_gid_event_num.size());
  rqio->GetHeader(0);

  // transfer the evt livetime to the rq block 3 
  rqio->livetime.AddRQ("livetime_latch_samples", "uint64", "1");
  rqio->livetime.AddRQ("livetime_end_samples", "uint64", "1");
  rqio->livetime.SetNSequences(ef->lt_num_sequences);
  for (unsigned short i=0; i<ef->lt_num_sequences; i++) {
    rqio->GetLivetime(i);
    rqio->livetime.Get("livetime_latch_samples")->SetInt(ef->lh_timestamp_latch[i]);
    rqio->livetime.Get("livetime_end_samples")->SetInt(ef->lh_timestamp_end[i]);
  }  
  
  // transfer the event numbers from the evt file to the rq file 
  rqio->events.AddRQ("file_number", "uint32", "1");
  rqio->events.AddRQ("event_number", "uint32", "1");
  rqio->events.AddRQ("source_filename", "char", "39");
  rqio->events.AddRQ("event_timestamp_samples", "uint64", "1");
  rqio->events.AddRQ("luxstamp_samples", "uint64", "1");
  rqio->events.AddRQ("time_since_livetime_start_samples", "uint64", "1", -1);
  rqio->events.AddRQ("time_until_livetime_end_samples", "uint64", "1", -1);
  rqio->events.SetNSequences(ef->fh_gid_event_num.size());
  
  Rq* file_number = rqio->events.Get("file_number");
  Rq* event_number = rqio->events.Get("event_number");
  Rq* source_filename = rqio->events.Get("source_filename");
  Rq* event_timestamp = rqio->events.Get("event_timestamp_samples");
  Rq* lux_stamp = rqio->events.Get("luxstamp_samples");
  Rq* time_since = rqio->events.Get("time_since_livetime_start_samples");
  Rq* time_until = rqio->events.Get("time_until_livetime_end_samples");
  
  int fn = atoi(evt_filename.substr(evt_filename.rfind("_f")+2, evt_filename.rfind("_")).c_str() );
  uint64 timestamp;
  uint64 global_timestamp = TimeSinceRunStart(evt_filename.substr(0,19));
  size_t lt=0;
  size_t last_lt = ef->lh_timestamp_latch.size();
  
  for (size_t i=0; i<ef->fh_gid_event_num.size(); i++) {
    rqio->GetEvent(i);
    file_number->SetInt(fn);
    event_number->SetInt(ef->fh_gid_event_num[i]);
    source_filename->SetString(evt_filename);
    timestamp = ef->evt_events[i]->event_trigger_timestamp;
    event_timestamp->SetInt(timestamp);
    lux_stamp->SetInt(timestamp + global_timestamp);
    
    while (ef->lh_timestamp_end[lt] < timestamp && lt < last_lt) {
      lt++;
    }
    if (lt < last_lt) {
      time_since->SetInt(timestamp - ef->lh_timestamp_latch[lt]);
      time_until->SetInt(ef->lh_timestamp_end[lt] - timestamp);
    }
  }

  // Write the output rqs
  rqio->WriteFile(rq_dir+"/"+rq_filename);
	
  //clean up
  delete ef;
  delete rqio;

  return 0;
}


//------++++++------++++++------++++++------++++++------++++++------++++++------
//  Functions Implimentations
//------++++++------++++++------++++++------++++++------++++++------++++++------

ptree LoadXMLFile(string filename) {
  ptree pt;  
  ifstream file(filename.c_str());
  read_xml(file, pt);
  file.close();
  return pt;
}

//______________________________________________________________________________
ptree FindModule(ptree &fullPT, string module_name) {
  ptree module;
  string dp = "data_processing_settings";
  string tag = "module_name";
  BOOST_FOREACH(ptree::value_type &v, fullPT.get_child(dp)) {
    if (v.second.count(tag)==1) {
      if (v.second.get<string>(tag)==module_name){
        module = v.second;
      }
    }
  }
  return module;
}

//______________________________________________________________________________
uint64 TimeSinceRunStart(string datetime) {
  // Declare and nullify time structures
  tm time1 = { 0 };
  tm time2 = { 0 };
  
  datetime.erase(0,datetime.find('_')+1);
  datetime.erase(datetime.find('T'), 1);
  
  // Obtain date numbers from input string
  unsigned int year = atoi(datetime.substr(0, 4).c_str());
  unsigned int month = atoi(datetime.substr(4, 2).c_str());
  unsigned int day = atoi(datetime.substr(6, 2).c_str());
  unsigned int time_hour = atoi(datetime.substr(8, 2).c_str());
  unsigned int time_min = atoi(datetime.substr(10, 2).c_str());
  
  // Set time for desired start date (setting our "zero")
  time1.tm_year = 2011 - 1900; // years start from 1900 in this struct
  time1.tm_mon = 0; //Jan == 0, Feb == 1 ... Dec == 11
  time1.tm_mday = 1; // 1 - 31
  time1.tm_hour = 0; // 0 - 23
  time1.tm_min = 0; //0 - 59
  time1.tm_sec = 0; // 0 - 59
  
  // Set time based on input datetime string
  time2.tm_year = year - 1900; // years start from 1900 in this struct
  time2.tm_mon = month-1; //Jan == 0, Feb == 1 ... Dec == 11
  time2.tm_mday = day; // 1 - 31
  time2.tm_hour = time_hour; // 0 - 23
  time2.tm_min = time_min; //0 - 59
  time2.tm_sec = 0; // 0 - 59

  // Convert these time structures to numerical values in seconds
  time_t secs1=mktime(&time1), secs2=mktime(&time2); // convert into numbers (seconds)
  time_t difference_in_seconds = secs2 - secs1; // time in seconds since our "zero"
  uint64 difference_in_samples = (uint64)difference_in_seconds*100000000ULL;

  return difference_in_samples;
}

//______________________________________________________________________________
string GetDatasetName(string evtFileName) {
  if (evtFileName.rfind('/')) evtFileName.erase(0, evtFileName.rfind('/')+1);
  evtFileName.erase(evtFileName.rfind("_f"));
  return evtFileName;
}

//______________________________________________________________________________
string GenerateRQSeedFilename(string evtFullName) {
	evtFullName.erase(evtFullName.rfind(".evt"));
	evtFullName.append("_seed.rq");
	return evtFullName;
}






