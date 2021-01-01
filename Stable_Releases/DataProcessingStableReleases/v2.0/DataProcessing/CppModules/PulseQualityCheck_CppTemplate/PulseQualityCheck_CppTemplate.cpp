#include <iostream>
#include <fstream>
#include <string>

#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>

#include "../Utilities/LEvtEvent.h"
#include "../Utilities/LEvtFile.h"
#include "../Utilities/RQFile_IO.hh"

#define VERBOSITY 0

using std::vector;
using std::cout;
using std::endl;
using std::string;
using boost::property_tree::ptree;

// =================================================================== //
// Making a new module?                                                //
//                                                                     //
// SEE LINE ~153 FOR THE BEGINNING OF WHAT YOU SHOULD EDIT             //
// =================================================================== // 

void PrintOutFileContents(RQFileIO *rqio, unsigned int n_seqs_max=1000000) {

  cout << "\nXML\n" << rqio->xml << endl;

  cout << "\nHEADER block\n" << rqio->header.header_string << endl << endl;
  for (unsigned int h=0; h<rqio->header.GetNSequences(); h++) {
    if(h>n_seqs_max) break;
    rqio->GetHeader(h);
    rqio->header.PrintRQs();
  }

  cout << "\nEVENTS block\n" << rqio->events.header_string << endl << endl;
  for (unsigned int e=0; e<rqio->events.GetNSequences(); e++) {
    if(e>n_seqs_max) break;
    rqio->GetEvent(e);
    rqio->events.PrintRQs();
  }

  cout << "\nLIVETIME block\n" << rqio->livetime.header_string << endl << endl;  
  for (unsigned int l=0; l<rqio->livetime.GetNSequences(); l++) {
    if(l>n_seqs_max) break;
    rqio->GetLivetime(l);
    rqio->livetime.PrintRQs();
  }  
}

void load_mVns_per_phe(ptree& pt, vector<double> &mVns_per_phe){
  mVns_per_phe.clear();
  mVns_per_phe.resize(122);
  for(ptree::iterator iter=pt.begin();iter!=pt.end();iter++){
    if (iter->second.count("global")==0) continue;
    if (iter->second.get_child("global").count("iq_type")==0) continue;
    if (iter->second.get_child("global").get<string>("iq_type").compare("pmt_gains")) continue;
    if (iter->second.count("fit")==0) continue;
    
    BOOST_FOREACH(ptree::value_type &v,iter->second.get_child("fit")){
      if(v.first=="channel"){
        int ind = v.second.get<int>("PMT_number")-1;
        mVns_per_phe[ind]=v.second.get<double>("mVns_per_phe");
      }
    }
  }
}
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

#if VERBOSITY>0
  cout<<evt_dir+"/"+evt_filename<<"\n"<<rq_dir+"/"+rq_filename<<"\n"<<dp_xml_full_path<<"\n"<<iq_xml_full_path<<endl;
#endif

  // Read input xml files into BOOST ptrees____________________________________
  ptree dp_settings, iqs;

  ifstream dp_settings_ifs(dp_xml_full_path.c_str());
  //ifstream dp_settings_ifs(dp_xml_full_path, ifstream::in);
  read_xml(dp_settings_ifs,dp_settings);
  dp_settings_ifs.close();

  ifstream iq_settings_ifs(iq_xml_full_path.c_str());
  read_xml(iq_settings_ifs,iqs);
  iq_settings_ifs.close();

  // Fill variables from xml __________________________________________________
  vector<double> mVns_per_phe;
  load_mVns_per_phe(iqs,mVns_per_phe);

  ptree module_settings;
  FindModule(dp_settings, module_settings, "PulseQualityCheck_CppTemplate");
  cout<<module_settings.get<string>("relative_path")<<endl;

#if VERBOSITY>0
  for(int i=0;i<122;i++){
    cout<<"PMT number"<<i+1<<":\t"<<mVns_per_phe[i]<<" mVns per phe at anode"<<endl;
  }
#endif

  //Load data into RQFileIO object_____________________________________________  
  RQFileIO *rqio = new RQFileIO;
  rqio->ReadFile(rq_dir+"/"+rq_filename);

  //if the rq file contains zero events, do nothing and quit
  if(rqio->events.GetNSequences()==0) return 0;

  //Load data into LEvtFile object_____________________________________________
  //   LEvtFile* ef = new LEvtFile(evt_dir,evt_filename);
  //   ef->CalibrateChannels(mVns_per_phe);


#if VERBOSITY>0 
  //What event are we up to?
  if(ef->fh_number_events>0) std::cout<<"***Last Event: "<<ef->fh_gid_event_num.back()<<"@"<<ef->fh_gid_byte_loc.back()<<std::endl;
#endif

#if VERBOSITY>1
  PrintOutFileContents(rqio,3);
#endif
  // =================================================================== //
  //    |                                                          |     //
  //    |     THE MEAT OF A MODULE IS BELOW THIS POINT             |     //
  //    |                                                          |     //
  //    |                                                          |     //
  //    |                                                          |     //
  //   \./                                                        \./    //
  // =================================================================== // 

  // Add the RQs that this module is going to generate_______________________
  size_t pulse_dim = dp_settings.get<int>("data_processing_settings.global.max_num_pulses");
  size_t channel_dim = 122;
  string pulse_dim_str=dp_settings.get<string>("data_processing_settings.global.max_num_pulses");
  rqio->events.AddRQ("pulse_quality_flag", "double", pulse_dim_str);
  //to do : read sizes from RQIO

  // Grab the (new and pre-existing) RQs needed for algorithm__________________ 
  Rq *number = rqio->events.Get("event_number");
  RqD *pulse_quality_flag = dynamic_cast<RqD*>(rqio->events.Get("pulse_quality_flag"));
  RqD *peak_area_phe = dynamic_cast<RqD*>(rqio->events.Get("peak_area_phe"));
  RqD *peak_area_negative_phe = dynamic_cast<RqD*>(rqio->events.Get("peak_area_negative_phe"));
  RqD *pulse_length_samples = dynamic_cast<RqD*>(rqio->events.Get("pulse_length_samples"));

  //Make concise names of some settings for the event loop
  double bipolar_yint = module_settings.get<double>("parameters.negative_v_positive_area_yint");
  double bipolar_slope = module_settings.get<double>("parameters.negative_v_positive_area_slope");
  double bipolar_thresh = module_settings.get<double>("parameters.negative_area_thresh");
  int max_length = module_settings.get<int>("parameters.maximum_pulse_samples");


  //-----EVENT PROCESSING LOOP BEGINS-----
  for (unsigned int e=0; e<rqio->events.GetNSequences(); e++) {
    rqio->GetEvent(e);
#if VERBOSITY>0
    cout<<"Processing event number: "<< number->GetInt()<<endl;
#endif

    //     if(0){
    //       LEvtEvent* ee=ef->evt_events[e];
    //       ee->Sum();
    //     }

    // Flag obviously dud pulses
    for (size_t i=0; i<pulse_dim; i++){
      double flag=1;
      double pulse_area_negative_phe=0;
      double pulse_area_positive_phe=0;
      for(size_t ii=0;ii<channel_dim;ii++){
	 pulse_area_negative_phe+=peak_area_negative_phe->GetDouble(i,ii);
	 pulse_area_positive_phe+=(peak_area_phe->GetDouble(i,ii)-peak_area_negative_phe->GetDouble(i,ii));
      }
      if( ( pulse_area_negative_phe<bipolar_thresh && pulse_area_negative_phe<(bipolar_yint+bipolar_slope*pulse_area_positive_phe)) || (*pulse_length_samples)(i)>max_length) {
#if VERBOSITY>0
        cout<<"Flagging pulse "<<i+1<<" in event "<<number->GetInt()<<endl;
        cout<<"(*pulse_area_negative_phe)(i)"<<pulse_area_negative_phe<<endl;
        cout<<"(*pulse_area_positive_phe)(i)"<<pulse_area_positive_phe<<endl;
        cout<<"(*pulse_length_samples)(i)"<<(*pulse_length_samples)(i)<<endl;
#endif
        flag=0;
      }
      pulse_quality_flag->SetDouble(flag,i);
    }
  }
  //-----EVENT PROCESSING LOOP ENDS-----

  // Write the output rqs
  rqio->WriteFile(rq_dir+"/"+rq_filename);
  //rqio->WriteFile("/media/hdv0/testData/LocalRQ1/lux10_20121219T1601_cp00055/"+rq_filename);

  //clean up___________________________________________________________________
  //   delete ef;
  delete rqio;

  return 0;
}
