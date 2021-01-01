#include <iostream>
#include <fstream>
#include <string>
#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "../Utilities/RQFile_IO.hh"
#include "CorCentroid.h"

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
  
  string evt_dir = string(argv[2]);
  string evt_filename = string(argv[1]);
  string rq_dir = string(argv[4]);
  string rq_filename =string(argv[3]);
  string dp_xml_full_path = string(argv[6]);
  string iq_xml_full_path = string(argv[7]);
 
  // Read input xml files into BOOST ptrees____________________________________
  ptree dp_settings, iqs;
 
	ifstream dp_settings_ifs(dp_xml_full_path.c_str());
	read_xml(dp_settings_ifs,dp_settings);
	dp_settings_ifs.close();
  
  ptree module_settings;
  FindModule(dp_settings, module_settings, "PositionReconstruction_CorrCentroid");  
  cout<<module_settings.get<string>("relative_path")<<endl;
  
  //Load data into RQFilrIO object_____________________________________________  
  RQFileIO *rqio = new RQFileIO;
  rqio->ReadFile(rq_dir+"/"+rq_filename);
  
    // Add the RQs that this module is going to generate_______________________
  size_t pulse_dim = dp_settings.get<int>("data_processing_settings.global.max_num_pulses");
  string pulse_dim_str=dp_settings.get<string>("data_processing_settings.global.max_num_pulses");
  
  rqio->events.AddRQ("cor_x_cm", "float", pulse_dim_str, -99);
  rqio->events.AddRQ("cor_y_cm", "float", pulse_dim_str, -99);
  
  // Grab the (new and pre-existing) RQs needed for algorithm__________________ 
  Rq* peak_area = rqio->events.Get("peak_area_phe");
  Rq* xcm = rqio->events.Get("cor_x_cm");
  Rq* ycm = rqio->events.Get("cor_y_cm");
  vector<double> areas(122);
  double x,y;
  //-----EVENT PROCESSING LOOP BEGINS-----
  for (unsigned int e=0; e<rqio->events.GetNSequences(); e++) {
    rqio->GetEvent(e);
    
	for (size_t i=0; i<pulse_dim; i++){
		for (size_t p=0; p<122; p++){

			areas[p] = peak_area->GetDouble(i, p);
		}
		CorrectedCentroidPosRecon(areas, x, y);
			xcm->SetDouble(x, i);
			ycm->SetDouble(y, i);
    }
  }

 //-----EVENT PROCESSINiG LOOP ENDS-----
  
  // Write the output rqs
  rqio->WriteFile(rq_dir+"/"+rq_filename);

  //clean up___________________________________________________________________
  delete rqio;

  return 0;
}
