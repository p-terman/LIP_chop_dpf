//
//  main.cpp
//  PositionReconstruction
//
//  Created by Chang Lee on 3/20/13.
//  Copyright (c) 2013 Chang Lee. All rights reserved.
//

#include <iostream>
#include <numeric>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>

#include "Pulse.h"
#include "PBomb.h"

using namespace std;
using boost::property_tree::ptree;

//__________________________________________________________________________
void get_subtree(ptree& full_tree, ptree& sub_tree, string type, string idkey, string idvalue){ //e.g. type is iq, idkey is global.iq_type, idvalue is LRF
	BOOST_FOREACH(ptree::value_type &v, full_tree.get_child(type)) {
		if(v.second.count(idkey)==1&&v.second.get<string>(idkey)==idvalue){
			sub_tree = v.second;
		}
	}
}
//__________________________________________________________________________
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
//__________________________________________________________________________
void load_sigma_per_phe(ptree& pt, vector<double> &v){
	v.clear();
	v.resize(122);
	for(ptree::iterator iter=pt.begin();iter!=pt.end();iter++){
		if(iter->second.count("iq_type")==0||iter->second.get<string>("iq_type")!="pmt_gains") continue;
		if(iter->second.count("fit")==0) continue;
 
		BOOST_FOREACH(ptree::value_type &w,iter->second.get_child("fit")){
			if(w.first=="channel"){
				int ind = w.second.get<int>("PMT_number")-1;
				v[ind]=w.second.get<double>("sigma_mVns_per_phe");
			}
		}
	}
 }
//__________________________________________________________________________

int main(int argc, const char * argv[])
{
    #if VERBOSITY>0
    cout << "Starting" << endl;
    cout << "I have " << argc << " arguments" << endl;
    for(int i=0; i<argc; i++) 
        cout << "\t" << i << ") " << argv[i] << endl;
    #endif

    //parse input arguments
    if( argc !=8  ) {
        cerr << "Please specify:\n-evt_filename, evt_dir, rq_filename, "
		<<"rq_dir, module no., dp_xml_full_path, iq_xml_full_path" << endl;
        return -1;
    }
	clock_t t = clock(); //tic

    string module_path = string(argv[0]);
    string evt_dir = string(argv[2]);
    string evt_filename = string(argv[1]);
    string rq_dir = string(argv[4]);
    string rq_filename =string(argv[3]);
    string dp_xml_full_path = string(argv[6]);
    string iq_xml_full_path = string(argv[7]);
    
	#if VERBOSITY>0
    cout<<evt_dir+"/"+evt_filename<<"\n"<<rq_dir+"/"+rq_filename<<"\n"
	<<dp_xml_full_path<<"\n"<<iq_xml_full_path<<endl;
	#endif
    
    // Read input xml files into BOOST ptrees
    ptree dp_settings, iqs;
   
    ifstream dp_settings_ifs;
	dp_settings_ifs.open(dp_xml_full_path.c_str());
    read_xml(dp_settings_ifs,dp_settings);
    dp_settings_ifs.close();
    
    ifstream iq_settings_ifs;
	iq_settings_ifs.open(iq_xml_full_path.c_str());
    read_xml(iq_settings_ifs,iqs);
    iq_settings_ifs.close();

	// Fill variables from xml 
	ptree module_settings;
	get_subtree(dp_settings,module_settings,"data_processing_settings",
				"module_name","PositionReconstruction_Hitmap");
	#if VERBOSITY>0
		cout<<module_settings.get<string>("relative_path")<<endl;
	#endif

	size_t keyType = module_settings.get<size_t>("parameters.keyType");
	size_t statMethod = module_settings.get<size_t>("parameters.chiSq0LH1");
	string hit_map_path_relative = module_settings.get<string>("parameters.hitmapPath");
	string disabledPMTStr = module_settings.get<string>("parameters.disabledPMTs");
	vector<size_t> disabledPMTs = Dimensions(disabledPMTStr);
	size_t nPBombInterpolate = module_settings.get<size_t>("parameters.nPBombInterpolate");
	size_t nKey = module_settings.get<size_t>("parameters.nKey");
	size_t pulse_dim = dp_settings.get<size_t>("data_processing_settings.global.max_num_pulses");
	string pulse_dim_str=dp_settings.get<string>("data_processing_settings.global.max_num_pulses");
	
	// handle exceptions in dp module settings
	if ((keyType!=0) & (keyType!=1)){
		cerr<<"Choose key generation method: 0 for scarce map, 1 for hottest PMT";
		return -1;
	}
	if ((statMethod!=0) & (statMethod!=1)){
		cerr<<"Choose statistics method: 0 for chi-squared, 1 for likelihood";
		return -1;
	}
	if (nKey>=nTopPMT){
		cerr<<"Choose nKey smaller than 61";
		return -1;
	}
	
	//Load data into RQFileIO object
	RQFileIO *rqio = new RQFileIO;
	if(rqio->ReadFile(rq_dir+"/"+rq_filename)==false)
	{cerr<<"Failed to read RQ file"<<endl; return -1;}
	
	// Add the RQs that this module is going to generate
	rqio->events.AddRQ("x_cm_tmplt", "float", pulse_dim_str, -99);
	rqio->events.AddRQ("y_cm_tmplt", "float", pulse_dim_str, -99);
	rqio->events.AddRQ("xy_sigma_cm", "float", pulse_dim_str);
	rqio->events.AddRQ("xy_chisq_p", "float", pulse_dim_str);

    // Load Rq pointers
	Rq* event_number = rqio->events.Get("event_number");
	Rq* areas = rqio->events.Get("peak_area_phe");
	Rq* pulse_start_samples = rqio->events.Get("pulse_start_samples");
	Rq* pulse_end_samples = rqio->events.Get("pulse_end_samples");
	Rq* tbasym = rqio->events.Get("top_bottom_asymmetry");	
	Rq* pulseClass = rqio->events.Get("pulse_classification");	
	Rq* x_cm_tmplt = rqio->events.Get("x_cm_tmplt");
	Rq* y_cm_tmplt = rqio->events.Get("y_cm_tmplt");
	Rq* xy_sigma_cm = rqio->events.Get("xy_sigma_cm");
	Rq* xy_chisq_p = rqio->events.Get("xy_chisq_p");

	//	Read HitMap
	HitLib * Hitmap = new HitLib;
	string hit_map_path = module_path.substr(0,module_path.rfind('/'))+'/'+hit_map_path_relative;
	if(Hitmap->BuildFrom(hit_map_path, disabledPMTs, keyType)!=0)
	{cerr<<"Failed to build a Hit Map"<<endl; return -1;}
	
	//	Load gains and erros, then calculate the gain uncertainty
	vector<double> mVns_per_phe, sigma_per_phe;
	vector<float> ChiUnc;
	load_mVns_per_phe(iqs,mVns_per_phe);
	load_sigma_per_phe(iqs,sigma_per_phe);
	for (size_t un=0; un<mVns_per_phe.size(); un++)
		ChiUnc.push_back(1+pow(sigma_per_phe[un]/mVns_per_phe[un],2));

	//	Prepare variables for the loops
	long int pulse_width;
	size_t nPulse = 0;
	
	// Loop through the events
	for (size_t e=0; e<rqio->events.GetNSequences(); e++) {
		
		// Load an event
		rqio->GetEvent(e);
		#if VERBOSITY>0
		cout<<endl<<"Event: "<<event_number->GetInt();
		#endif
		
		// Loop through the pulses
		for (size_t p=0; p<pulse_dim; ++p) {
			
			// Create an obj with known info
			Pulse* pls = new Pulse(event_number->GetInt(), p, Hitmap->pmtIdx, Hitmap->liveTopPMT, statMethod);

			// read pulse duration, and suum top peak areas, skip all empty pulses
			pulse_width = pls->WidthCheck(pulse_end_samples,
										  pulse_start_samples);
			
			// if the pulse is not filled, jump to the next event
			if (pulse_width<0){
				pls->FillNaN(x_cm_tmplt,y_cm_tmplt);
				delete pls;
				break;	
			}
			
			// if empty, jump to the next pulse
			if (pulse_width==0 | 
				pls->ReadTopAreaAndCheck(areas,10)!=0 | 
				tbasym->GetDouble(p)==-1) {
				pls->FillNaN(x_cm_tmplt,y_cm_tmplt);
				delete pls;
				continue; 
			}

			// Once the pulse passes all exceptions, raise count
			++nPulse;
			
			// generate keys for narrow search
			if (keyType==0)	// scarce map key
				pls->GenerateScarceMapKey(nKey, Hitmap, ChiUnc);

			if (keyType==1) // hottest PMT key
				pls->GeneratePMTKey(nKey);
			
			pls->SelectPBombsFromKey(Hitmap);			
			
			// loop through the selected PBombs
			pls->PBombLoop(ChiUnc);
			
			// Interpolate the best PBombs
			pls->InterpolateXY(nPBombInterpolate);
			
			// record the final values
			x_cm_tmplt->SetDouble(pls->ReadFinalX(),p);
			y_cm_tmplt->SetDouble(pls->ReadFinalY(),p);
			if (0==statMethod)
				xy_chisq_p->SetDouble(pls->ReadChiSqP(),p);
			if (1==statMethod)
				xy_sigma_cm->SetDouble(pls->ReadSigmaRad(),p);
	
			// Clean the pulse object
			pls->Clear();
			delete pls;
			
		} // Pulse loop
	} // Event loop
	
    // Write the output rqs
	rqio->WriteFile(rq_dir+"/"+rq_filename);
	
	// Print out the results
	t=clock()-t; //toc
	cout<<endl<<"Position Reconstruction module processed "
		<<rqio->events.GetNSequences()<<" events and "<<nPulse
		<<" pulses in "<<((float)t)/CLOCKS_PER_SEC<<" seconds."<<endl;
	
	//clean up
	rqio->Clear();
	delete rqio;
	Hitmap->Delete();
	delete Hitmap;
    
    return 0;
}
