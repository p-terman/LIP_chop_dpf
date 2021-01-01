#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "TTreeCacheUnzip.h"
#include "TH1I.h"
#include "TCanvas.h"

#include "Lux_EVT_Header.h"
#include "Lux_EVT_Channel.h"
#include "Lux_EVT_Event.h"
#include "Lux_EVT_Pulse.h"

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <ctime>
#include <vector>

// photon energy per mV convertion
#define Gain 2e6*5*1.5
#define R 25

using namespace std;

// a list of the gain each channel has
double channelgain[4]= {4.04, 4.13, 4.07, 4.09};
////////////////////////////////////////////////////////////////////////////////////////////////////
struct Peak{  // stored peak information
	int nchan;		// channel number
	long long truestart;	// start of peak (samples)
	long long trueend;		// end of peak (samples)
	long long start;   // start after compensating for noise
	long long end;
	long long length;
	double area;
	double height;
	double baseline;
	bool sat;
	double *data;				// histogram data
	
	// positions
	int p_height;
	int p_start;
	int p_end;
};
struct SimplePulse{ // stucture created to prevent some allocation error with vector pointers
	long long truestart;
	long long trueend;
	long long start; // start time
	long long end; // end time
	int length;
	vector<Peak> peak;
};
////////////////////////////////////////////////////////////////////////////////////////////////////
class rq1Pulse{
	// class that holds a reduced quantity pulse information
	private:
		long long start; 						// start of the pulse
		long long end;
		long long length;
		int p_start; 								// start position of the pulse
		int p_end; 									// end position of the pulse
		int p_height; 							// pulse max height position
		
		double *pdata; 							// will hold the bin data of the pulse
		double *ptime; 							// will hold the times
		
		void findChannelInfo(Peak peak, long long *startpos, long long *endpos, int index);
		void findPercents();
		void findArea();
		
	public:
		rq1Pulse(int nchan);  // empty constructor: where all value are set to 0
		rq1Pulse(SimplePulse pulse, int nchan ,double timeRes); 
		
		int npeaks; 		// number of peaks
		double t0;			// start of the pulse // note that "long long" is the same as "Long_t"
		double t10l;		// 10% rise time of the pulse
		double t50l;		// 50% rise time of the pulse
		double t1;			// time to pulse max
		double t50r;		// 50% fall time of the pulse
		double t10r;		// 10% fall time of the pulse
		double t2;			// end of the pulse
		double prompt_fraction; 			// area of the pulse between t0 and t1, divided by the total area
		double pulse_height;			// max height of the pulse
		double pulse_area; 
		
		// peak information; 
		int *channel;							// which channels							
		double *peak_area;				// peak area
		double *peak_height;			// maximum peak height
		double *baseline;					// baseline value as recorded by struck(mV)		
};
////////////////////////////////////////////////////////////////////////////////////////////////////
Peak PeakFill(Lux_EVT_Pulse *sig);
vector<SimplePulse> CombinedSignals(Lux_EVT_Event* event);
////////////////////////////////////////////////////////////////////////////////////////////////////
void writeRQ1File(int argc, char* argv[]) // real main function
{
	// get the file
	char *inFileName = argv[1];
	TFile f(inFileName);
	TTree *T = (TTree*) f.Get("Lux_EVT_Tree");
	cout << "Loaded file " << inFileName;
	
	// allocate memory for Lux_EVT class objects
	Lux_EVT_Event *event = new Lux_EVT_Event();
	Lux_EVT_Header *header = new Lux_EVT_Header();
	// set branch addresses
	T->SetBranchAddress("event", &event);
	T->SetBranchAddress("header", &header);
	T->GetEntry(0);
	int nevents = (int)T->GetEntries();
	int nchannels = event->Get_numChannels();
	double timeRes = event->Get_channel(0)->Get_timeRes();
	cout << ", loaded tree, and assigned branches\n";
	
	// data holding variables
	vector<SimplePulse> pulses;
	int *channelNumber = new int[nchannels];
	for (int i=0; i<nchannels; i++)	channelNumber[i] = event->Get_channel(0)->Get_channelNumber();
	bool **sat = new bool*[nevents];
	for (int i=0; i<nevents; i++){
		sat[i] = new bool[nchannels];
		for (int j=0; j<nchannels; j++)
			sat[i][j] = false;
	}
	vector<rq1Pulse> rq1pulses[nevents];
	long long *event_number = new long long[nevents];
	unsigned long long *event_timestamp = new unsigned long long[nevents];
	int max_pulses = 0;
	
	cout << nevents << " events in file\n";
	for (int ev=0; ev<nevents; ev++) // loop through events
	{
		cout << "\nEvent number " << event->Get_eventNum() << ": ";
		T->GetEntry(ev);
		pulses = CombinedSignals(event); // find the pulses in the event
		
		for (unsigned int p=0; p<pulses.size(); p++){
			// get the pulse characteristics of the pulse
			rq1Pulse temp(pulses[p], nchannels, timeRes);
			rq1pulses[ev].push_back( temp );
			
			// load the saturation information
			for (int ch=0; ch<nchannels; ch++){ 
				for (unsigned int i=0; i<pulses[p].peak.size(); i++)
					if (channelNumber[ch] == pulses[p].peak[i].nchan and !sat[ev][ch])
						sat[ev][ch] = pulses[p].peak[i].sat;
			}
		}
		if ( rq1pulses[ev].size() > max_pulses )
			max_pulses = rq1pulses[ev].size();
		
		event_number[ev] = event->Get_eventNum();
		event_timestamp[ev] = event->Get_timeStamp();
	}
	cout << "\nMaximum # of pulses: " << max_pulses << endl;
	cout << "Number of channels: " << nchannels << endl;
	cout << "Writing File\n";
	
	// WRITE TO FILE
	// create file output name
	string outFileName = inFileName;
	size_t position = outFileName.rfind(".root");
	outFileName.erase(position, 5);
	outFileName += ".rq1";
	
	ofstream file(outFileName.c_str(), ios::out | ios::binary);
	int endian = 0x01020304;
	file.write(reinterpret_cast<char*>(&endian),sizeof(int));
	
	// CREATE the header block
	// declare varibles needed for header
	int headLen;
	char headString[2000];
	int nlines;
	
	// find Dataset_Name 
	position = outFileName.rfind("_f");
	string dataset = outFileName;
	dataset.erase(position);
	
	// create header's header string
	sprintf(headString,"Dataset_Name;char;%i;Number_of_Events;ULong_t;1;First_Event;ULong_t;1;",
					(int)dataset.size());
	
	// write header block
	headLen = strlen(headString); nlines = 1;
	file.write(reinterpret_cast<char*>(&headLen),sizeof(int)); // length of header string
	file.write(headString, sizeof(char)); // the header string
	file.write(reinterpret_cast<char*>(&nlines),sizeof(int)); // number of lines	
	
	// data
	file.write(reinterpret_cast<char*>(&dataset),sizeof(char));
	ULong_t temp = header->Get_numEvents();
	file.write(reinterpret_cast<char*>(&temp),sizeof(ULong_t));
	temp = header->Get_firstEvent();
	file.write(reinterpret_cast<char*>(&temp),sizeof(ULong_t));
	cout << "Wrote header block\n";
	
	// CREATE the rq1 block
	sprintf(headString, "t0;double;%i;t10l;double;%i;t50l;double;%i;t1;double;%i;\
	t50r;double;%i;t10r;double;%i;t2;double;%i;prompt_fraction;double;%i;num_peaks;int;%i;\
	chan_number;int;%i,%i;peak_area;double;%i,%i;peak_height;double;%i,%i;\
	baseline;double;%i,%i;sat;bool;%i;evt_num;ULong_t;1;evt_timestamp;ULong64_t;1;"
	,max_pulses, max_pulses, max_pulses, max_pulses, max_pulses, max_pulses, max_pulses, max_pulses,
	 max_pulses, max_pulses, nchannels, max_pulses, nchannels, max_pulses, nchannels, 
	 max_pulses, nchannels, nchannels);
	
	// write rq1 block
	headLen = strlen(headString); nlines = nevents;
	file.write(reinterpret_cast<char*>(&headLen),sizeof(int)); // length of header string
	file.write(headString, sizeof(char)); // the header string
	file.write(reinterpret_cast<char*>(&nlines),sizeof(int)); // number of lines	
	cout << "Starting RQ1 block\t";
	for (int ev=0; ev<nevents; ev++){   // output data
		// fill all the events with the same number of pulse by adding empty pulses
		for (unsigned int i=rq1pulses[ev].size(); i<max_pulses; i++)
			rq1pulses[ev].push_back( rq1Pulse(nchannels) );
		
		for (int i=0; i<max_pulses; i++) // t0
			file.write(reinterpret_cast<char*>(&rq1pulses[ev][i].t0),sizeof(double));
		
		for (int i=0; i<max_pulses; i++) // t10l
			file.write(reinterpret_cast<char*>(&rq1pulses[ev][i].t10l),sizeof(double));
		
		for (int i=0; i<max_pulses; i++) // t50l
			file.write(reinterpret_cast<char*>(&rq1pulses[ev][i].t50l),sizeof(double));
		
		for (int i=0; i<max_pulses; i++) // t1
			file.write(reinterpret_cast<char*>(&rq1pulses[ev][i].t1),sizeof(double));
		
		for (int i=0; i<max_pulses; i++) // t50r
			file.write(reinterpret_cast<char*>(&rq1pulses[ev][i].t50r),sizeof(double));
		
		for (int i=0; i<max_pulses; i++) // t10r
			file.write(reinterpret_cast<char*>(&rq1pulses[ev][i].t10r),sizeof(double));
		
		for (int i=0; i<max_pulses; i++) // t2
			file.write(reinterpret_cast<char*>(&rq1pulses[ev][i].t2),sizeof(double));
		
		for (int i=0; i<max_pulses; i++) // prompt_fraction
			file.write(reinterpret_cast<char*>(&rq1pulses[ev][i].prompt_fraction),sizeof(double));
		
		for (int i=0; i<max_pulses; i++) // number of peaks
			file.write(reinterpret_cast<char*>(&rq1pulses[ev][i].npeaks),sizeof(int));
		
		for (unsigned int i=0; i<max_pulses; i++) // channel number
			for (int j=0; j<nchannels; j++)
				file.write(reinterpret_cast<char*>(&rq1pulses[ev][i].channel[j]),sizeof(int));
		
		for (int i=0; i<max_pulses; i++) // peak_area
			for (int j=0; j<nchannels; j++)
				file.write(reinterpret_cast<char*>(&rq1pulses[ev][i].peak_area[j]),sizeof(double));
		
		for (int i=0; i<max_pulses; i++) // peak_height
			for (int j=0; j<nchannels; j++)
				file.write(reinterpret_cast<char*>(&rq1pulses[ev][i].peak_height[j]),sizeof(double));
		
		for (int i=0; i<max_pulses; i++) // baseline
			for (int j=0; j<nchannels; j++)
				file.write(reinterpret_cast<char*>(&rq1pulses[ev][i].baseline[j]),sizeof(double));
		
		for (int i=0; i<nchannels; i++) // saturation
			file.write(reinterpret_cast<char*>(&sat[ev][i]),sizeof(bool));
		
		file.write(reinterpret_cast<char*>(&event_number[ev]),sizeof(long long)); // event number
		file.write(reinterpret_cast<char*>(&event_timestamp[ev]),sizeof(unsigned long long)); // timestamp
	}
	cout << "finished RQ1 block\n";
	
	// CREATE Livetime block
	sprintf(headString,"timeStampLatch;ULong_t;1;timeStampEnd;ULong_t;1");
	headLen = strlen(headString); nlines = header->Get_numSeq();
	file.write(reinterpret_cast<char*>(&headLen),sizeof(int)); // length of header string
	file.write(headString, sizeof(char)); // the header string
	file.write(reinterpret_cast<char*>(&nlines),sizeof(int)); // number of lines
	
	for (int i=0; i<nlines; i++){
		temp = header->Get_timeStampLatch(i);
		file.write(reinterpret_cast<char*>(&temp),sizeof(ULong_t));
		temp = header->Get_timeStampEnd(i);
		file.write(reinterpret_cast<char*>(&temp),sizeof(ULong_t));
	}
	cout << "wrote Livetime block\n";
	// clean up
	file.close();
	delete event;
	delete header;
	delete[] channelNumber;
	delete[] event_number;
	delete[] event_timestamp;
	for (int i=0; i<nevents; i++) delete[] sat[i];
	delete[] sat;
}
////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
	writeRQ1File(argc, argv);
	return 0;
}
////////////////////////////////////////////////////////////////////////////////////////////////////
rq1Pulse::rq1Pulse(int nchan)
{ 
	// empty constructor used for cases if no pulse if found in an event
	t0=0;
	t10l=0;
	t50l=0;
	t1=0;
	t50r=0;
	t10r=0;
	t2=0;
	prompt_fraction = 0;
	npeaks = 0;
	channel = new int[nchan];
	peak_area = new double[nchan];
	peak_height = new double[nchan];
	baseline = new double[nchan];
	for (int i=0; i<nchan; i++){
		channel[i] = -1;
		peak_area[i] = 0;
		peak_height[i] = 0;
		baseline[i] = 0;
	}
}
//**************************************************************************************************
rq1Pulse::rq1Pulse(SimplePulse pulse, int nchan, double timeRes)
{
	// This function uses the raw pulse data to find the characteristics 
	// of the pulse: t0, t10l, t50l, t1, etc...
	
	// Transfer some the basic information for ease
	start = pulse.truestart;
	length = pulse.length;
	npeaks = pulse.peak.size();
	
	// allocate memory for each channel dependant variables and fill
	channel = new int[nchan];
	peak_area = new double[nchan];
	peak_height = new double[nchan];
	baseline = new double[nchan];
	for (int i=0; i<nchan; i++){
		if (i<npeaks){
			channel[i] = pulse.peak[i].nchan;
			peak_area[i] = pulse.peak[i].area;
			peak_height[i] = pulse.peak[i].height;
			baseline[i] = pulse.peak[i].baseline;
		}
		else{
			channel[i] = -1;
			peak_area[i] = 0;
			peak_height[i] = 0;
			baseline[i] = 0;
		}
	}
	// create temporary array to combined peak data
	pdata = new double[length];
	ptime = new double[length];
	//initilize start and end pointers
	p_start = 0; p_end = length;
	for (long long i=0; i<length; i++){
		pdata[i]=0; //initilize to 0
		ptime[i]=timeRes*(start+i);//+4970;
		
		// add up contributions from the peaks
		for (int j=0; j<npeaks; j++){
			long long k = start+i-pulse.peak[j].truestart;
			if ( k >= 0 and k < pulse.peak[j].length )
				pdata[i] += pulse.peak[j].data[k];
		}
		
		if (pulse.start == pulse.truestart + i) // finds t0
			p_start = i;
		if (pulse.end == pulse.truestart + i+1) // finds t2
			p_end = i;
		if (pdata[i] > pulse_height){ // finds pulse_height and t1
			pulse_height = pdata[i];
			p_height = i;
			t1 = ptime[i];
		}
		if (i==1){
			pulse_height = pdata[i];
			p_height = i;
		}
	}
	do{ // fix starting position
		if (pdata[p_start]>= pdata[p_start+1])
			p_start++;
		else
			break;
	}while(p_start<p_height-1);
	do{ // fix ending position
		if (pdata[p_end]>= pdata[p_end-1])
			p_end--;
		else
			break;
	}while(p_end>p_height+1);
	// assign t0, t2
	t0 = ptime[p_start];
	t2 = ptime[p_end];
	
	findPercents();
	findArea();
	
	// delete allocated memory
	delete[] pdata;
	delete[] ptime;
}
//**************************************************************************************************
void rq1Pulse::findPercents()
{
	// FIND t10l, t10r, t50l, t50r
	// declare variable used for interpolation
	double x1[2], x3[2], y1[2]={ptime[0],ptime[0]}, y3[2]={0,0};
	
	// find the left side first
	int t10pos; // used for a check near the end
	for (int i=p_height;i>p_start && i>0;i--){	//int i=p_start; i<p_height; i++){
		if ( pdata[i]>=.10*pulse_height and pdata[i-1]<=.10*pulse_height ){
			x3[0] = pdata[i]; x1[0] = pdata[i-1]; y3[0] = ptime[i]; y1[0] = ptime[i-1];
			t10pos = i-1;
		}
		if ( pdata[i]>=.50*pulse_height and pdata[i-1]<=.50*pulse_height ){
			x3[1] = pdata[i]; x1[1] = pdata[i-1]; y3[1] = ptime[i]; y1[1] = ptime[i-1];
			//break;
		}
	}
	if (!y3[0]){ // didn't find t10l: look in the previous 4 bins
		for (int i=p_start; i>p_start-4 && i>0; i--)	//int i=p_start-4; i<p_start; i++)
			if ( pdata[i]>=.10*pulse_height and pdata[i-1]<=.10*pulse_height ){
				x3[0] = pdata[i]; x1[0] = pdata[i-1]; y3[0] = ptime[i]; y1[0] = ptime[i-1];
				t10pos = i-1;
				break;
			}
	}
	if (!y3[1]){ // didn't find t50l: look in the previous 4 bins
		for (int i=p_start; i>p_start-4 && i>0; i--)	//int i=p_start-4; i<p_start; i++)
			if ( pdata[i]>=.50*pulse_height and pdata[i-1]<=.50*pulse_height ){
				x3[1] = pdata[i]; x1[1] = pdata[i-1]; y3[1] = ptime[i]; y1[1] = ptime[i-1];
				break;
			}
	}
	if (!y3[0]){ // still didn't find t10l: look from the beginning
		for (int i=1; i<p_start-4; i++)
			if ( pdata[i]>=.10*pulse_height and pdata[i-1]<=.10*pulse_height ){
				x3[0] = pdata[i]; x1[0] = pdata[i-1]; y3[0] = ptime[i]; y1[0] = ptime[i-1];
				t10pos = i-1;
				break;
			}
	}
	if (!y3[1]){ // still didn't find t50l: look from the beginning
		for (int i=1; i<p_start-4; i++)
			if ( pdata[i]>=.50*pulse_height and pdata[i-1]<=.50*pulse_height ){
				x3[1] = pdata[i]; x1[1] = pdata[i-1]; y3[1] = ptime[i]; y1[1] = ptime[i-1];
				break;
			}
	}//cout << y3[0] << " " << y3[1] << "\t";
	// interpolate to get t10l and t50l
	t10l = (.10*pulse_height-x1[0])*(y3[0]-y1[0])/(x3[0]-x1[0]) + y1[0];
	t50l = (.50*pulse_height-x1[1])*(y3[1]-y1[1])/(x3[1]-x1[1]) + y1[1];
	// fix start position if t10l come earlier
	if (t10l < t0)	{ t0 = ptime[t10pos]; p_start = t10pos;}
	
	// find the right side: similiar algorithm as earlier
	y3[0] = 0; y3[1] = 0; y1[0] = ptime[0]; y1[1] = ptime[0]; //reset
	for (int i=p_height; i<p_end; i++){
		if ( pdata[i]<=.10*pulse_height and pdata[i-1]>=.10*pulse_height ){
			x3[0] = pdata[i]; x1[0] = pdata[i-1]; y3[0] = ptime[i]; y1[0] = ptime[i-1];
			t10pos = i;
			break;
		}
		if ( pdata[i]<=.50*pulse_height and pdata[i-1]>=.50*pulse_height ){
			x3[1] = pdata[i]; x1[1] = pdata[i-1]; y3[1] = ptime[i]; y1[1] = ptime[i-1];
		}
	}
	if (!y3[0]){ // didn't find t10r: look in the next 4 bins
		for (int i=p_end; i<p_end+4 && i<length; i++)
			if ( pdata[i]<=.10*pulse_height and pdata[i-1]>=.10*pulse_height ){
				x3[0] = pdata[i]; x1[0] = pdata[i-1]; y3[0] = ptime[i]; y1[0] = ptime[i-1];
				t10pos = i;
				break;
			}
	}
	if (!y3[1]){ // didn't find t50r: look in the next 4 bins
		for (int i=p_end; i<p_end+4 && i<length; i++)
			if ( pdata[i]<=.50*pulse_height and pdata[i-1]>=.50*pulse_height ){
				x3[1] = pdata[i]; x1[1] = pdata[i-1]; y3[1] = ptime[i]; y1[1] = ptime[i-1];
				break;
			}
	}
	if (!y3[0]){ // still didn't find t10r: look all the way to the end
		for (int i=p_end+4; i<length; i++)
			if ( pdata[i]<=.10*pulse_height and pdata[i-1]>=.10*pulse_height ){
				x3[0] = pdata[i]; x1[0] = pdata[i-1]; y3[0] = ptime[i]; y1[0] = ptime[i-1];
				t10pos = i;
				break;
			}
	};
	if (!y3[1]){ // still didn't find t50r: look all the way to the end
		for (int i=p_end+4; i<length; i++)
			if ( pdata[i]<=.50*pulse_height and pdata[i-1]>=.50*pulse_height ){
				x3[1] = pdata[i]; x1[1] = pdata[i-1]; y3[1] = ptime[i]; y1[1] = ptime[i-1];
				break;
			}
	}//cout << y3[0] << " " << y3[1] << "\t";
	// interpolate to get t10l and t50l
	t10r = (.10*pulse_height-x1[0])*(y3[0]-y1[0])/(x3[0]-x1[0]) + y1[0];
	t50r = (.50*pulse_height-x1[1])*(y3[1]-y1[1])/(x3[1]-x1[1]) + y1[1];
	// fix end position if t10r is larger
	if (t10r > t2)	{ t2 = ptime[t10pos]; p_end = t10pos;}
}
//**************************************************************************************************
void rq1Pulse::findArea()
{
	// FIND the pulse area, prompt fraction
	// initilize
	prompt_fraction = 0;
	pulse_area = 0;
	for (int i=p_start; i<p_end && i<length; i++)
	{
		pulse_area += .5*(pdata[i]+pdata[i+1]);
		if (i < p_height )
			prompt_fraction += .5*(pdata[i]+pdata[i+1]); // using trapezoidal rule
	}
	prompt_fraction /= pulse_area;
}
////////////////////////////////////////////////////////////////////////////////////////////////////
vector<SimplePulse> CombinedSignals(Lux_EVT_Event* event)
{
	// this function finds the peaks associated with a certain pulse by matching peaks comparing the
	// individual peaks start and end times and groups peaks together who's start and end times
	// overlap.  Also loads some basic data.
	
	vector<SimplePulse> pulse;  // hold my information of which signals come from which pulse
	Lux_EVT_Pulse *sig;
	int numberofpulses = 0;
	bool New;
	cout << " " << event->Get_numPulses() << " number of peaks; ";
	for (ULong_t i=0; i<event->Get_numPulses(); i++){ // loop through all the signals
		sig = event->Get_pulse(i); 
		New = true; // assume we have a new pulse
		
		Peak peak = PeakFill(sig);
		
		for (int j=0; j<numberofpulses; j++){ // go through all other pulse information to find matches
			bool SameChannel = false;
			for (int k=0; k<pulse[j].peak.size(); k++)
				if (peak.nchan == pulse[j].peak[k].nchan)
					SameChannel = true;			
			
			if ( pulse[j].start<=peak.start and pulse[j].end>peak.start
					 and !SameChannel){
			// if peak starts within the range 
				New = false;
				pulse[j].peak.push_back(peak); // add peak to pulse
				// check pulse boundaries
				if (peak.start < pulse[j].start) // finds the minimum start
					pulse[j].start = peak.start;
				if (peak.end > pulse[j].end) // finds the maximum end
					pulse[j].end = peak.end;
				if (peak.trueend > pulse[j].trueend)
					pulse[j].trueend = peak.trueend;
				if (peak.truestart < pulse[j].truestart)
					pulse[j].truestart = peak.truestart;
			}
			else if ( pulse[j].start<peak.end and pulse[j].end>=peak.end
								and !SameChannel){
			// if peak ends within the range
				New = false;
				pulse[j].peak.push_back(peak);
				// check pulse boundaries
				if (peak.start < pulse[j].start) // finds the minimum start
					pulse[j].start = peak.start;
				if (peak.end > pulse[j].end) // finds the maximum end
					pulse[j].end = peak.end;
				if (peak.trueend > pulse[j].trueend)
					pulse[j].trueend = peak.trueend;
				if (peak.truestart < pulse[j].truestart)
					pulse[j].truestart = peak.truestart;
			}
			else if ( pulse[j].start>=peak.start and pulse[j].end<=peak.end
								and !SameChannel){
			// if the range is within the peak's range
				New = false;
				pulse[j].peak.push_back(peak);
				pulse[j].start = peak.start;
				pulse[j].end = peak.end;
				// check pulse boundaries
				if (peak.start < pulse[j].start) // finds the minimum start
					pulse[j].start = peak.start;
				if (peak.end > pulse[j].end) // finds the maximum end
					pulse[j].end = peak.end;
				if (peak.trueend > pulse[j].trueend)
					pulse[j].trueend = peak.trueend;
				if (peak.truestart < pulse[j].truestart)
					pulse[j].truestart = peak.truestart;
			}
			pulse[j].length = pulse[j].trueend - pulse[j].truestart; // fix the length
		} // ends peak matches loop (j)
		
		if (New){ // have a new pulse
			pulse.resize(++numberofpulses);
			pulse[numberofpulses-1].start = peak.start;
			pulse[numberofpulses-1].end = peak.end;
			pulse[numberofpulses-1].truestart = peak.truestart;
			pulse[numberofpulses-1].trueend = peak.trueend;
			pulse[numberofpulses-1].peak.push_back(peak);
		}
	} // ends pulses loop (i)
	cout << numberofpulses << " pulse(s) found:\t";
	return pulse;
}
////////////////////////////////////////////////////////////////////////////////////////////////////
Peak PeakFill(Lux_EVT_Pulse *sig)
{	
	Peak peak;
	// transfer peak data
	peak.nchan = sig->Get_channelNumber();
	peak.baseline = sig->Get_pulseBaseline();
	peak.truestart = sig->Get_pulseStarts();
	peak.length = sig->Get_pulseLength();
	peak.trueend = peak.truestart + peak.length;
	peak.data = new double[peak.length];
	// inverse peak.dataset
	for (int i=0; i<peak.length; i++)
		peak.data[i] = -2000./pow(2,14)*(sig->Get_pulseData()->GetBinContent(i+1))*9.538/( 1.6e-7*Gain*R);
	peak.data[0] = 0; // some artifact that keeps getting in the way
	
	// declaration of variables
	int tailmean = 0 , tailsigma = 0; 
	int tailpoints = 14; // how many points to use from the tailmean
	double threshold; // distinguish between noise and signal
	
	// FIND the peak height and the mean value of the tail
	peak.height = peak.data[1];
	peak.p_height = 1;
	for (int i=2; i<peak.length; i++){
		if ( i >= peak.length-tailpoints ) // average over the tail
			tailmean += peak.data[i];
		if ( peak.height < peak.data[i]){ // find the peak height and its location
			peak.height = peak.data[i];
			peak.p_height = i;
		}
	}
	tailmean /= tailpoints;
	
	// find standard deviation of the tail
	for (int i=peak.length-tailpoints; i<peak.length; i++)
		tailsigma += pow(peak.data[i]-tailmean, 2);
	tailsigma = sqrt(tailsigma/(tailpoints-1));

	// FIND the start and end of a peak
	// set threshold value dependant of peak height
	if (peak.height >= 10) 	threshold = 2;
	else	threshold = .2;
	// initilize start and end
	peak.start = peak.truestart; peak.p_start=3;
	peak.end = peak.trueend-tailpoints; peak.p_end = peak.length-tailpoints;
	
	for (int i=4; i<peak.p_height+1; i++) // find rough position of start position		
		if ( peak.data[i+1]>=threshold
			 	 and peak.data[i]>=threshold
			 	 and peak.data[i-1]>=threshold ){
			peak.start = i+peak.truestart-1; peak.p_start = i-1;
			break;
		}
	
	for (int i=peak.p_start-3; i<peak.p_height+1; i++) // refine the start position
		if (peak.data[i]>tailmean+5*tailsigma  and peak.data[i]>0.00001)
		{	peak.start = i+peak.truestart-1; peak.p_start = i-1; break;}
	
	// set starting position such that it ignores the tail end
	int j = peak.length - tailpoints;
	if (peak.height < 5)	j -= peak.p_start;
	
	for (; j>peak.p_height-1; j--){ // find rough position of end position
		if ( peak.data[j+1]>=threshold 
				 and peak.data[j]>=threshold 
				 and peak.data[j-1]>=threshold ){
			peak.end = j+peak.truestart+1; peak.p_end = j+1;
			break;
		}
	}
	
	for (j=peak.p_end+4; j>peak.p_height-1; j--) // refine the end position
		if (peak.data[j]>tailmean+5*tailsigma and peak.data[j]>0.00001)
		{	peak.end = j+peak.truestart+1; peak.p_end = j+1; break;	}
		
	// find peak area and saturation
	peak.area = 0.;
	peak.sat = false;
	for (int i=peak.p_start; i<peak.p_end; i++){
		peak.area += .5*(peak.data[i]+peak.data[i+1]);
		if ( peak.data[i] == peak.data[i+1] and peak.data[i] == peak.data[i+2]
				 and peak.data[i] == peak.height) // checks for saturation
			peak.sat = true;
	}
	//peak.area *= 9.538/( 1.6e-7*25*2e6 );
	
	return peak;
}
