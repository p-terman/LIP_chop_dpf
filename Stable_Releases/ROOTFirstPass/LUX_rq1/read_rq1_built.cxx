// New intent is to store data in a root file which is updated on darksurfer
// for analysis purposes.  This should be automatic but isn't.
//
// K.Clark - 22/02/10
// ken.clark@case.edu

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <ctype.h>
#include <vector>
#include <string>

#include <TROOT.h>
#include <TObject.h>
#include <TTree.h>
#include <TFile.h>
#include <TGraph.h>
#include <TCanvas.h>

#include "LUX_rq1_event.h"
#include "LUX_rq1_pulse.h"
#include "LUX_rq1_header.h"
#include "LUX_rq1_channel.h"

#define full_endian_swap(x) endian_swap((char *) &x,sizeof(x))

using namespace std;

void endian_swap(char * b, int n){ // taken from http://www.codeproject.com/KB/cpp/endianness.aspx
  register int i = 0;
  register int j = n-1;
  while (i<j){
    std::swap(b[i], b[j]);
    i++, j--;
  }
}

Int_t stupid_method(char *type){

  if (!strcmp(type,"char") || !strcmp(type,"int8") || !strcmp(type,"uint8")) return 1;
  else if (!strcmp(type,"int16") || !strcmp(type,"uint16")) return 2;
  else if (!strcmp(type,"int32") || !strcmp(type,"uint32") || !strcmp(type,"float") || !strcmp(type,"single")) return 4;
  else if (!strcmp(type,"int64") || !strcmp(type,"uint64") || !strcmp(type,"double")) return 8;
  else return 16;

}

void read_rq1(char *file) {

  TCanvas c1;
  Int_t split=2;
  Int_t bsize=51200000;

  ifstream infile;
  printf("Trying to open file %s\n",file);
  infile.open(file,ios::in | ios::binary);
  if (infile.is_open()){
    strtok(file,".");
    strcat(file,".root");
    TFile *rootfile = new TFile(file,"recreate");
    TTree *LUX_rq1_tree = new TTree("LUX_rq1_tree","Tree of rq1 data");
    LUX_rq1_header *lrq1_header = new LUX_rq1_header();
    LUX_rq1_event *lrq1_event = new LUX_rq1_event();
    LUX_rq1_pulse *lrq1_pulse = new LUX_rq1_pulse();
    LUX_rq1_channel *lrq1_channel = new LUX_rq1_channel();
    vector<LUX_rq1_channel> channel_vector;
    LUX_rq1_tree->Branch("Header","LUX_rq1_header",&lrq1_header,bsize,1);
    LUX_rq1_tree->Branch("Event","LUX_rq1_event",&lrq1_event,bsize,split);
    
    UInt_t temp_headerlength;
    Int_t headerlength;
    Char_t *temp_entry;
    Int_t num_vars;
    Int_t semicolon_ct=0;
    Int_t num_lines;
    Char_t *header;
    Int_t channel_int_ct=0;
    Int_t channel_dbl_ct=0;
    Int_t pulse_int_ct=0;
    Int_t pulse_dbl_ct=0;
    Int_t event_int_ct=0;
    Int_t event_dbl_ct=0;
 
    UInt_t **temp_int;
    Double_t **temp_dbl;
    UInt_t ***temp_int_array;
    Double_t ***temp_dbl_array;
    Char_t **var_names;
    Char_t **var_types;
    Int_t *var_count;
    Int_t *var_indices;
    Int_t *var_type_size;
    Char_t **temp_char;
    string temp_string;
    Char_t *strchr_val;
    Int_t num_len;
    Char_t *sub_str;
    Int_t n_pulses=0;
    Int_t n_channels=0;
    Bool_t swap_needed=0;
    Bool_t int_array_def_flag=0;
    Bool_t dbl_array_def_flag=0;
    Int_t index;

    infile.seekg(0,ios::beg);
    UInt_t endianness;
    infile.read((char*)&endianness,sizeof(endianness));
    if (endianness != 0x01020304) swap_needed=1;

    //Now it's time to read out the blocks of data.  There are two.
    for (Int_t block=0;block<3;block++){
      if (block == 0){ //ascii header
      infile.read((char*)&temp_headerlength,sizeof(UInt_t));
      if (swap_needed) full_endian_swap(temp_headerlength);
      headerlength = (Int_t)temp_headerlength;
      }
      else { //file header
	UShort_t temp_file_headerlength;
	infile.read((char*)&temp_file_headerlength,sizeof(UShort_t));
	if (swap_needed) full_endian_swap(temp_file_headerlength);
	headerlength = (Int_t)temp_file_headerlength;
      }
      
      header = new Char_t[headerlength];
      infile.read((Char_t*)header,(headerlength)*sizeof(Char_t));
      semicolon_ct=0;
      for (Int_t i=0;i<headerlength;i++) if (header[i]==';'){semicolon_ct++;}
      num_vars = semicolon_ct/3;
      if(num_vars==0){continue;}
      
      var_names = new Char_t*[num_vars];
      var_types = new Char_t*[num_vars];
      var_count = new Int_t[num_vars];
      var_indices = new Int_t[num_vars];
      var_type_size = new Int_t[num_vars];
      
      for (Int_t i=0;i<num_vars;i++){
	if (i==0) temp_entry = strtok(header,";");
	else temp_entry = strtok(NULL,";");;
	if (swap_needed) full_endian_swap(temp_entry);
	var_names[i]=temp_entry;
	//printf("variable %i: %s\n",i,temp_entry);
	temp_entry = strtok(NULL,";");
	//if (swap_needed) full_endian_swap(temp_entry);
	var_types[i]=temp_entry;
	var_type_size[i] = stupid_method(temp_entry);
	temp_entry = strtok(NULL,";");
	if ((strchr_val = strchr(temp_entry,','))!=NULL){
	  sub_str = new Char_t;
	  num_len = strchr_val - temp_entry;
	  memcpy(sub_str,temp_entry,num_len);
	  //if (swap_needed) full_endian_swap(sub_str);
	  var_count[i]=atoi(sub_str);
	  memcpy(sub_str,temp_entry+2,num_len+2);
	  //if (swap_needed) full_endian_swap(sub_str);
	  var_indices[i]=atoi(sub_str);
	  delete[] sub_str;
	}
	else{
	  //if (swap_needed) full_endian_swap(temp_entry);
	  var_count[i]=atoi(temp_entry);
	  var_indices[i]=1;
	}
	//printf("var %i %s: %i by %i\n",i,var_names[i],var_count[i],var_indices[i]);
	
      }
      infile.read((char*)&num_lines,sizeof(Int_t));
      if (swap_needed) full_endian_swap(num_lines);
      printf("There are %i lines in the file.\n",num_lines);
      temp_char = new Char_t*[num_vars];
      
      //FIXME : using dimensions of 50 to determine the number of pulses per event
      //is a shitty kludge that could fail at any time.
      //temp_int = new UInt_t*[var_count[1]];
      //temp_dbl = new Double_t*[var_count[1]];
      temp_int = new UInt_t*[50];
      temp_dbl = new Double_t*[50];
      temp_int_array = new UInt_t**[50];
      temp_dbl_array = new Double_t**[50];
      for (Int_t i=0;i<50;i++){
	temp_int[i] = new UInt_t[num_vars];
	temp_dbl[i] = new Double_t[num_vars];
	temp_dbl_array[i] = new Double_t*[50];
	for (Int_t m=0;m<50;m++) temp_dbl_array[i][m] = new Double_t[num_vars];
      }
      //************************************

      for (Int_t i=0;i<num_lines;i++){
	//printf("Event %i\n",i);
	if (block==2) if (i%50000==0) printf("%3.1f%% finished\n",(100*Double_t(i)/Double_t(num_lines)));
	pulse_dbl_ct=0;
	pulse_int_ct=0;
	event_dbl_ct=0;
	event_int_ct=0;
	channel_dbl_ct=0;
	channel_int_ct=0;
	for (Int_t j=0;j<num_vars;j++){
	  //printf("Var #%i\n",j);
	  temp_char[j] = new Char_t[var_count[j]];
	  //printf("Declarations finished for variable %i\n",j);
	  if (!strcmp(var_types[j],"char")) {
	    infile.read((char*)temp_char[j],var_type_size[j]*var_count[j]);
	    cout << "This is file " << temp_char[j] << endl;
	    /*if(block >1){
	      if (swap_needed) full_endian_swap(temp_char[j]);
	      }*/ // I don't think char's ever need this swap so removing it (could fail at random PHP 2010-10-16
	  }
	  else if (!strcmp(var_types[j],"uint32")){
	    if (var_indices[j]==1) {
	      for (Int_t k=0;k<var_count[j];k++){
		infile.read((char*)&temp_int[k][pulse_int_ct],var_type_size[j]);

		if (swap_needed) full_endian_swap(temp_int[k][pulse_int_ct]); //we should always need this swap since data should be binary
		temp_int[k][pulse_int_ct]=UInt_t(temp_int[k][pulse_int_ct]);
		//printf("read temp_int[%i][%i] = %i\n",k,pulse_int_ct,temp_int[k][pulse_int_ct]);
	      }
	      pulse_int_ct++;
	    }
	    else {
	      //printf("integer array array start\n");
	      int_array_def_flag = 1;
	      temp_int_array[channel_int_ct] = new UInt_t*[var_indices[j]];
	      for (Int_t m=0;m<var_count[j];m++){
		//temp_int_array[channel_int_ct][m] = new UInt_t[var_indices[j]];
		for (Int_t k=0;k<var_indices[j];k++) {
		  temp_int_array[channel_int_ct][k] = new UInt_t[num_vars];
		  infile.read((char*)&temp_int_array[channel_int_ct][k][m],var_type_size[j]);
		  if (swap_needed) full_endian_swap(temp_int_array[channel_int_ct][k][m]);
		  temp_int_array[channel_int_ct][k][m]=UInt_t(temp_int_array[channel_int_ct][k][m]);
		}
		n_pulses = var_count[j];
		if (var_indices[j]!=1) n_channels = var_indices[j];
	      }
	      channel_int_ct++;
	      //printf("integer array array done\n");
	    } 
	  }
	  else if (!strcmp(var_types[j],"double")){
	    if (var_indices[j]==1) {
	      //if (var_count[j]>5) var_count[j]=5;
	      for (Int_t k=0;k<var_count[j];k++){
		infile.read((char*)&temp_dbl[k][pulse_dbl_ct],var_type_size[j]);
		if (swap_needed) full_endian_swap(temp_dbl[k][pulse_dbl_ct]); //Should need this swap
		temp_dbl[k][pulse_dbl_ct]=Double_t(temp_dbl[k][pulse_dbl_ct]);
		//printf("read temp_dbl[%i][%i] = %.2f\n",k,pulse_dbl_ct,temp_dbl[k][pulse_dbl_ct]);
	      }
	      pulse_dbl_ct++;
	    }
	    else {
	      //printf("inside the double array array portion\n");
	      dbl_array_def_flag = 1;
	      for (Int_t m=0;m<var_indices[j];m++){
		for (Int_t k=0;k<var_count[j];k++) {
		  infile.read((char*)&temp_dbl_array[channel_dbl_ct][m][k],var_type_size[j]);
		  if (swap_needed) full_endian_swap(temp_dbl_array[channel_dbl_ct][m][k]);
		  temp_dbl_array[channel_dbl_ct][m][k]=Double_t(temp_dbl_array[channel_dbl_ct][m][k]);
		  //printf("read temp_dbl_array[%i][%i][%i] = %.2f\n",channel_dbl_ct,m,k,temp_dbl_array[channel_dbl_ct][m][k]);

		}		 
		n_pulses = var_count[j];
		n_channels = var_indices[j];
	      }
	      channel_dbl_ct++;
	    }
	  }
	  else {
	    printf("I don't recognize the type %s!\n",var_types[j]);
	    printf("You're screwed.  Good luck with all of that.\n");
	  }
	}
	//printf("Through the reading, into the writing of event %i\n",i);
       
	if (block==1) lrq1_header->Fill(temp_char,temp_int[0],temp_dbl[0]);
	else if (block==2) {
	  if (n_pulses>5) printf("Event %i has %i pulses!\n",i,n_pulses);
	  
	  lrq1_event->Set_num_pulses(n_pulses);
	  lrq1_event->Set_num_channels(n_channels);
	  //printf("There are %i pulses\n",n_pulses);
	  for (Int_t i=0;i<n_pulses;i++){
	    lrq1_pulse->Fill(temp_char,temp_int[i],temp_dbl[i]);
	    lrq1_event->AddPulse(lrq1_pulse,i);
	    for (Int_t j=0;j<n_channels;j++){
	      //fill channel here...
	      if (!int_array_def_flag) {
		//printf("Need to define temp_int_array\n");
		temp_int_array[j]=new UInt_t*[1];
		temp_int_array[j][i]=new UInt_t[1];
		temp_int_array[j][i][0]=0;
		//printf("temp_int_array defined\n");
	      }
	      if (!dbl_array_def_flag){
		//printf("Need to define temp_dbl_array\n");
		temp_dbl_array[j]=new Double_t*[1];
		temp_dbl_array[j][i]=new Double_t[1];
		temp_dbl_array[j][i][0]=0;
	      }
	      index = i*n_channels+j;
	      lrq1_channel->Fill(temp_char,temp_int_array[j][i],temp_dbl_array[j][i]);
	      lrq1_event->AddChannel(lrq1_channel,index);
	      lrq1_channel->Clear();
	    }
	  }
	}
	//printf("Before filling tree\n");
	LUX_rq1_tree->Fill();
	//printf("done filling tree\n");
	lrq1_event->Clear();
      }
    }
    rootfile->cd();
    //printf("Before writing ROOT file\n");
    rootfile->Write();
    printf("Before closing ROOT file\n");
    rootfile->Close();  
    printf("ROOT file written and closed\n");
    
    delete[] header;
    delete[] var_names;
    delete[] var_types;
    delete[] var_count;
    delete[] var_indices;
    delete[] var_type_size;
    delete[] temp_char;
    delete[] temp_int;
    delete[] temp_int_array;
    delete[] temp_dbl_array;  
    delete[] sub_str;
    
  }
  else {
    cout << "File could not be opened" << endl;
  }
  
  return;
}

int main(int argc, char **argv)
{
  int result=0;
  if (argc > 1){
    read_rq1(argv[1]);
  }   
  else{
    printf("Usage: read_rq1_built <binary file>");
    printf(" where <binary_file> is the name of the merged file");
    printf(" to be read\n");
  }
  return (result==0)?EXIT_SUCCESS:EXIT_FAILURE;

}
