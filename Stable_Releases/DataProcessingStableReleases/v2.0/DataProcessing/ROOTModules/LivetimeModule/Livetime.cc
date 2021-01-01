/////////////////////////////////////////////////////////////////////
// This is a very basic module that will probably be superceded by
// a full "basic rq" module. It's purpose to determine the livetime
// of an evt file and output it to screen.
//
// LOG
// 2012: Dec 4 - MW - Creation
// 2012: Dec 18 - JRV - Appended .lvtm to filename
//
/////////////////////////////////////////////////////////////////////


#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <cstring>
#include <string>

// A shortcut to the endian byteswapping
#define EndianByteSwap(x) ByteSwap((unsigned char *) &x,sizeof(x))

// The following two help read in binary data without cluttering the screen.
#define readval(x) read((char*)&x,sizeof(x));if(byteswap)EndianByteSwap(x)
#define readstr(x,y) read((char*)x,y);if(byteswap)EndianByteSwap(x)

using namespace std;

//////////////////////////////////////////////////////////////////////////////
//              void ByteSwap(unsigned char * b, int n)
//////////////////////////////////////////////////////////////////////////////
void ByteSwap(unsigned char * b, int n){
   // This function swithes between big endianness to little endianness and 
   //vise versa taken from http://www.codeproject.com/KB/cpp/endianness.aspx

   register int i = 0;
   register int j = n-1;
   while (i<j){
      std::swap(b[i], b[j]);
      i++, j--;
   }
}

int main(int argc, char* argv[]) {

    if( argc == 1 ) {
      cerr << "Please specify a .evt file" << endl;
      return -1;
    }

    // Open evt file 
    string infilename_path = argv[1];
    string infilename = argv[2];
    infilename = infilename_path + infilename;
    ifstream infile;
    infile.open (infilename.data(), ios::binary);
    if(!infile.is_open()) {
        cerr << "Error opening evt file " << infilename << endl;
        return 0;
    }



    // Check the Endianness of the file.
    unsigned int Endianness;
    infile.read( (char*)&Endianness, sizeof(Endianness) );
    // Set the byteswap boolean if needed.
    bool byteswap = false;
    if (Endianness == 0x01020304) byteswap = false;
    if (Endianness == 0x04030201) byteswap = true;

    // Read XML Settings (length of the string, then full text).
    unsigned int xml_settings_string_length;
    infile.readval(xml_settings_string_length);
    char* xml_settings_string = new char[xml_settings_string_length+1];
    infile.readstr(xml_settings_string, xml_settings_string_length);
    xml_settings_string[xml_settings_string_length] = '\0';
    delete[] xml_settings_string;

    ///////////////////////////
    // Read in File Header[5]
    ///////////////////////////
    // Read in endianness, again.
    infile.readval(Endianness);
    // Read in date and time of run.
    unsigned int datetime;
    infile.readval(datetime);
    // Read in location.
    unsigned int location;
    infile.readval(location);
    // Read in number of events in the file.
    unsigned int num_events_in_file;
    infile.readval(num_events_in_file);


    // Read in byte locations of the file.
    for(unsigned int i=0; i<num_events_in_file; i++) {
      // Read in the event number.
      unsigned int event_number;
      unsigned int byte_location;
      infile.readval(event_number);
      infile.readval(byte_location);
    }

    ///////////////////////////
    // Read in Live Time Header
    ///////////////////////////
    // Read in number of sequences in file.
    unsigned short num_seqs;
    infile.readval(num_seqs);

    // Read in each timestamp latch.
    unsigned long long total_livetime = 0;
    for(unsigned int i=0; i<num_seqs; i++) {
      unsigned long long timestamp_latch;
      unsigned long long timestamp_end;
      infile.readval(timestamp_latch);
      infile.readval(timestamp_end);
      total_livetime += timestamp_end - timestamp_latch;
    }

  cout  << total_livetime;

  string outfile_name = argv[4];
  outfile_name += argv[3];
  outfile_name.append(".lvtm");

  ofstream out_file;
  out_file.open(outfile_name.data(), ios::out);
  out_file << total_livetime << endl;
  out_file.close();
  return 0;
} 

