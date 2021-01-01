//////////////////////////////////////////////////////////////////////////////
// Dump an rq1 file to ascii.
//
//  Michael Woods
//
//////////////////////////////////////////////////////////////////////////////


#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <cstring>
#include <string>
//#include <TFile.h>
//#include <TTree.h>
//#include <TString.h>
//#include <Riostream.h>

using namespace std;


// A shortcut to the endian byteswapping
#define EndianByteSwap(x) ByteSwap((unsigned char *) &x,sizeof(x))

// The following two help read in binary data without cluttering the screen.
#define readval(x) read((char*)&x,sizeof(x));if(byteswap)EndianByteSwap(x)
#define readstr(x,y) read((char*)x,y);if(byteswap)EndianByteSwap(x)

// Make a useful function for colorizing outout via cout (depends on BLD and
// NRM being defined in the code, though.
#define COL(i) BLD << i << NRM

// Easy two dimensional arrays
#define CREATEARRAY_2D(TYPE,NAME,XSIZE,YSIZE)TYPE **NAME = new TYPE*[(XSIZE)]; for (int tempi = 0; tempi < (XSIZE); tempi++) NAME[tempi] = new TYPE[(YSIZE)];
#define DELETEARRAY_2D(NAME,XSIZE) for (int tempi = 0; tempi < (XSIZE); tempi++) delete[] NAME[tempi]; delete[] NAME;

//////////////////////////////////////////////////////////////////////////////
//              void syntax()
//////////////////////////////////////////////////////////////////////////////
void syntax() {
    cout << "Syntax:" << endl;
    cout << "   ./rq1_to_ascii [-c] [-e #] lux10_filename.rq1" << endl;
    cout << "" << endl;
    cout << "   -c      Turn coloring of output off." << endl;
    cout << "   -e {#}  Print only a given event number." << endl;
    cout << "   -h      Show this help." << endl;
}

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
//////////////////////////////////////////////////////////////////////////////
//             vector<string> DelimitStringToVec
//////////////////////////////////////////////////////////////////////////////
// Turn "area;double;1;name;char;15" into:
//      vector<string> { "area", "double", "1", "name", "char", "15" }
// Thanks Surge.
//////////////////////////////////////////////////////////////////////////////
void DelimitStringToVec(string str, vector<string>* blocks, 
                        string separator=";"){
    blocks->clear();
    string headerstring = str;

    size_t pos1 = 0, pos2 = 0;
    size_t length = headerstring.length();
    while ( pos1 < length && pos2 < headerstring.npos ){
            pos2 = headerstring.find_first_of(separator, pos1);
            if (pos2 != pos1){
                string sub = headerstring.substr(pos1, pos2-pos1);
                blocks->push_back( sub );
            }
            pos1 = pos2+1;
    }

    return ;
}
//////////////////////////////////////////////////////////////////////////////
// Print out an array given its size.
//////////////////////////////////////////////////////////////////////////////
template <typename AnyType>
void PrintRQ(ifstream &infile, AnyType rq_datum, vector<int> dimensions,
             bool byteswap, bool display, bool color=true) {
  char BLD[10] = "";
  char NRM[10] = "";

  if(color) {
    sprintf(BLD, "%c[1m", 0x1B);
    sprintf(NRM, "%c[0m", 0x1B);
  }  
  
  int length = 1;
  for (size_t d=0; d<dimensions.size(); d++) {
    length *= dimensions[d];
    if (display)
      cout << COL("[") << COL(dimensions[d]) << COL("]");
  }
  if (dimensions.size() == 1 && display && dimensions[0] == 1) cout << " = ";
  else if (display) cout << endl;
  
  for (int i=1; i<=length; i++) {
    infile.readval(rq_datum);
    if (!display) continue;
    
    cout.width(8);
    cout << boolalpha;
    cout << rq_datum << "\t";    
    for (size_t d=0; d<dimensions.size(); d++) {
      int mod = dimensions[d];
      for (size_t dd=d+1; dd<dimensions.size(); dd++)
        mod *= dimensions[dd];
      if (i % mod == 0) cout << "\n";
    }
  }  
}
void PrintRQ(ifstream &infile, char rq_datum, vector<int> dimensions,
                   bool byteswap, bool display, bool color=true) {
  char BLD[10] = "";
  char NRM[10] = "";

  if(color) {
    sprintf(BLD, "%c[1m", 0x1B);
    sprintf(NRM, "%c[0m", 0x1B);
  }  
  
  int length = 1;
  for (size_t d=0; d<dimensions.size(); d++) {
    length *= dimensions[d];
    if (display)
      cout << COL("[") << COL(dimensions[d]) << COL("]");
  }
  if (dimensions.size() == 1 && display) cout << " = ";
  else if (display) cout << endl;
  
  for (int i=1; i<=length; i++) {
    infile.readval(rq_datum);
    if (!display) continue;
    
    cout << rq_datum;    
    for (size_t d=0; d<dimensions.size(); d++) {
      int mod = dimensions[d];
      for (size_t dd=d+1; dd<dimensions.size(); dd++)
        mod *= dimensions[dd];
      if (i % mod == 0) cout << "\n";
    }
  }
}
//////////////////////////////////////////////////////////////////////////////
// Function that finds data types and dimensions before calling the printer
//////////////////////////////////////////////////////////////////////////////
void RQPrintSetup(ifstream &infile, string name, string type, string dim, 
                  bool byteswap, bool display, bool color=true) {
  
  // get the dimension vector 
  vector<string> temp_dim;
  vector<int> dimensions;
  DelimitStringToVec(dim, &temp_dim, ",");
  for (size_t d=0; d<temp_dim.size(); d++) {
    dimensions.push_back( atoi(temp_dim[d].c_str()) );
  }
  
  if (display) cout << name;  
  // find the type && call PrintRQ to print out the values 
  if(type=="char" || type=="uint8") {
    char rq_datum = 0;      
    PrintRQ(infile, rq_datum, dimensions, byteswap, display, color);
  }
  else if(type =="int" || type=="int32") {
    int rq_datum = 0; 
    PrintRQ(infile, rq_datum, dimensions, byteswap, display, color);
  }
  else if(type=="unsigned int" || type=="uint32") {
    unsigned int rq_datum = 0;  
    PrintRQ(infile, rq_datum, dimensions, byteswap, display, color);
  }
  else if(type=="long long" || type=="int64") {
    long long rq_datum = 0; 
    PrintRQ(infile, rq_datum, dimensions, byteswap, display, color);
  }	    
  else if(type=="unsigned long long" || type == "uint64") {
    unsigned long long rq_datum = 0;  
    PrintRQ(infile, rq_datum, dimensions, byteswap, display, color);
  }
  else if(type=="double") {
    double rq_datum = 0;   
    PrintRQ(infile, rq_datum, dimensions, byteswap, display, color);
  }
  else if(type=="float" || type == "single") {
    float rq_datum = 0;   
    PrintRQ(infile, rq_datum, dimensions, byteswap, display, color);
  }
  else if(type=="bool") {
    bool rq_datum = 0;   
    PrintRQ(infile, rq_datum, dimensions, byteswap, display, color);
  }
  else if(type=="short" || type=="int16") {
    short rq_datum = 0;   
    PrintRQ(infile, rq_datum, dimensions, byteswap, display, color);
  }
  else if(type=="unsigned short" || type=="uint16") {
    unsigned short rq_datum = 0;   
    PrintRQ(infile, rq_datum, dimensions, byteswap, display, color);
  }
  else {
    cerr << "Error: Unknown type '" << type << "' in rq file\n";
    return;
  }
}

//////////////////////////////////////////////////////////////////////////////
// Our main printing program.
//////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv){
    
    // Set color string for parsing. BLD=bold, NRM=no color/formatting.
    char BLD[10] = "";
    char NRM[10] = "";
    bool color = false;

    // Let's parse some runtime arguments. The instructions should be:
    //  
    //  Syntax:
    //          ./rq1_to_ascii [-c] [-e #] lux10_0000000000000_f000000000.rq1
    //
    //          -c      Turn coloring of output off.
    //          -e {#}  Print only a given event number.
    //          -h      Show this help.

    int only_one_event = -1;
    for(int arg=0; arg<argc; arg++) {
        // Color
        if(!strcmp(argv[arg], "-c")) {
            color = true;
            sprintf(BLD, "%c[1m", 0x1B);
            sprintf(NRM, "%c[0m", 0x1B);
        }
        // Help
        if(!strcmp(argv[arg], "-h")) {
            syntax();
            return 0;
        }
        // Syntax
        if(argc == 1 || argc > 5) {
            syntax();
            return 0;
        }
        // Event number
        if(!strcmp(argv[arg], "-e")) {
            only_one_event = atoi(argv[arg+1]);
            arg++;
        }
    }
    
    // Open rq1 file 
    char* infilename = argv[argc-1];
    ifstream infile;    
    infile.open (infilename, ios::binary);
    if(!infile.is_open()) {
        cerr << "Error opening rq1 file." << endl;
        return 0;
    }
    
    // Proclaim what you're doing.
    cout << BLD << "You are going to be getting a large ascii dump of the contents";
    cout << "of the following file: " << endl << infilename << NRM << endl;
    cout << COL("/------------------------------------\\") << endl;
    cout << COL("|              FILE HEADINGS         |") << endl;
    cout << COL("\\------------------------------------/") << endl;

    
    // Check the Endianness of the file.
    unsigned int Endianness;
    infile.read( (char*)&Endianness, sizeof(Endianness) );        
    // Set the byteswap boolean if needed.
    bool byteswap = false;
    //if (Endianness == 0x01020304) byteswap = false;
    cout << COL("Endianness: \t") << Endianness << endl;
    
    // Read XML Settings (length of the string, then full text).
    unsigned int xml_settings_string_length;
    infile.readval(xml_settings_string_length); 
    char* xml_settings_string = new char[xml_settings_string_length+1];
    infile.readstr(xml_settings_string, xml_settings_string_length);
    xml_settings_string[xml_settings_string_length] = '\0';
    cout << COL("---------------XML Information--------------------") << endl;
    cout << COL("XML Settings String Length: \t");
    cout << xml_settings_string_length << endl;
    cout << COL("XML Settings String:        \n");
    cout << xml_settings_string << endl;
    cout << COL("------------------- END XML ----------------------") << endl;
    delete[] xml_settings_string;
    
    
    // RQ Block Format Contaniers
    unsigned short header_string_size = 0;
    char * header_cstring;
    string header_string;
    int header_number_of_lines = 0;
    vector<string> header_string_vec;
    
    //////////////////////////////
    // Read in the RQ Header Block
    //////////////////////////////
    cout << COL("/------------------------------------\\") << endl;
    cout << COL("|           RQ HEADER BLOCK           |") << endl;
    cout << COL("\\------------------------------------/") << endl;
    
    // Read in the header string
    infile.readval(header_string_size);
    header_cstring = new char[header_string_size+1];
    infile.readstr(header_cstring, header_string_size);
    header_cstring[header_string_size] = '\0';
    header_string = header_cstring;
    delete [] header_cstring;
    infile.readval(header_number_of_lines);
    cout << COL("Header string length: \t") << header_string_size << endl;
    cout << COL("Header string:        \t") << header_string << endl;
    cout << COL("Number of data lines: \t") << header_number_of_lines << endl;

    // Parse the header string so we know what we're reading in.    
    DelimitStringToVec(header_string, &header_string_vec);
    if (header_string_vec.size() % 3 != 0) {  
      cerr << "Oh crap, multidimensional header data broke me." << endl;
      cerr << "Please let Mike know." << endl;
      return 0;
    }
    // Loop over each third entry in the above vector.
    for(size_t i=0; i<header_string_vec.size() ; i+=3) {
      RQPrintSetup(infile, header_string_vec[i], header_string_vec[i+1],
                   header_string_vec[i+2], byteswap, true, color);
    }
    //-----------------------------------------------------------------

    //////////////////////////////
    // Read in the RQ Events Block
    //////////////////////////////
    // Ok, now let's sit down and have a cup of coffee, a glass of whiskey,
    // and talk about how the rq's are going to be read and displayed to
    // screen. Let's discuss what kind of rq's can have what dimensions.
    // 0D RQ:   event_number - One number, no array needed.
    // 1D RQ:   t10l -
    //      The time it takes for the summed waveform of all channels to go
    //      above 10% on the left side of the max height of the pulse. Need a
    //      1D for each pulse in an event (e.g. 2 S1's and 2 S2's needs a
    //      1D array with a length of four.
    // 2D RQ:   baseline_mV -
    //      A baseline noise value for each channel (1D), for each peak in
    //      each channel (1D). Two dimensions with lengths (just as an exampl)
    //      of [4,122] (4 peaks, 122 chans).
    // 3D RQ:   channel_fits -
    //      You fit each channel (1D) for each of its peaks (1D) a gaussian
    //      with three parameters (1D). [4,122,3]. I don't know how often
    //      we'll see one of these.
    //
    // Ok, now these are going to be printed out with the 1st dimension going
    // across the screen horizontally, the second dimension going down the
    // screen vertically, and (if need be but not yet implemented) the third
    // dimension just being a repeat of the 2D grid from previously.
    //
    // The variable array_size is a 3D integer array that starts off as simply
    // {1,1,1} which represents a 0D array. Based on the parsing the
    // dimensions we need are saved. This is done with try/catch statements

    cout << COL("/------------------------------------\\") << endl;
    cout << COL("|           RQ EVENTS BLOCK           |") << endl;
    cout << COL("\\------------------------------------/") << endl;
    
    // Read in the header string
    infile.readval(header_string_size);
    header_cstring = new char[header_string_size+1];
    infile.readstr(header_cstring, header_string_size);
    header_cstring[header_string_size] = '\0';
    header_string = header_cstring;
    delete [] header_cstring;
    infile.readval(header_number_of_lines);
    cout << COL("Header string length: \t") << header_string_size << endl;
    cout << COL("Header string:        \t") << header_string << endl;
    cout << COL("Number of data lines: \t") << header_number_of_lines << endl;
    
    // Parse the header string so we know what we're reading in.    
    DelimitStringToVec(header_string, &header_string_vec);
    if (header_string_vec.size() % 3 != 0) {  
      cerr << "Oh crap, multidimensional header data broke me." << endl;
      cerr << "Please let Mike know." << endl;
      return 0;
    }
    
    // Print out the rqs 
    for(int event_num=1; event_num<=header_number_of_lines; event_num++) {
      char title_buf[256];
      bool display = true;
      if( only_one_event != -1 && event_num != only_one_event)
        display = false;
      sprintf(title_buf, "==Event %d from first_evt_in_file==", event_num);
      if (display) cout << COL(title_buf) << endl;
      
      for(size_t i=0; i<header_string_vec.size() ; i+=3)
        RQPrintSetup(infile, header_string_vec[i], header_string_vec[i+1],
                    header_string_vec[i+2], byteswap, display, color);
    }

    //-----------------------------------------------------------------

    //////////////////////////////
    // Read in the Livetime Block
    //////////////////////////////    
    cout << COL("/------------------------------------\\") << endl;
    cout << COL("|         RQ LIVETIME BLOCK           |") << endl;
    cout << COL("\\------------------------------------/") << endl;
    
    // Read in the header string
    infile.readval(header_string_size);
    header_cstring = new char[header_string_size+1];
    infile.readstr(header_cstring, header_string_size);
    header_cstring[header_string_size] = '\0';
    header_string = header_cstring;
    delete [] header_cstring;
    infile.readval(header_number_of_lines);
    cout << COL("Header string length: \t") << header_string_size << endl;
    cout << COL("Header string:        \t") << header_string << endl;
    cout << COL("Number of data lines: \t") << header_number_of_lines << endl;

    // Parse the header string so we know what we're reading in.    
    DelimitStringToVec(header_string, &header_string_vec);
    if (header_string_vec.size() % 3 != 0) {  
      cerr << "Oh crap, multidimensional header data broke me." << endl;
      cerr << "Please let Mike know." << endl;
      return 0;
    }
    
    // Print out the rqs 
    for(int e=1; e<=header_number_of_lines ; e++) {      
      for(size_t i=0; i<header_string_vec.size() ; i+=3)
        RQPrintSetup(infile, header_string_vec[i], header_string_vec[i+1],
                    header_string_vec[i+2], byteswap, true, color);
    }
    return 0;    
}

