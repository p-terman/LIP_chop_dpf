/*
This program takes the .rq1 file produce by the LUX_First_Pass and converts
the files into the .rq1.root format. This program is a stand alone and not meant
to run interactively in the root environment. Formerly called "rq1_to_root"

2013-04-15 AC Updated to be called as a module in the DP framework, making a 
directory called rootfiles if needed saving to it

LUX Collaboration
Group: ROOT Analysis Chain
Author: Sergey Uvarov
Version: 2.1
*/

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <TObject.h>
#include <TObjString.h>
#include <TTree.h>
#include <TFile.h>
#include <TSystem.h>

using namespace std;

#define readval(x) read((char*)&x,sizeof(x));if(byteswap)EndianByteSwap(x)
#define readstr(x,y) read((char*)x,y);if(byteswap)ByteSwap((unsigned char*)x,y)
#define EndianByteSwap(x) ByteSwap((unsigned char *) &x,sizeof(x))

void ByteSwap(unsigned char * b, int n) { 
  // This function swithes between big endianness to little endianness and vise
  // versa taken from http://www.codeproject.com/KB/cpp/endianness.aspx    
  register int i = 0;
  register int j = n-1;
  while (i<j){
    std::swap(b[i], b[j]);
    i++, j--;
  }
}
//______________________________________________________________________________
vector<string> HeaderBlocks(string str, const char* separator) {
	vector<string> blocks;
	size_t pos1 = 0, pos2 = 0;
	do {
	  pos2 = str.find(separator, pos1);
	  if (pos1 != string::npos and pos1 < str.size()) {
	    string sub = str.substr(pos1, pos2-pos1);
	    blocks.push_back(sub);
	  }
	  pos1 = pos2+1;
	}while (pos1 != string::npos and pos2 != string::npos);
	return blocks;
}
//______________________________________________________________________________
int RowToColumn(int linear_col, vector<int> dims) {
  // gets the linearized row-major index from the linearized column-major array
  int linear_row = 0;
  for (size_t d=0; d<dims.size(); d++) {
    int divide_dim = 1;
    for (size_t i=0; i<d; i++) divide_dim *= dims[i];
    int multiply_dim = 1;
    for (size_t i=d+1; i<dims.size(); i++) multiply_dim *= dims[i];    
    linear_row += ((linear_col/divide_dim) % dims[d]) * multiply_dim;
  }
  return linear_row;
}
//______________________________________________________________________________
//                             Block Class
//------------------------------------------------------------------------------
// This class handles one rq value. It determines its data type, its dimensions,
// and name. It will also read in its values and added the address to the TTree.
//______________________________________________________________________________
class Block {
	private:
		char* raw;
		void DetermineType(string );
		void DetermineDim(string );
	
	public:		
		Block();
		Block(string, string, string);
		~Block();
		
		string name;
		char type;
		int bytes;
		int lineardim;
		vector<int> dims;	
		
		void ReadFile(istream &, bool);
		void Display() const;
		void AddBranchToTree( TTree* );
};
//______________________________________________________________________________
Block::Block() {
	// empty constructor
	name=""; 
	type = 'B';
	lineardim = 0;
	bytes = 0;
	raw = NULL;
}
//______________________________________________________________________________
Block::Block(string n, string t, string d) {
	// constructor that builds the rq information
	//cout << "Created RQ: ";
	name = n; 
	//cout << name << " | ";
	type = 'B';
	lineardim = 0;
	bytes = 0;
	raw = NULL;
	DetermineType(t);
	//cout << type << " | " << bytes << " | ";
	DetermineDim(d);
	//for (size_t i=0; i<dims.size(); i++)  cout <<  dims[i] << " ";
	//cout << "| " << lineardim << " ";
	raw = new char[lineardim*bytes+1];
	raw[lineardim*bytes] = '\0';
	//cout << "done\n";
}
//______________________________________________________________________________
Block::~Block() {
	// destructor
	//cout << "delete RQ " << name << " line 112\n";
	//if (raw and lineardim) delete[] raw;
}
//______________________________________________________________________________
void Block::Display() const {
	// This function prints to screen the rq information 
	char *temp = new char[bytes+1];
	temp[bytes] = '\0';
	
	cout << name << " | " << type << " | " << bytes << " | ";
	for (size_t i=0; i<dims.size(); i++) cout << dims[i] << " ";
	cout << endl;
	
	for (int i=0; i<lineardim; i++){
		// copy over data
		int offset = i * bytes;
		for (int j=0; j<bytes; j++)
			temp[j] = raw[offset+j];

		if (type == 'O')
			cout << *(bool*)temp << " ";
		else if (type == 'C')
			cout << (const char*)temp << " ";
		else if (type == 'B')
			cout << *(char*)temp << " ";
		else if (type == 'b')
			cout << *(unsigned char*)temp << " ";
		else if (type == 'S')
			cout << *(short*)temp << " ";
		else if (type == 's')
			cout << *(unsigned short*)temp << " ";
		else if (type == 'I')
			cout << *(int*)temp << " ";
		else if (type == 'i')
			cout << *(unsigned int*)temp << " ";
		else if (type == 'L')
			cout << *(long long*)temp << " ";
		else if (type == 'l')
			cout << *(unsigned long long*)temp << " ";
		else if (type == 'F')
			cout << *(float*)temp << " ";
		else if (type == 'D')
			cout << *(double*)temp << " ";
		else 
			cout << (const char*)(raw) << endl;
	}
	delete[] temp;
	cout << endl;
}
//______________________________________________________________________________
void Block::DetermineType(string str) {
	if ( !str.compare("bool") || !str.compare("Bool_t")
			|| !str.compare("logical") ){
		type = 'O';
		bytes = 1;
	}
	else if ( !str.compare("string")) {
	  type = 'C';
	  bytes = 1;
	}
	else if ( !str.compare("char") || !str.compare("Char_t") ){
		type = 'B';
		bytes = 1;
	}
	else if ( !str.compare("unsigned char") || !str.compare("UChar_t")
			|| !str.compare("uint8")){
		type = 'b';
		bytes = 1;
	}
	else if ( !str.compare("short") || !str.compare("Short_t")
			|| !str.compare("int16") ){
		type = 'S';
		bytes = 2;
	}
	else if ( !str.compare("unsigned short") || !str.compare("UShort_t")
			|| !str.compare("uint16") ){
		type = 's';
		bytes = 2;
	}
	else if ( !str.compare("int") || !str.compare("Int_t")
			|| !str.compare("int32")){
		type = 'I';
		bytes = 4;
	}
	else if ( !str.compare("unsigned int") || !str.compare("uint")
						|| !str.compare("uint32")){
		type = 'i';
		bytes = 4;
	}
	else if ( !str.compare("long long") || !str.compare("Long_t")
			|| !str.compare("int64") || !str.compare("long long int")){
		type = 'L';
		bytes = 8;
	}
	else if ( !str.compare("unsigned long long") || !str.compare("ULong_t")
						|| !str.compare("uint64") ){
		type = 'l';
		bytes = 8;
	}
	else if ( !str.compare("float") || !str.compare("Float_t")
			|| !str.compare("single")){
		type = 'F';
		bytes = 4;
	}
	else if ( !str.compare("double") || !str.compare("Double_t")){
		type = 'D';
		bytes = 8;
	}
	else{
		type = 'B';
		bytes = 1;
	}
}
//______________________________________________________________________________
void Block::DetermineDim(string str) {
  // Gets 
  dims.clear();
  vector<string> dimensions = HeaderBlocks(str, ",");
  lineardim = 1;
  for (size_t i=0; i<dimensions.size(); i++) {
    dims.push_back( atoi(dimensions[i].c_str()) );
    lineardim *= dims[i];
  }
}
//______________________________________________________________________________
void Block::ReadFile(istream &file, bool byteswap) {
	char* temp = new char[bytes];
	for (int i=0; i<lineardim; i++) {
		file.read( temp, bytes );
		if (byteswap)	ByteSwap( (unsigned char*)temp, bytes);
		int offset = RowToColumn(i, dims)*bytes;
		//cout << "RQ offsets " << i << "->" << offset << " ";
		for (int j=0; j<bytes; j++){			
			raw[offset+j] = temp[j];
		}
	}
	//cout << endl;
	delete[] temp;
}
//______________________________________________________________________________
void Block::AddBranchToTree(TTree* tree) {
	char str[100];
	string temp = name;
	for (size_t i=0; i<dims.size(); i++) {
	  if (dims[i] != 1) {
	    sprintf(str, "[%i]", dims[i]);
	    temp += str;
	  }
	}
	temp += "/";
	temp += type;
	temp += '\0';
	tree->Branch(name.data(), raw, temp.data());
}
//______________________________________________________________________________
//                             RQ Reader class 
//------------------------------------------------------------------------------
// This class handles one rq block. This includes reading the header block 
// information to get all the rqs, and writting those rqs to the TTree. 
//______________________________________________________________________________
class rqReader {
	protected:
		vector<Block> blocks;
		string headerstring;
		TTree *tree;
	public:
		rqReader(string s, string tree_name);
		rqReader(istream &file, string tree_name, bool swap=false);
		~rqReader();
				
		void ReadFile(istream &file, bool swap=false);
		void Write() { tree->Write(); }		
};
//______________________________________________________________________________
rqReader::rqReader(string str, string tree_name) {
	tree = NULL;
	headerstring = str;
	vector<string> headerblock = HeaderBlocks(str, ";" );
	
	// sanity check
	if ( (headerblock.size() % 3) ){ 
		cerr << "Error: Incorrect header string structure\n";
		cerr << "Header string size is " << headerblock.size() << endl;
		exit(0);
	}
	
	for (unsigned int i=0; i<headerblock.size(); i+=3 ){
		Block temp( headerblock[i], headerblock[i+1], headerblock[i+2] 
);
		blocks.push_back(temp);
	}
	// create TTree
	tree = new TTree(tree_name.c_str(), "");
	for (size_t i=0; i<blocks.size(); i++)
		blocks[i].AddBranchToTree(tree);
	
	headerblock.clear();
}
//______________________________________________________________________________
rqReader::rqReader(istream &file, string tree_name, bool byteswap) {
	// set defaults 
	tree = NULL;
	
	// Read the header information
	unsigned short header_string_size;
  file.readval(header_string_size);
  char* header_cstring = new char[header_string_size+1];
  file.readstr(header_cstring, header_string_size);
  header_cstring[header_string_size] = '\0';
  unsigned int header_number_of_lines;
  file.readval(header_number_of_lines);
	
	// get possible rqs 
	headerstring = header_cstring;
	vector<string> headerblock = HeaderBlocks(headerstring, ";");
	
	// sanity check
	if ( (headerblock.size() % 3) ){ 
		cerr << "Error: Incorrect header string structure\n";
		cerr << "Header string size is " << headerblock.size() << endl;
		exit(0);
	}
	
	// build rqs 
	for (unsigned int i=0; i<headerblock.size(); i+=3 ){
		Block temp( headerblock[i], headerblock[i+1], headerblock[i+2] 
);
		blocks.push_back(temp);
	}
	
	//headerblock.clear();
	
	// create TTree and add rq links to TTree 
	tree = new TTree(tree_name.c_str(), "");
	for (size_t i=0; i<blocks.size(); i++)
		blocks[i].AddBranchToTree(tree);
	
	// load rq data into TTree 
	for ( unsigned int i=0; i<header_number_of_lines; i++)
		this->ReadFile(file, byteswap);
	
	// write to root file 
	this->Write();
}
//______________________________________________________________________________
rqReader::~rqReader() {
	// destructor 
	blocks.clear();
	if (tree) delete tree;
	tree = NULL;
}
//______________________________________________________________________________
void rqReader::ReadFile(istream &file, bool byteswap) {
	// load the data
	for (size_t block=0; block<blocks.size(); block++ ){
		blocks[block].ReadFile(file, byteswap);
		//blocks[block].Display();
	}
	this->tree->Fill();	
}
//______________________________________________________________________________
//                            Main Function              
//______________________________________________________________________________
int main(int argc, char *argv[]) {
	
	for (int arg=1; arg<argc; arg++) {
	  string filename;
	  string outdir;
	  bool inDpFramework=0;
	  if (argc==8){
      //looks like you're running DPFramework and so file,dir are args 3,4 
      //when someone tries to rootify 8 files and complains, check other args
	    inDpFramework=1;
	    arg=argc; // only one file to be processed...
	    filename=string(argv[4])+string("/")+string(argv[3]); // ... this one.
	    outdir=string(argv[4])+string("/rootfiles");//.roots go here
	    gSystem->mkdir(outdir.c_str()); // if not made, make a rootfiles dir
	    outdir.append("/");
	  }
	  else{
	    //looks like you're not running in DP framework
	    filename = argv[arg];
	  }
	  // get the input file 
    string extension = filename;
    extension.erase(0, extension.rfind(".rq"));
    if (extension.compare(".rq1") and extension.compare(".rq2")
      and extension.compare(".rq")){
      cerr << "Error: "<< filename  << " is not a .rq file\n";  
      continue;
    }
    
    // open rq1 file 
	  ifstream infile(filename.c_str(), ios::binary);
	  cout<<filename<<endl;
	
	  // create root file 
	  string outfilename = filename;
	  if(inDpFramework){
	    outfilename=outdir+string(argv[3]);
	    cout<<outfilename<<endl;
	  }
	  outfilename.append(".root");
	  TFile rootfile(outfilename.c_str(), "recreate");
	
	  // Read the Endianness and set the byteswap
	  unsigned int endianness;
	  infile.read( (char*)&endianness, sizeof(endianness) );        
	  // Set the byteswap boolean if needed.
	  bool byteswap = false;
	  if (endianness != 0x01020304 and endianness != 0) byteswap = true;
	  
	  // read and store the xml string 
	  unsigned int xml_settings_string_length;
    infile.readval(xml_settings_string_length); 
    char* xml_settings_string = new char[xml_settings_string_length+1];
    infile.readstr(xml_settings_string, xml_settings_string_length);
    xml_settings_string[xml_settings_string_length] = '\0';
	  TObjString xml( (const char*)xml_settings_string );
	  xml.Write("xml");
	  
	  // read and write the header block
	  string tree_name = "header";
	  rqReader *header = new rqReader(infile, tree_name, byteswap);
	  delete header;
	  
	  // read and write the rq1 block 
	  tree_name = "events";
	  rqReader *events = new rqReader(infile, tree_name, byteswap);
	  delete events;
	  
	  // read and write the live time block 
	  tree_name = "livetime";
	  rqReader *livetime = new rqReader(infile, tree_name, byteswap);
	  delete livetime;
	  
	  // close the files
	  infile.close();
	  rootfile.Close();    
	}
	return 0;
}
