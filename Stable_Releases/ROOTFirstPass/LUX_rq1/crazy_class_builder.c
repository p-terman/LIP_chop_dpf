// So can we actually build a ROOT class with no a priori knowledge
// of the structure of the memory structure.
// This reads from the binary file and sets up the tree structure.  There's
// another code that does the filling of the tree.
// This is a totally weird way to go about it, but we'll see if it works.
//
// K.Clark - 07/04/10
// ken.clark@case.edu

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <ctype.h>
#include <vector>
#include <string>
#include <sstream>

#include <TROOT.h>
#include <TObject.h>
#include <TTree.h>
#include <TFile.h>
#include <TGraph.h>
#include <TCanvas.h>

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

Int_t stupid_method(const char *type){

  if (!strcmp(type,"char") || !strcmp(type,"int8") || !strcmp(type,"uint8")) return 1;
  else if (!strcmp(type,"int16") || !strcmp(type,"uint16")) return 2;
  else if (!strcmp(type,"int32") || !strcmp(type,"uint32") || !strcmp(type,"float") || !strcmp(type,"single")) return 4;
  else if (!strcmp(type,"int64") || !strcmp(type,"uint64") || !strcmp(type,"double")) return 8;
  else return 16;

}

const Char_t* stupider_method(const char *type){

  if (!strcmp(type,"char")) return "string";
  else if (!strcmp(type,"uint32")) return "UInt_t";
  else if (!strcmp(type,"double")) return "Double_t";
  else return "Screwed";
}

int set_up_header_files(const char *file) {
  
  ifstream infile;
  infile.open(file,ios::in | ios::binary);
  if (infile.is_open()){
    Int_t lng_sz=sizeof(UInt_t);
    UShort_t temp_headerlength;
    Int_t headerlength;
    Int_t headerlength_swap;
    Char_t *temp_entry;
    Int_t num_vars;
    Int_t semicolon_ct=0;
    Int_t num_lines;
    Char_t *header;
    Int_t header_bytes=0;
    Bool_t swap_needed=0;

    Char_t **var_names;
    Char_t **var_types;
    Int_t *var_count;
    Int_t *var_indices;
    Char_t **pulse_var_names;
    Char_t **pulse_var_types;
    Int_t *pulse_var_indices;
    Int_t *pulse_var_count;
    Char_t **event_var_names;
    Char_t **event_var_types;
    Int_t *event_var_count;
    Char_t **header_var_names;
    Char_t **header_var_types;
    Int_t *header_var_count;
    Char_t **channel_var_names;
    Char_t **channel_var_types;
    Int_t *channel_var_indices;
    Int_t *channel_var_count;

    Int_t *var_type_size;
    Char_t **temp_char;
    UInt_t *temp_int;
    Double_t *temp_dbl;
    UInt_t **temp_int_array;
    Double_t **temp_dbl_array;
    string temp_string;
    Char_t *strchr_val;
    Int_t num_len;
    Char_t *sub_str;
    Int_t event_var_ct=0;
    Int_t pulse_var_ct=0;
    Int_t header_var_ct=0;
    Int_t channel_var_ct=0;
    Int_t startpos=0;
    Int_t header_char_ct=0;
    Int_t header_int_ct=0;
    Int_t header_dbl_ct=0;
    Int_t event_char_ct=0;
    Int_t event_int_ct=0;
    Int_t event_dbl_ct=0;
    Int_t pulse_char_ct=0;
    Int_t pulse_int_ct=0;
    Int_t pulse_dbl_ct=0;

    Int_t read_next_xml=0;
    Double_t infile_vers_num;
    string read_junk, read_junk_2;
    Double_t header_vers_num;
    //Need to adjust the reader to read the xml header
    //and then decide if the classes need to be rebuilt.
    
    infile.seekg(0,ios::beg);
    UInt_t endianness;
    infile.read((char*)&endianness,sizeof(endianness));
    if (endianness != 0x01020304) swap_needed=1;

    //Lines below should be commented out for files with no xml headers
    //since I couldn't think of a really clever way to distinguish those with
    //from those without headers. - KJC 06/05/10
    
    UInt_t settings_chars;
    infile.read((char*)&settings_chars,sizeof(settings_chars));
    //need to read in settings_chars characters, can't read in line by line
    //or else should the characters be counted as they are read in?
    if (swap_needed) full_endian_swap(settings_chars);
    printf("There are %i characters in the header.\n",settings_chars);
    //printf("There are %i characters in the xml settings files\n",settings_chars);
    Char_t temp_xml[settings_chars];
    printf("Before reading the xml header\n");
    infile.read((char*)&temp_xml,settings_chars*sizeof(Char_t));
    printf("Before deciding about xml header\n");

    // Begin comment block to remove strtok method
    /*if (strcmp(&temp_xml[0],"<")){
      printf("Found an xml header, I'll read that first.\n");
      Char_t *xml_piece;
      //xml_piece = strtok(temp_xml,"<>");
      if (swap_needed) full_endian_swap(xml_piece);
      //xml_piece = strtok(NULL,"<>");
      //printf("%s\n",xml_piece);
      while (xml_piece != NULL){
	if (!strcmp(xml_piece,"analysis_version")) {read_next_xml=1;}
	xml_piece = strtok(NULL,"<>");
	if (swap_needed) full_endian_swap(xml_piece);
	if (read_next_xml) {
	  infile_vers_num = strtod(xml_piece,NULL);
	  //printf("Analysis version: %.2f\n",infile_vers_num);
	  read_next_xml=0;
	}
      }
    */ // End comment block to get rid of strtok method PHP - 2010/10/15

    if(strcmp(&temp_xml[0],"<")){
      printf("Found an xml header, will read that first. \n");
      Char_t* analysis_version_parse_str1 = "<analysis_version>";
      Char_t* analysis_version_parse_str2 = "</analysis_version>";
      size_t pos2, pos3;
      pos2=(string(temp_xml)).find(analysis_version_parse_str1);
      pos3=(string(temp_xml)).find(analysis_version_parse_str2);
      string analysis_version_result = (string(temp_xml)).substr(pos2+(string(analysis_version_parse_str1)).length(),pos3-pos2-(string(analysis_version_parse_str1)).length());
      istringstream analysis_version_buffer;
      analysis_version_buffer.str(analysis_version_result);
      analysis_version_buffer >> infile_vers_num;
      cout << "Analysis version: " << infile_vers_num << endl;
     
 
      ifstream event_head;
      event_head.open("LUX_rq1_event.h");
      if(event_head){
	getline(event_head,read_junk);
	getline(event_head,read_junk_2);
	//cout << "This is read_junk " << read_junk << endl;
	//cout << "This is read_junk_2 " << read_junk_2 << endl;
	if(!read_junk.empty()){
	  header_vers_num = strtod((read_junk_2.substr(read_junk_2.find_first_not_of("/ "))).c_str(),NULL);
	}
	else{header_vers_num = 0;}
	//printf("Header version is %.3f\n",header_vers_num);
	//printf("Data file version is %.3f\n",infile_vers_num);
	if (header_vers_num==infile_vers_num){
	  printf("The version numbers agree, I don't have to rebuild the classes.\n");
	  exit(1);
	}
	else printf("This data needs new classes, I'll get to work rebuilding them.\n");
      }
    }
    else printf("No xml header in this file.\n");
    
    // Commenting out of the lines (for files without xml files)
    // should end here.

    ofstream rq1_header_h("LUX_rq1_header.h");
    ofstream rq1_header_cxx("LUX_rq1_header.cxx");
    ofstream rq1_event_h("LUX_rq1_event.h");
    ofstream rq1_event_cxx("LUX_rq1_event.cxx");
    ofstream rq1_pulse_h("LUX_rq1_pulse.h");
    ofstream rq1_pulse_cxx("LUX_rq1_pulse.cxx");
    ofstream rq1_channel_h("LUX_rq1_channel.h");
    ofstream rq1_channel_cxx("LUX_rq1_channel.cxx");

    rq1_header_h << "#ifndef LUX_RQ1_HEADER_H\n";
    rq1_header_h << "#define LUX_RQ1_HEADER_H\n\n";
    rq1_header_h << "#include \"TObject.h\"\n\n";
    rq1_header_h << "#include <iostream>\n";
    rq1_header_h << "#include <string>\n";
    rq1_header_h << "using namespace std;\n\n";
    rq1_header_h << "class LUX_rq1_event;\n";
    rq1_header_h << "class LUX_rq1_pulse;\n";
    rq1_header_h << "class LUX_rq1_channel;\n";

    rq1_header_cxx << "#include <iostream>\n";
    rq1_header_cxx << "#include <string>\n";
    rq1_header_cxx << "#include \"LUX_rq1_pulse.h\"\n";
    rq1_header_cxx << "#include \"LUX_rq1_header.h\"\n";
    rq1_header_cxx << "#include \"LUX_rq1_channel.h\"\n";
    rq1_header_cxx << "#include \"LUX_rq1_event.h\"\n\n";
    rq1_header_cxx << "using namespace std;\n\n";
    rq1_header_cxx << "ClassImp(LUX_rq1_header)\n\n";

    rq1_event_h << "// Analsysis version this class was created with\n";
    rq1_event_h << showpoint << "// " << infile_vers_num << "\n";
    rq1_event_h << "// Don't change the above lines!!\n";
    rq1_event_h << "#ifndef LUX_RQ1_EVENT_H\n";
    rq1_event_h << "#define LUX_RQ1_EVENT_H\n";
    rq1_event_h << "#include \"TObject.h\"\n";
    rq1_event_h << "#include \"TClonesArray.h\"\n";
    rq1_event_h << "#include <vector>\n";
    rq1_event_h << "class LUX_rq1_pulse;\n";
    rq1_event_h << "class LUX_rq1_header;\n";
    rq1_event_h << "class LUX_rq1_channel;\n";

    rq1_event_cxx << "#include <iostream>\n";
    rq1_event_cxx << "#include \"LUX_rq1_pulse.h\"\n";
    rq1_event_cxx << "#include \"LUX_rq1_header.h\"\n";
    rq1_event_cxx << "#include \"LUX_rq1_event.h\"\n";
    rq1_event_cxx << "#include \"LUX_rq1_channel.h\"\n";
    rq1_event_cxx << "using namespace std;\n\n";
    rq1_event_cxx << "ClassImp(LUX_rq1_event)\n\n";
 
    rq1_pulse_h << "#ifndef LUX_RQ1_PULSE_H\n";
    rq1_pulse_h << "#define LUX_RQ1_PULSE_H\n";
    rq1_pulse_h << "#include \"TObject.h\"\n";
    rq1_pulse_h << "\n";
    rq1_pulse_h << "#include <iostream>\n";
    rq1_pulse_h << "using namespace std;\n";
    rq1_pulse_h << "\n";
    rq1_pulse_h << "class LUX_rq1_event;\n";
    rq1_pulse_h << "class LUX_rq1_header;\n";
    rq1_pulse_h << "class LUX_rq1_channel;\n\n";
    rq1_pulse_h << "class LUX_rq1_pulse : public TObject{\n";
    rq1_pulse_h << " public:\n";
    rq1_pulse_h << "  LUX_rq1_pulse();\n";
    
    rq1_pulse_cxx << "#include <iostream>\n";
    rq1_pulse_cxx << "#include \"LUX_rq1_event.h\"\n";
    rq1_pulse_cxx << "#include \"LUX_rq1_pulse.h\"\n";
    rq1_pulse_cxx << "#include \"LUX_rq1_header.h\"\n";
    rq1_pulse_cxx << "#include \"LUX_rq1_channel.h\"\n";
    rq1_pulse_cxx << "using namespace std;\n\n";
    rq1_pulse_cxx << "ClassImp(LUX_rq1_pulse);\n\n";
    
    rq1_channel_h << "#ifndef LUX_RQ1_CHANNEL_H\n";
    rq1_channel_h << "#define LUX_RQ1_CHANNEL_H\n";
    rq1_channel_h << "#include \"TObject.h\"\n";
    rq1_channel_h << "\n";
    rq1_channel_h << "#include <iostream>\n";
    rq1_channel_h << "using namespace std;\n";
    rq1_channel_h << "\n";
    rq1_channel_h << "class LUX_rq1_event;\n";
    rq1_channel_h << "class LUX_rq1_header;\n";
    rq1_channel_h << "class LUX_rq1_pulse;\n";
    rq1_channel_h << "class LUX_rq1_channel : public TObject{\n";
    rq1_channel_h << " public:\n";
    rq1_channel_h << "  LUX_rq1_channel();\n";
    rq1_channel_h << "  LUX_rq1_channel(LUX_rq1_channel *lrq1c);\n";
    rq1_channel_h << "  ~LUX_rq1_channel();\n";
    
    rq1_channel_cxx << "#include <iostream>\n";
    rq1_channel_cxx << "#include \"LUX_rq1_event.h\"\n";
    rq1_channel_cxx << "#include \"LUX_rq1_pulse.h\"\n";
    rq1_channel_cxx << "#include \"LUX_rq1_header.h\"\n";
    rq1_channel_cxx << "#include \"LUX_rq1_channel.h\"\n";
    rq1_channel_cxx << "using namespace std;\n\n";
    rq1_channel_cxx << "ClassImp(LUX_rq1_channel);\n\n";

    //Now it's time to read out the blocks of data.  There are two.
    //The first is the overall header and sets up the reading of the second.

    for (Int_t block=0;block<2;block++){
      infile.read((char*)&temp_headerlength,sizeof(UShort_t));
      if (swap_needed) full_endian_swap(temp_headerlength);
      headerlength = (Int_t)temp_headerlength;
      //printf("The header is %i characters long\n",headerlength);
      header = new Char_t[headerlength];
      infile.read((Char_t*)header,headerlength*sizeof(Char_t));
      //if (swap_needed) full_endian_swap(header);
      semicolon_ct=0;
      for (Int_t i=0;i<headerlength;i++) if (header[i]==';'){semicolon_ct++;}
      num_vars = semicolon_ct/3;
      
      var_names = new Char_t*[num_vars];
      var_types = new Char_t*[num_vars];
      var_count = new Int_t[num_vars];
      var_indices = new Int_t[num_vars];
      var_type_size = new Int_t[num_vars];

      if (block==1){
	channel_var_names = new Char_t*[num_vars];
	channel_var_types = new Char_t*[num_vars];
	channel_var_count = new Int_t[num_vars];
	channel_var_indices = new Int_t[num_vars];
	pulse_var_names = new Char_t*[num_vars];
	pulse_var_types = new Char_t*[num_vars];
	pulse_var_count = new Int_t[num_vars];
      }
      else if (block==0){
	header_var_names = new Char_t*[num_vars];
	header_var_types = new Char_t*[num_vars];
	header_var_count = new Int_t[num_vars];
      }
      
      for (Int_t i=0;i<num_vars;i++){
	if (i==0) temp_entry = strtok(header,";");
	else temp_entry = strtok(NULL,";");
	//printf("%s\n",temp_entry);
	//if (swap_needed) full_endian_swap(temp_entry); -> strings shouldn't need swaps.
	var_names[i]=temp_entry;
	temp_entry = strtok(NULL,";");
	//if (swap_needed) full_endian_swap(temp_entry);
	var_types[i]=temp_entry;
	var_type_size[i] = stupid_method(temp_entry);
	temp_entry = strtok(NULL,";");
	//Don't know if this swap actually goes here or not...
	//if (swap_needed) full_endian_swap(temp_entry);

	if ((strchr_val = strchr(temp_entry,','))!=NULL){
	  sub_str = new Char_t;
	  num_len = strchr_val - temp_entry;
	  memcpy(sub_str,temp_entry,num_len);\
	  //if (swap_needed) full_endian_swap(sub_str);
	  var_count[i]=atoi(sub_str);
	  memcpy(sub_str,temp_entry+2,num_len+2);
	  //if (swap_needed) full_endian_swap(sub_str);
	  var_indices[i]=atoi(sub_str);
	  delete[] sub_str;
	  if (block==1){
	    channel_var_names[channel_var_ct]=var_names[i];
	    channel_var_types[channel_var_ct]=var_types[i];
	    channel_var_indices[channel_var_ct]=var_indices[i];
	    channel_var_count[channel_var_ct]=var_count[i];
	    //printf("channel var: %s\n",channel_var_names[pulse_var_
	    channel_var_ct++;
	  }
	}
	else{
	  //if (swap_needed) full_endian_swap(temp_entry);
	  var_count[i]=atoi(temp_entry);
	  var_indices[i]=1;
	  if (block==0){
	    header_var_names[header_var_ct]=var_names[i];
	    header_var_types[header_var_ct]=var_types[i];
	    header_bytes += Int_t(var_count[i]*var_type_size[i]);
	    header_var_count[header_var_ct] = var_count[i];
	    header_var_ct++;
	  }
	  if (block==1){
	    pulse_var_names[pulse_var_ct]=var_names[i];
	    pulse_var_types[pulse_var_ct]=var_types[i];
	    pulse_var_count[pulse_var_ct]=var_count[i];
	    pulse_var_ct++;
	  }
	}
      }
      if (block==0){
	rq1_header_h << "class LUX_rq1_header : public TObject{\n\n";
	rq1_header_h << "public:\n";
	rq1_header_h << "  LUX_rq1_header();\n";
	rq1_header_h << "  LUX_rq1_header(LUX_rq1_header *lrq1h);\n";
	rq1_header_h << "  ~LUX_rq1_header();\n";
	rq1_header_h << "  LUX_rq1_header(";
	for (Int_t i=0;i<header_var_ct-1;i++) rq1_header_h << stupider_method(header_var_types[i]) << ",";
	rq1_header_h << stupider_method(header_var_types[header_var_ct-1]) << ");\n\n";
	for (Int_t i=0;i<header_var_ct;i++) 
	  if (!strcmp(header_var_types[i],"char")) rq1_header_h << "  void Set_" << header_var_names[i] << "(string x){" << header_var_names[i] << " = x;};\n";	
	else rq1_header_h << "  void Set_" << header_var_names[i] << "(" << stupider_method(header_var_types[i]) << " x){" << header_var_names[i] << " = x;};\n";
	rq1_header_h << "\n";
	for (Int_t i=0;i<header_var_ct;i++) {
	  if (!strcmp(header_var_types[i],"char")) rq1_header_h << "  string Get_" << header_var_names[i] << "(){return " <<  header_var_names[i] << ";};\n";
	  else rq1_header_h << "  " << stupider_method(header_var_types[i]) << " Get_" << header_var_names[i] << "(){return " <<  header_var_names[i] << ";};\n";
	}
	rq1_header_h << "\n";
	rq1_header_h << "  void Fill(";
	for (Int_t i=0;i<header_var_ct-1;i++) rq1_header_h << stupider_method(header_var_types[i]) << ",";
	rq1_header_h << stupider_method(header_var_types[header_var_ct-1]) << ");\n";
	rq1_header_h << "  void Fill(string str);\n";
	rq1_header_h << "  void Fill(Char_t** a, UInt_t* b, Double_t* c);\n";
	rq1_header_h << "  void Clear(const Option_t* option=\"\");\n\n";
	rq1_header_h << "private:\n";
	for (Int_t i=0;i<header_var_ct;i++){
	  if (!strcmp(header_var_types[i],"char")) {
	    rq1_header_h << "  " << "string " << header_var_names[i] << ";\n";
	  }
	  else {
	    rq1_header_h << "  " << stupider_method(header_var_types[i]) << " " << header_var_names[i] << ";\n";
	  }
	}
	rq1_header_h << "\n";
	rq1_header_h << "ClassDef(LUX_rq1_header,1);\n";
	rq1_header_h << "};\n";
	rq1_header_h << "#endif\n";	
	rq1_header_h.close();
	
	rq1_header_cxx << "LUX_rq1_header::LUX_rq1_header(){\n";
	for (Int_t i=0;i<header_var_ct;i++) {
	  if (!strcmp(header_var_types[i],"char")) rq1_header_cxx << "  " << header_var_names[i] << " = \"\";\n";
	  else rq1_header_cxx << "  " << header_var_names[i] << " = 0;\n";
	}	
	rq1_header_cxx << "}\n\n";
	Char_t name1='a';
	rq1_header_cxx << "LUX_rq1_header::LUX_rq1_header(";
	for (Int_t i=0;i<header_var_ct-1;i++) rq1_header_cxx << stupider_method(header_var_types[i]) << " " << Char_t(name1+i) <<  ", ";
	rq1_header_cxx << stupider_method(header_var_types[header_var_ct-1]) << " " << Char_t(name1+header_var_ct-1) << "){\n";
	for (Int_t i=0;i<header_var_ct;i++) rq1_header_cxx << "  " << header_var_names[i] << " = " << Char_t(name1+i) << ";\n";
	rq1_header_cxx << "}\n\n";
	rq1_header_cxx << "LUX_rq1_header::LUX_rq1_header(LUX_rq1_header *lrq1h){\n";
	for (Int_t i=0;i<header_var_ct;i++) rq1_header_cxx << "  " << header_var_names[i] << " = lrq1h->Get_" << header_var_names[i] << "();\n";
	rq1_header_cxx << "}\n\n";
	rq1_header_cxx << "LUX_rq1_header::~LUX_rq1_header(){}\n\n";
	rq1_header_cxx << "void LUX_rq1_header::Fill(";
	for (Int_t i=0;i<header_var_ct-1;i++) rq1_header_cxx << stupider_method(header_var_types[i]) << " " << Char_t(name1+i) <<  ", ";
	rq1_header_cxx << stupider_method(header_var_types[header_var_ct-1]) << " " << Char_t(name1+header_var_ct-1) << "){\n";
	for (Int_t i=0;i<header_var_ct;i++) rq1_header_cxx << "  " << header_var_names[i] << " = " << Char_t(name1+i) << ";\n";
	rq1_header_cxx << "}\n\n";
	
	rq1_header_cxx << "void LUX_rq1_header::Fill(string str){\n";
	for (Int_t i=0;i<header_var_ct;i++) {
	  if (!strcmp(header_var_types[i],"char")) {
	    rq1_header_cxx << "  " << header_var_names[i] << " = str.substr(" << startpos << "," << startpos+(header_var_count[i]*stupid_method(header_var_types[i])) << ").c_str();\n";
	  }
	  else if (!strcmp(header_var_types[i],"uint32"))  {
	    rq1_header_cxx << "  " << header_var_names[i] << " = atoi(str.substr(" << startpos << "," << startpos+(header_var_count[i]*stupid_method(header_var_types[i])-1) << ").c_str());\n";
	  }
	  else if (!strcmp(header_var_types[i],"double"))  {
	    rq1_header_cxx << "  " << header_var_names[i] << " = strtod(str.substr(" << startpos << "," << startpos+(header_var_count[i]*stupid_method(header_var_types[i])) << ").c_str());\n";
	  }
	  startpos = startpos+(header_var_count[i]*stupid_method(header_var_types[i]));
	}
	rq1_header_cxx << "}\n\n";
	
	rq1_header_cxx << "void LUX_rq1_header::Fill(Char_t** a, UInt_t* b, Double_t* c){\n";
	for (Int_t i=0;i<header_var_ct;i++){
	  if (!strcmp(header_var_types[i],"char")) {
	    rq1_header_cxx << "  " << header_var_names[i] << " = a[" << header_char_ct << "];\n";
	    header_char_ct++;
	  }
	  else if (!strcmp(header_var_types[i],"uint32")) {
	    rq1_header_cxx << "  " << header_var_names[i] << " = b[" << header_int_ct << "];\n";
	    header_int_ct++;
	  }
	  else if (!strcmp(header_var_types[i],"double")) {
	    rq1_header_cxx << "  " << header_var_names[i] << " = c[" << header_dbl_ct << "];\n";
	    header_dbl_ct++;
	  }
	}
	rq1_header_cxx << "}\n\n";
	rq1_header_cxx << "void LUX_rq1_header::Clear(const Option_t* option){\n";
	for (Int_t i=0;i<header_var_ct;i++) {
	  if (!strcmp(header_var_types[i],"char")) rq1_header_cxx << "  " << header_var_names[i] << " = \"\";\n";
	  else rq1_header_cxx << "  " << header_var_names[i] << " = 0;\n";
	}
	rq1_header_cxx << "}\n";
	rq1_header_cxx.close();
	infile.read((char*)&num_lines,sizeof(Int_t));
	if (swap_needed) full_endian_swap(num_lines);
	infile.ignore(header_bytes*num_lines);
      }
      //printf("Got to the building of the event, pulse and channel files\n");
      //here build the event, pulse and channel header and cxx files
      if (block==1){
	rq1_channel_h << "  LUX_rq1_channel(";
	for (Int_t i=0;i<channel_var_ct-1;i++) rq1_channel_h << stupider_method(channel_var_types[i]) << ", ";
	rq1_channel_h << stupider_method(channel_var_types[channel_var_ct-1]) << ");\n\n";
	for (Int_t i=0;i<channel_var_ct;i++) rq1_channel_h << "  void Set_" << channel_var_names[i] << "(" << stupider_method(channel_var_types[i]) << " x){" << channel_var_names[i] << "=x;};\n";
	rq1_channel_h << "\n";
	for (Int_t i=0;i<channel_var_ct;i++) rq1_channel_h << "  " << stupider_method(channel_var_types[i]) << " Get_" << channel_var_names[i] << "(){return " << channel_var_names[i] << ";};\n";
	rq1_channel_h << "\n  void Fill(";
	for (Int_t i=0;i<channel_var_ct-1;i++) rq1_channel_h << stupider_method(channel_var_types[i]) << ", ";
	rq1_channel_h << stupider_method(channel_var_types[channel_var_ct-1]) << ");\n";
	
	rq1_channel_h << "  void Fill(Char_t** a, UInt_t* b, Double_t* c);\n";
	rq1_channel_h << "  void Clear(const Option_t* option=\"\");\n\n";
	rq1_channel_h << " private:\n";
	for (Int_t i=0;i<channel_var_ct;i++) rq1_channel_h << "  " << stupider_method(channel_var_types[i]) << " " << channel_var_names[i] << ";\n";
	rq1_channel_h << "  ClassDef(LUX_rq1_channel,1);\n";
	rq1_channel_h << "};\n";
	rq1_channel_h << "#endif\n";
	rq1_channel_h.close();

	rq1_channel_cxx << "LUX_rq1_channel::LUX_rq1_channel(){\n";
	for (Int_t i=0;i<channel_var_ct;i++) rq1_channel_cxx << "  " << channel_var_names[i] << " = 0;\n";
	rq1_channel_cxx << "}\n\n";
	rq1_channel_cxx << "LUX_rq1_channel::LUX_rq1_channel(";
	Char_t name='a';
	for (Int_t i=0;i<channel_var_ct-1;i++){
	rq1_channel_cxx << stupider_method(channel_var_types[i]) << " " << Char_t(name+i) << ", ";
	
	}
	
	
	
	rq1_channel_cxx << stupider_method(channel_var_types[channel_var_ct-1]) << " " << Char_t(name+channel_var_ct-1) << "){\n";
	for (Int_t i=0;i<channel_var_ct;i++) rq1_channel_cxx << "  " << channel_var_names[i] << " = " << Char_t(name+i) << ";\n";
	rq1_channel_cxx << "}\n\n";
	rq1_channel_cxx << "LUX_rq1_channel::LUX_rq1_channel(LUX_rq1_channel *lrq1c){\n";
	for (Int_t i=0;i<channel_var_ct;i++) rq1_channel_cxx << "  " << channel_var_names[i] << " = lrq1c->Get_" << channel_var_names[i] << "();\n";
	rq1_channel_cxx << "}\n\n";
	rq1_channel_cxx << "LUX_rq1_channel::~LUX_rq1_channel(){}\n\n";
	rq1_channel_cxx << "void LUX_rq1_channel::Fill(";
	for (Int_t i=0;i<channel_var_ct-1;i++) rq1_channel_cxx << stupider_method(channel_var_types[i]) << " " << Char_t(name+i) << ", ";
	rq1_channel_cxx << stupider_method(channel_var_types[channel_var_ct-1]) << " " << Char_t(name+channel_var_ct-1) << "){\n";
	for (Int_t i=0;i<channel_var_ct;i++) rq1_channel_cxx << "  " << channel_var_names[i] << " = " << Char_t(name+i) << ";\n";
	rq1_channel_cxx << "}\n\n";
	rq1_channel_cxx << "void LUX_rq1_channel::Fill(Char_t** a, UInt_t* b, Double_t* c){\n";
	for (Int_t i=0;i<channel_var_ct;i++){
	  if (!strcmp(channel_var_types[i],"char")) rq1_channel_cxx << "  " << channel_var_names[i] << " = a[" << i << "];\n";
	  else if (!strcmp(channel_var_types[i],"uint32")) rq1_channel_cxx << "  " << channel_var_names[i] << " = b[" << i << "];\n";
	  else if (!strcmp(channel_var_types[i],"double")) rq1_channel_cxx << "  " << channel_var_names[i] << " = c[" << i << "];\n";
	}
	rq1_channel_cxx << "}\n\n";
	
	rq1_channel_cxx << "void LUX_rq1_channel::Clear(const Option_t* option){\n";
	for (Int_t i=0;i<channel_var_ct;i++) rq1_channel_cxx << "  " << channel_var_names[i] << " = 0;\n";
	rq1_channel_cxx << "}\n\n";
	rq1_channel_cxx.close();

 	rq1_event_h << "class LUX_rq1_event : public TObject{\n\n";
	rq1_event_h << "public :\n\n";
	rq1_event_h << "  LUX_rq1_event();\n";
	rq1_event_h << "  LUX_rq1_event(LUX_rq1_event *lrq1e);\n";
	rq1_event_h << "  LUX_rq1_event(UInt_t a, UInt_t b);\n";
	rq1_event_h << "  ~LUX_rq1_event();\n\n";
	rq1_event_h << "  void Set_num_pulses(UInt_t x){num_pulses=x;};\n\n";
	rq1_event_h << "  UInt_t Get_num_pulses(){return num_pulses;};\n";
	rq1_event_h << "  LUX_rq1_pulse* Get_pulse(UInt_t num);\n\n";
	rq1_event_h << "  void AddPulse(LUX_rq1_pulse *lrq1p, Int_t num);\n";
	rq1_event_h << "  void Clear(const Option_t* option=\"\");\n\n";
	rq1_event_h << "  void Set_num_channels(UInt_t x){num_channels=x;};\n\n";
	rq1_event_h << "  UInt_t Get_num_channels(){return num_channels;};\n";
	rq1_event_h << "  LUX_rq1_channel* Get_channel(UInt_t n_pulse, UInt_t n_channel);\n\n";
	rq1_event_h << "  void AddChannel(LUX_rq1_channel *lrq1c, Int_t num);\n";
	rq1_event_h << "private :\n\n";
	rq1_event_h << "  TClonesArray *pulse;\n";
	rq1_event_h << "  UInt_t num_pulses;\n\n";
	rq1_event_h << "  TClonesArray *channel;\n";
	rq1_event_h << "  UInt_t num_channels;\n\n";
	rq1_event_h << "  ClassDef(LUX_rq1_event,1);\n";
	rq1_event_h << "};\n\n";
	rq1_event_h << "#endif\n";
	rq1_event_h.close();
	
	rq1_event_cxx << "LUX_rq1_event::LUX_rq1_event(){\n";
	rq1_event_cxx << "  num_pulses = 0;\n";
	rq1_event_cxx << "  pulse = new TClonesArray(\"LUX_rq1_pulse\");\n";
	rq1_event_cxx << "  pulse->BypassStreamer();\n";
	rq1_event_cxx << "  num_channels = 0;\n";
	rq1_event_cxx << "  channel = new TClonesArray(\"LUX_rq1_channel\");\n";
	rq1_event_cxx << "  channel->BypassStreamer();\n";
	rq1_event_cxx << "}\n\n";
	rq1_event_cxx << "void LUX_rq1_event::Clear(const Option_t* option){\n";
	rq1_event_cxx << "  num_pulses=0;\n";
	rq1_event_cxx << "  num_channels=0;\n";
	rq1_event_cxx << "}\n\n";
	rq1_event_cxx << "LUX_rq1_event::LUX_rq1_event(UInt_t a, UInt_t b){\n";
	rq1_event_cxx << "  num_pulses=a;\n";
	rq1_event_cxx << "  pulse = new TClonesArray(\"LUX_rq1_pulse\");\n";
	rq1_event_cxx << "  pulse->BypassStreamer();\n";
	rq1_event_cxx << "  num_channels=b;\n";
	rq1_event_cxx << "  channel = new TClonesArray(\"LUX_rq1_channel;\");\n";
	rq1_event_cxx << "  channel->BypassStreamer();\n";
	rq1_event_cxx << "}\n\n";
	rq1_event_cxx << "LUX_rq1_event::LUX_rq1_event(LUX_rq1_event *lrq1e){\n";
	rq1_event_cxx << "  num_pulses = lrq1e->Get_num_pulses();\n";
	rq1_event_cxx << "  num_channels = lrq1e->Get_num_channels();\n";
	rq1_event_cxx << "}\n\n";
	rq1_event_cxx << "LUX_rq1_event::~LUX_rq1_event(){delete pulse;delete channel;}\n\n";
	rq1_event_cxx << "void LUX_rq1_event::AddPulse(LUX_rq1_pulse *lrq1p, Int_t index){\n";
	rq1_event_cxx << "  TClonesArray &pls = *pulse;\n";
	rq1_event_cxx << "  pls[index] = new(pls[index]) LUX_rq1_pulse(lrq1p);\n";
	rq1_event_cxx << "}\n\n";
	rq1_event_cxx << "LUX_rq1_pulse* LUX_rq1_event::Get_pulse(UInt_t num){\n";
	rq1_event_cxx << "  if (num<num_pulses) return (LUX_rq1_pulse*)pulse->At(num);\n";
	rq1_event_cxx << "  else {\n";
	rq1_event_cxx << "    cerr << \"Pulse index out of range, returning 0\"<<endl;\n";
	rq1_event_cxx << "    return 0;\n";
	rq1_event_cxx << "  }\n";
	rq1_event_cxx << "}\n\n";
	rq1_event_cxx << "void LUX_rq1_event::AddChannel(LUX_rq1_channel *lrq1c, Int_t index){\n";
	rq1_event_cxx << "  TClonesArray &chan = *channel;\n";
	rq1_event_cxx << "  chan[index] = new(chan[index]) LUX_rq1_channel(lrq1c);\n";
	rq1_event_cxx << "}\n\n";
	rq1_event_cxx << "LUX_rq1_channel* LUX_rq1_event::Get_channel(UInt_t n_pulse, UInt_t n_channel){\n";
	rq1_event_cxx << "  UInt_t chan_number = n_pulse*num_channels+n_channel;\n";
	rq1_event_cxx << "  UInt_t total_chans = num_channels*num_pulses;\n";
	rq1_event_cxx << "  if (chan_number<total_chans) return (LUX_rq1_channel*)channel->At(chan_number);\n";
	rq1_event_cxx << "  else {\n";
	rq1_event_cxx << "    cerr << \"Channel index out of range, returning 0\"<<endl;\n";
	rq1_event_cxx << "    return 0;\n";
	rq1_event_cxx << "  }\n";
	rq1_event_cxx << "}\n";
	rq1_event_cxx.close();

	rq1_pulse_h << "  LUX_rq1_pulse(";
	for (Int_t i=0;i<pulse_var_ct-1;i++) rq1_pulse_h << stupider_method(pulse_var_types[i]) << ", ";
	rq1_pulse_h << stupider_method(pulse_var_types[pulse_var_ct-1]) << ");\n";
	rq1_pulse_h << "  LUX_rq1_pulse(LUX_rq1_pulse *lrq1p);\n";
	rq1_pulse_h << "  ~LUX_rq1_pulse();\n\n";
	for (Int_t i=0;i<pulse_var_ct;i++) rq1_pulse_h << "  void Set_" << pulse_var_names[i] << "(" << stupider_method(pulse_var_types[i]) << " x){" << pulse_var_names[i] << " = x;};\n";
	rq1_pulse_h << "\n";
	for (Int_t i=0;i<pulse_var_ct;i++) rq1_pulse_h << "  " << stupider_method(pulse_var_types[i]) << " Get_" << pulse_var_names[i] << "(){return " << pulse_var_names[i] << ";};\n";
	rq1_pulse_h << "\n";
	rq1_pulse_h << "  void Fill(";
	for (Int_t i=0;i<pulse_var_ct-1;i++) rq1_pulse_h << stupider_method(pulse_var_types[i]) << ", ";
	rq1_pulse_h << stupider_method(pulse_var_types[pulse_var_ct-1]) << ");\n";
	rq1_pulse_h << "  void Fill(Char_t** a, UInt_t* b, Double_t* c);\n";
	rq1_pulse_h << "  void Clear(const Option_t* option=\"\");\n\n";
	rq1_pulse_h << "private : \n\n";
	for (Int_t i=0;i<pulse_var_ct;i++) rq1_pulse_h << "  " << stupider_method(pulse_var_types[i]) << " " << pulse_var_names[i] << ";\n";
	rq1_pulse_h << "  ClassDef(LUX_rq1_pulse,2);\n";
	rq1_pulse_h << "};\n\n";
	rq1_pulse_h << "#endif\n";
	rq1_pulse_h.close();

	rq1_pulse_cxx << "LUX_rq1_pulse::LUX_rq1_pulse(){\n";
	for (Int_t i=0;i<pulse_var_ct;i++) rq1_pulse_cxx << "  " << pulse_var_names[i] << " = 0;\n";
	rq1_pulse_cxx << "}\n\n";
	rq1_pulse_cxx << "void LUX_rq1_pulse::Clear(const Option_t* option){\n";
	for (Int_t i=0;i<pulse_var_ct;i++) rq1_pulse_cxx << "  " << pulse_var_names[i] << " = 0;\n";
	rq1_pulse_cxx << "}\n";
	rq1_pulse_cxx << "LUX_rq1_pulse::LUX_rq1_pulse(";
	for (Int_t i=0;i<pulse_var_ct-1;i++){ 
	  rq1_pulse_cxx << stupider_method(pulse_var_types[i]) << " " << name << i << ", ";
	}
	  rq1_pulse_cxx << stupider_method(pulse_var_types[pulse_var_ct-1]) << " " << name << pulse_var_ct-1 << "){\n";
	for (Int_t i=0;i<pulse_var_ct;i++){
	  rq1_pulse_cxx << "  " << pulse_var_names[i] << " = " << name << i << ";\n";
	}
	rq1_pulse_cxx << "}\n\n";
	rq1_pulse_cxx << "LUX_rq1_pulse::LUX_rq1_pulse(LUX_rq1_pulse *lrq1p){\n";
	for (Int_t i=0;i<pulse_var_ct;i++) rq1_pulse_cxx << "  " << pulse_var_names[i] << " = lrq1p->Get_" << pulse_var_names[i] << "();\n";
	rq1_pulse_cxx << "}\n\n";
	rq1_pulse_cxx << "LUX_rq1_pulse::~LUX_rq1_pulse(){}\n\n";
	rq1_pulse_cxx << "void LUX_rq1_pulse::Fill(";
	for (Int_t i=0;i<pulse_var_ct-1;i++){
	  rq1_pulse_cxx << stupider_method(pulse_var_types[i]) << " " << name << i <<", ";
	}
	rq1_pulse_cxx << stupider_method(pulse_var_types[pulse_var_ct-1]) << " " << name << (pulse_var_ct-1) << "){\n";
	for (Int_t i=0;i<pulse_var_ct;i++){
	  rq1_pulse_cxx << "  " << pulse_var_names[i] << " = " << name << i << ";\n";
	}
	//for (Int_t i=0;i<pulse_var_ct;i++) rq1_pulse_cxx << "  " << pulse_var_names[i] << " = " << Char_t(name+i) << ";\n";
	rq1_pulse_cxx << "}\n\n";
	rq1_pulse_cxx << "void LUX_rq1_pulse::Fill(Char_t**a, UInt_t* b, Double_t* c){\n";
	for (Int_t i=0;i<pulse_var_ct;i++){
	  if (!strcmp(pulse_var_types[i],"char")) rq1_pulse_cxx << "  " << pulse_var_names[i] << " = a[" << i << "];\n";
	  else if (!strcmp(pulse_var_types[i],"uint32")) rq1_pulse_cxx << "  " << pulse_var_names[i] << " = b[" << i << "];\n";
	  else if (!strcmp(pulse_var_types[i],"double")) rq1_pulse_cxx << "  " << pulse_var_names[i] << " = c[" << i << "];\n";
	}
	rq1_pulse_cxx << "}";
	rq1_pulse_cxx.close();
      }
      delete[] header;
    }
  }
  else cout << "File could not be opened" << endl;
}

int main(int argc, char **argv)
{
  Int_t result = 0;
  if (argc > 0)
    result = set_up_header_files(argv[1]);
  else{
    printf("Usage: crazy_class_builder <binary file>");
    printf(" where <binary_file> is the name of the merged file");
    printf(" to be read\n");
  }
}
