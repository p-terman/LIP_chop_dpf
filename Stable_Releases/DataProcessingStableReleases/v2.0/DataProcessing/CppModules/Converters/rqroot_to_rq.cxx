/*
This program takes the .rq1.root file produce by the LUX_First_Pass and converts
the files into the .rq1 format. This program is a stand alone and not meant to 
run interactively in the root environment. 

LUX Collaboration
Group: ROOT Analysis Chain
Author: Sergey Uvarov
Version: 0.1
*/
#include <string>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <iostream>
#include "TObject.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TLeaf.h"
#include "TBranch.h"
#include "TTree.h"
#include "TFile.h"

using namespace std;

typedef struct {
  string name;
  string type;
  vector<string> strdims;
  vector<int> dims;
}Leaf;
//______________________________________________________________________________
string GetType(string objtype, bool Unsigned) {  
  // Returns c++ types from TLeaf types.
  string type = "unknown";
  if      (!objtype.compare("TLeafC"))  type = "string";
  else if (!objtype.compare("TLeafB"))  type = "char";
  else if (!objtype.compare("TLeafS"))  type = "short";
  else if (!objtype.compare("TLeafI"))  type = "int";
  else if (!objtype.compare("TLeafF"))  type = "float";
  else if (!objtype.compare("TLeafD"))  type = "double";
  else if (!objtype.compare("TLeafL"))  type = "long long";
  else if (!objtype.compare("TLeafO"))  type = "bool";
  if (Unsigned) {
    type.insert(0, "unsigned ");
  }
  return type; 
}
//______________________________________________________________________________
void GetDimensions(string name, vector<string> &dimensions) {
  // finds the dimensions written in string form in the name. 
  dimensions.clear();
  size_t b_pos1 = 0, b_pos2 = 0;
  do {
    b_pos1 = name.find('[', b_pos2);
    b_pos2 = name.find(']', b_pos1);
    if (b_pos1 == string::npos or b_pos2 == string::npos) break;
    b_pos1++;
    string dim = name.substr(b_pos1, b_pos2-b_pos1);
    if (dim.size() > 0) dimensions.push_back(dim);
  }while(b_pos1 != string::npos and b_pos2 != string::npos);
}
//______________________________________________________________________________
void BuildLeaves(TTree *tree, vector<Leaf> &leaves) {
  // Transfers/interprits data from the TLeafs in the TTree. It gets the names,
  // data types, and dimensions. 
  leaves.clear();  
  TObjArray *treeleafs = tree->GetListOfLeaves();
  for (int i=0; i<treeleafs->GetEntries(); i++) {
    TLeaf *tleaf = (TLeaf*)(*treeleafs)[i];
    Leaf leaf;
    leaf.name = tleaf->GetName();
    leaf.type = GetType(tleaf->ClassName(), tleaf->IsUnsigned());
    GetDimensions(tleaf->GetTitle(), leaf.strdims);
    if (!leaf.strdims.size()) leaf.strdims.push_back("1");
    leaf.dims.resize(leaf.strdims.size());
    leaves.push_back(leaf);
  }  
  for (size_t j=0; j<leaves.size(); j++) {
    for (size_t k=0; k<leaves[j].strdims.size(); k++) {
      int temp = atoi( leaves[j].strdims[k].c_str() );      
      if (temp == 0 and leaves[j].strdims[k].compare("0"))
        for (int i=0; i<treeleafs->GetEntries(); i++) {
          TLeaf *tleaf = (TLeaf*)(*treeleafs)[i];
          if (!leaves[j].strdims[k].compare(tleaf->GetName()))
            temp = tleaf->GetMaximum();
        }
      leaves[j].dims[k] = temp;
    }
  }
}
//______________________________________________________________________________
string CreateHeader(vector<Leaf> leaves) {
  // Returns a header block string in the .rq1 files 
  string header = "";
  for (size_t i=0; i<leaves.size(); i++) {
    string dims = "";
    for (size_t j=0; j<leaves[i].dims.size(); j++) {
      char dim[50];
      if (j+1 == leaves[i].dims.size())
        sprintf(dim, "%i", leaves[i].dims[j]);
      else
        sprintf(dim, "%i,", leaves[i].dims[j]);
      dims.append(dim);
    }    
    char variable[200];
    sprintf(variable, "%s;%s;%s;", leaves[i].name.c_str(), 
                                   leaves[i].type.c_str(), 
                                   dims.c_str());
    header.append(variable);
  }
  return header;
}
//______________________________________________________________________________
int ColumnToRow(int linear_col, vector<int> dims) {
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
void WriteEntry(TTree *tree, vector<Leaf> leaves, ostream &os) {
  // writes an entry from the tree into the ostream with column-major order. 
  TObjArray *treeleafs = tree->GetListOfLeaves();
  for (int l=0; l<treeleafs->GetEntries(); l++) {
    TLeaf *tleaf = (TLeaf*)(*treeleafs)[l];
    //os << leaves[l].name << " | "
    //   << leaves[l].type << " | "
    //   << tleaf->GetLenType() << " | ";    
    
    int length = tleaf->GetNdata();
    int maxlength = 1;
    for (int i=0; i<(int)leaves[l].dims.size(); i++){
      //os << leaves[l].dims[i] << " ";
      maxlength *= leaves[l].dims[i];
    } 
    //os << "\n";
    string type = tleaf->ClassName();
    if (!type.compare("TLeafC")) {
      const char *value = (const char*) tleaf->GetValuePointer();
      size_t len = strlen(value);
      //os << value << " ";
      os.write(value, len);
    }
    else if (!type.compare("TLeafB")) {
      char* pointer = (char*) tleaf->GetValuePointer();
      for (int i=0; i<maxlength; i++) {
        int pos = ColumnToRow(i, leaves[l].dims);
        char value = 0;
        if (pos < length)
          value = pointer[pos];
        //os << value << " ";
        os.write((char*)&value, tleaf->GetLenType());
      }
    }
    else if (!type.compare("TLeafS")) {
      short* pointer = (short*) tleaf->GetValuePointer();
      for (int i=0; i<maxlength; i++) {
        int pos = ColumnToRow(i, leaves[l].dims);
        short value = 0;
        if (pos < length)
          value = pointer[pos];
        //os << value << " ";
        os.write((char*)&value, tleaf->GetLenType());
      }
    }
    else if (!type.compare("TLeafI")) {
      int* pointer = (int*) tleaf->GetValuePointer();
      for (int i=0; i<maxlength; i++) {
        int pos = ColumnToRow(i, leaves[l].dims);
        int value = 0;
        if (pos < length)
          value = pointer[pos];
        //os << value << " ";
        os.write((char*)&value, tleaf->GetLenType());
      }
    }
    else if (!type.compare("TLeafF")) {
      float* pointer = (float*) tleaf->GetValuePointer();
      for (int i=0; i<maxlength; i++) {
        int pos = ColumnToRow(i, leaves[l].dims);
        float value = 0;
        if (pos < length)
          value = pointer[pos];
        //os << value << " ";
        os.write((char*)&value, tleaf->GetLenType());
      }
    }
    else if (!type.compare("TLeafD")) {
      double* pointer = (double*) tleaf->GetValuePointer();
      for (int i=0; i<maxlength; i++) {
        int pos = ColumnToRow(i, leaves[l].dims);
        double value = 0;
        if (pos < length)
          value = pointer[pos];
        //os << value << " ";
        os.write((char*)&value, tleaf->GetLenType());
      }
    }
    else if (!type.compare("TLeafL")) {
      long long* pointer = (long long*) tleaf->GetValuePointer();
      for (int i=0; i<maxlength; i++) {
        int pos = ColumnToRow(i, leaves[l].dims);
        long long value = 0;
        if (pos < length)
          value = pointer[pos];
        //os << value << " ";
        os.write((char*)&value, tleaf->GetLenType());
      }
    }
    else if (!type.compare("TLeafO")) {
      bool* pointer = (bool*) tleaf->GetValuePointer();
      for (int i=0; i<maxlength; i++) {
        int pos = ColumnToRow(i, leaves[l].dims);
        bool value = 0;
        if (pos < length)
          value = pointer[pos];
        //os << value << " "; 
        os.write((char*)&value, tleaf->GetLenType());
      }
    }
    //os << "\n"; 
  }
}
//______________________________________________________________________________
int main(int argc, char *argv[]) {
  
  for (int i=1; i<argc; i++) {    
    // get filename and check if extension is correct
    string filename = argv[i];
    string extension = filename;
    extension.erase(0, extension.rfind(".rq"));
    if (extension.compare(".rq1.root") and extension.compare(".rq2.root")
      and extension.compare(".rq.root")){
      cerr << "Error: "<< filename  << " is not a .rq<>.root file\n";  
      continue;
    }
    
    // open root file and check if zombie
    TFile *file = NULL;
    file = new TFile(filename.c_str());
    if (!file or file->IsZombie()) {
      delete file;
      cerr << "Error: " << filename << " is empty or Zombie\n";
      continue;
    }
    
    // create necessary variables
    vector<Leaf> leaves;
    string headerstring;
    int xmlstringlength;
    unsigned short stringlength;
    int nlines;
    
    // create output file 
    string outfilename = filename;
    outfilename.erase(outfilename.find(".root"));
    ofstream outfile(outfilename.c_str(), ios_base::binary);
    
    // write endianness 
    unsigned int endian = 0x01020304;
    outfile.write((char*)&endian, sizeof(endian));
    
    // get and write xml 
    TObjString *xml = (TObjString*) file->Get("xml");
    string xmlstring = xml->GetString().Data();
    xmlstringlength = xmlstring.size();
    outfile.write((char*)&xmlstringlength, sizeof(xmlstringlength));
    outfile.write(xmlstring.c_str(), xmlstringlength);
    delete xml;
    xml = NULL;
    
    // get and write header string 
    TTree *header = (TTree*) file->Get("header");
    header->GetEntry(0);
    BuildLeaves(header, leaves);
    headerstring = CreateHeader(leaves);
    stringlength = headerstring.size();
    nlines = header->GetEntries();
    outfile.write((char*)&stringlength, sizeof(stringlength));
    outfile.write(headerstring.c_str(), stringlength);
    outfile.write((char*)&nlines, sizeof(nlines));
    
    // get and write header information
    for (int e=0; e<header->GetEntries(); e++) {
      header->GetEntry(e);
      WriteEntry(header, leaves, outfile);
    }
    delete header;
    header = NULL;
    
    // get and write events string 
    TTree *events = (TTree*) file->Get("events");
    events->GetEntry(0);
    BuildLeaves(events, leaves);
    headerstring = CreateHeader(leaves);
    stringlength = headerstring.size();
    nlines = events->GetEntries();
    outfile.write((char*)&stringlength, sizeof(stringlength));
    outfile.write(headerstring.c_str(), stringlength);
    outfile.write((char*)&nlines, sizeof(nlines));
    
    // get and write events information
    for (int e=0; e<events->GetEntries(); e++) {
      events->GetEntry(e);
      WriteEntry(events, leaves, outfile);
    }
    delete events;
    events = NULL;
    
    // get and write events string 
    TTree *livetime = (TTree*) file->Get("livetime");
    livetime->GetEntry(0);
    BuildLeaves(livetime, leaves);
    headerstring = CreateHeader(leaves);
    stringlength = headerstring.size();
    nlines = livetime->GetEntries();
    outfile.write((char*)&stringlength, sizeof(stringlength));
    outfile.write(headerstring.c_str(), stringlength);
    outfile.write((char*)&nlines, sizeof(nlines));
    
    // get and write events information
    for (int e=0; e<livetime->GetEntries(); e++) {
      livetime->GetEntry(e);
      WriteEntry(livetime, leaves, outfile);
    }
    delete livetime;
    livetime = NULL;
    
    // close the files 
    outfile.close();
    file->Close();      
  }  
  return 0;
}
