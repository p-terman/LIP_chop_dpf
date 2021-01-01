#ifndef __RQ_File_IO_HH__ 
#define __RQ_File_IO_HH__ 1
/*
Author: Sergey Uvarov
Date: 12/21/12
Version: 0.9

This RQ_File_IO class does 3 things: reads a RQ file, writes a RQ file, and maps
the RQ_File_IO framework to the LUX framework. The LUX framework is an 
intermediate framework where both an evt file and rq file. The LUX framework is 
where all the processing of the data will take place. 

General RQ File structure
  The rq file structure is broken down into 4 sections
    1st - xml string which contains all the input values to generate the set
    2nd - rq header block which contains very basic file information
    3rd - rq events block which contains all the events and pulses rqs. 
    4th - rq livetime block which contains the livetime of the file. 

What is a RQ Block 
  An rq block is written in binary format. A "header_string" describes
  the name, type, and dimension(s) of all the rqs written in the rq block. 
  Each "header_string" is a collection of "variable_string"s where the 
  "variable_string" format is as follows
          <name>;<type>;<dim1>,...,<dimN>;
  Parsing the "header_string" into "variable_string"s is needed before any 
  "data" can be read in. A mapping between the "type" to the corresponing 
  casting type and to know the number of "bytes". The dimension(s) subsection 
  of the a "variable_string" also needs to be parsed generically.    


*/
#include <fstream> 
#include <string>
#include <vector>
#include <iostream>
#include <cstring>

using std::ifstream;
using std::ofstream;
using std::filebuf;
using std::vector;
using std::string;
using std::memcpy;
using std::cout;
using std::endl;
using std::cerr;
using std::ostream;

//______________________________________________________________________________
class RQFileIO;
class RQBlock;
class Rq;
class RqC;
class RqO;
class RqS;
class RqUS;
class RqI;
class RqUI;
class RqLL;
class RqUL;
class RqF;
class RqD;
Rq* TypeMap(string rqtype);
vector<size_t> Dimensions(string rqdim);
//______________________________________________________________________________
class Rq {
  protected:
    friend class RQBlock;
    int RowToColumn(size_t i, size_t j);
    size_t bytes, total_bytes, linear_dim, type;
    vector<size_t> dimensions;
    double initial_value;
  public:
    Rq() { initial_value=0;}
    virtual ~Rq() {;}
    
    void Build(string n, string t, string d, double value=0);
    virtual void Clear() { return; }
    size_t GetBytes() const { return bytes; }
    size_t GetTotalBytes() const { return total_bytes; }
    size_t GetLinearDim() const { return linear_dim; }
    size_t GetDim1() const { 
      if (dimensions.size()>0) return dimensions[0];
      else return 0;
    }
    size_t GetDim2() const { 
      if (dimensions.size()>1) return dimensions[1];
      else return 0;
    }
    virtual string GetString() { return 0; }
    virtual long long GetInt(size_t i=0, size_t j=0) { return 0; }
    virtual double GetDouble(size_t i=0, size_t j=0) { return 0; }
    virtual void GetRaw(char* raw, size_t i) { return; }
    void Print(ostream &ss=cout);
    virtual void SetRaw(char* raw, size_t i) { return; }
    virtual void SetString(string s) { return; }
    virtual void SetInt(long long x, size_t i=0, size_t j=0) { return; }
    virtual void SetDouble(double x, size_t i=0, size_t j=0) { return; }
    
    string rqname, rqtype, rqdim;
};
//______________________________________________________________________________
class RQBlock {
  protected:
    friend class RQFileIO;   
    unsigned short header_length;
    unsigned int nsequences;
    vector<Rq*> rqs;
    size_t BuildHeader();
    size_t GetHeaderIndex(size_t bytes);
  public:
    RQBlock();
    ~RQBlock();
    
    size_t AddRQ(string rqname, string rqtype, string rqdim, double initialize=0);
    void Clear();
    Rq* Get(string rqname);
    unsigned int GetNSequences() const { return nsequences; }
    size_t NumberRQs() const { return rqs.size(); }
    void SetNSequences(unsigned int n) { nsequences=n; }
    void PrintRQs(ostream &ss=cout) {
      for (size_t i=0; i<rqs.size(); i++) 
        rqs[i]->Print(ss);
    }
    
    string header_string;
};
//______________________________________________________________________________
class RQFileIO {
  protected:
    vector< vector<char> > raw_header;
    vector< vector<char> > raw_events;
    vector< vector<char> > raw_livetime;
    size_t header_size, event_size, livetime_size;
    size_t current_header, current_event, current_livetime;
    bool first_header, first_event, first_livetime;
    vector<string> ParseString(string header);
  public:
    RQFileIO();
    ~RQFileIO();
    
    string xml;
    RQBlock header;
    RQBlock events;
    RQBlock livetime;
    
    void AddHeader(RQFileIO *from, size_t sequence=0);
    void AddEvents(RQFileIO *from, size_t sequence=0);
    void AddLivetime(RQFileIO *from, size_t sequence=0);
    void Copy(RQFileIO *from);
    void CopyStructure(RQFileIO *from);
    void CopyXML(RQFileIO* from);
    void Clear();
    bool GetHeader(size_t sequence=0);
    bool GetEvent(size_t sequence=0);
    bool GetLivetime(size_t sequence=0);
    bool ReadFile(string filename);
    bool WriteFile(string filename);    
};
//______________________________________________________________________________
class RqC : public Rq {
  protected:
    vector<char> data;
  public:
    RqC() { bytes=1; type=1; }
    ~RqC() { data.clear(); }
    void Clear() { data.clear(); data.resize(linear_dim, initial_value);  }  
    string GetString();
    long long GetInt(size_t i=0, size_t j=0);
    double GetDouble(size_t i=0, size_t j=0);
    void GetRaw(char* raw, size_t i) { raw[0] = data[i]; }
    void SetRaw(char* raw, size_t i) { data[i] = raw[0]; }
    void SetString(string s);
    void SetInt(long long x, size_t i=0, size_t j=0);
    void SetDouble(double x, size_t i=0, size_t j=0);
};
//______________________________________________________________________________
class RqO : public Rq {
  protected:
    vector<bool> data;
  public:
    RqO() { bytes=1; type=2; }
    ~RqO() { data.clear(); }
    void Clear() { data.clear(); data.resize(linear_dim, initial_value); }
    string GetString();
    long long GetInt(size_t i=0, size_t j=0);
    double GetDouble(size_t i=0, size_t j=0);
    void GetRaw(char* raw, size_t i) { raw[0] = (char)data[i]; }
    void SetRaw(char* raw, size_t i) { data[i] = (bool)raw[0]; }
    void SetString(string s);
    void SetInt(long long x, size_t i=0, size_t j=0);
    void SetDouble(double x, size_t i=0, size_t j=0);
};
//______________________________________________________________________________
class RqS : public Rq {
  protected:
    vector<short> data;
  public:
    RqS() { bytes=2; type=2; }
    ~RqS() { data.clear(); }
    void Clear() { data.clear(); data.resize(linear_dim, initial_value); }
    string GetString();
    long long GetInt(size_t i=0, size_t j=0);
    double GetDouble(size_t i=0, size_t j=0);
    void GetRaw(char* raw, size_t i) { memcpy(raw, &data[i], bytes); }
    void SetRaw(char* raw, size_t i) { memcpy(&data[i], raw, bytes); }
    void SetString(string s);
    void SetInt(long long x, size_t i=0, size_t j=0);
    void SetDouble(double x, size_t i=0, size_t j=0);
};
//______________________________________________________________________________
class RqUS : public Rq {
  protected:
    vector<unsigned short> data;
  public:
    RqUS() { bytes=2; type=2; }
    ~RqUS() { data.clear(); }
    void Clear() { data.clear(); data.resize(linear_dim, initial_value); }
    string GetString();
    long long GetInt(size_t i=0, size_t j=0);
    double GetDouble(size_t i=0, size_t j=0);
    void GetRaw(char* raw, size_t i) { memcpy(raw, &data[i], bytes); }
    void SetRaw(char* raw, size_t i) { memcpy(&data[i], raw, bytes); }
    void SetString(string s);
    void SetInt(long long x, size_t i=0, size_t j=0);
    void SetDouble(double x, size_t i=0, size_t j=0);
};
//______________________________________________________________________________
class RqI : public Rq {
  protected:
    vector<int> data;
  public:
    RqI() { bytes=4; type=2; }
    ~RqI() { data.clear(); }
    void Clear() { data.clear(); data.resize(linear_dim, initial_value); }
    string GetString();
    long long GetInt(size_t i=0, size_t j=0);
    double GetDouble(size_t i=0, size_t j=0);
    void GetRaw(char* raw, size_t i) { memcpy(raw, &data[i], bytes); }
    void SetRaw(char* raw, size_t i) { memcpy(&data[i], raw, bytes); }
    void SetString(string s);
    void SetInt(long long x, size_t i=0, size_t j=0);
    void SetDouble(double x, size_t i=0, size_t j=0);
};
//______________________________________________________________________________
class RqUI : public Rq {
  protected:
    vector<unsigned int> data;
  public:
    RqUI() { bytes=4; type=2; }
    ~RqUI() { data.clear(); }
    void Clear() { data.clear(); data.resize(linear_dim, initial_value); }
    string GetString();
    long long GetInt(size_t i=0, size_t j=0);
    double GetDouble(size_t i=0, size_t j=0);
    void GetRaw(char* raw, size_t i) { memcpy(raw, &data[i], bytes); }
    void SetRaw(char* raw, size_t i) { memcpy(&data[i], raw, bytes); }
    void SetString(string s);
    void SetInt(long long x, size_t i=0, size_t j=0);
    void SetDouble(double x, size_t i=0, size_t j=0);
};
//______________________________________________________________________________
class RqLL : public Rq {
  protected:
    vector<long long> data;
  public:
    RqLL() { bytes=8; type=2; }
    ~RqLL() { data.clear(); }
    void Clear() { data.clear(); data.resize(linear_dim, initial_value); }
    string GetString();
    long long GetInt(size_t i=0, size_t j=0);
    double GetDouble(size_t i=0, size_t j=0);
    void GetRaw(char* raw, size_t i) { memcpy(raw, &data[i], bytes); }
    void SetRaw(char* raw, size_t i) { memcpy(&data[i], raw, bytes); }
    void SetString(string s);
    void SetInt(long long x, size_t i=0, size_t j=0);
    void SetDouble(double x, size_t i=0, size_t j=0);
};
//______________________________________________________________________________
class RqUL : public Rq {
  protected:
    vector<unsigned long long> data;
  public:
    RqUL() { bytes=8; type=2; }
    ~RqUL() { data.clear(); }
    void Clear() { data.clear(); data.resize(linear_dim, initial_value); }
    string GetString();
    long long GetInt(size_t i=0, size_t j=0);
    double GetDouble(size_t i=0, size_t j=0);
    void GetRaw(char* raw, size_t i) { memcpy(raw, &data[i], bytes); }
    void SetRaw(char* raw, size_t i) { memcpy(&data[i], raw, bytes); }
    void SetString(string s);
    void SetInt(long long x, size_t i=0, size_t j=0);
    void SetDouble(double x, size_t i=0, size_t j=0);
};
//______________________________________________________________________________
class RqF : public Rq {
  protected:
    vector<float> data;
  public:
    RqF() { bytes=4; type=3; }
    ~RqF() { data.clear(); }
    void Clear() { data.clear(); data.resize(linear_dim, initial_value); }
    string GetString();
    long long GetInt(size_t i=0, size_t j=0);
    double GetDouble(size_t i=0, size_t j=0);
    void GetRaw(char* raw, size_t i) { memcpy(raw, &data[i], bytes); }
    void SetRaw(char* raw, size_t i) { memcpy(&data[i], raw, bytes); }
    void SetString(string s);
    void SetInt(long long x, size_t i=0, size_t j=0);
    void SetDouble(double x, size_t i=0, size_t j=0);
};
//______________________________________________________________________________
class RqD : public Rq {
  protected:
    vector<double> data;
  public:
    RqD() { bytes=8; type=3; }
    ~RqD() { data.clear(); }
    void Clear() { data.clear(); data.resize(linear_dim, initial_value); }
    string GetString();
    long long GetInt(size_t i=0, size_t j=0);
    double GetDouble(size_t i=0, size_t j=0);
    const double  operator() (size_t i=0, size_t j=0) {return  GetDouble(i,j);}
    void GetRaw(char* raw, size_t i) { memcpy(raw, &data[i], bytes); }
    void SetRaw(char* raw, size_t i) { memcpy(&data[i], raw, bytes); }
    void SetString(string s);
    void SetInt(long long x, size_t i=0, size_t j=0);
    void SetDouble(double x, size_t i=0, size_t j=0);
};
#endif
