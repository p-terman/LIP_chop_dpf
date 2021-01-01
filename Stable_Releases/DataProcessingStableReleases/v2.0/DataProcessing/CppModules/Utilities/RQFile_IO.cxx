#include "RQFile_IO.hh"
#include <cstdlib>

using std::atoi;

#define EndianByteSwap(x) ByteSwap((char*)&x, sizeof(x))
#define readbinary(x) read((char*)&x, sizeof(x))
#define readstring(x,y) read((char*)x, y)
#define writebinary(x) write((char*)&x, sizeof(x))
#define writestring(x,y) write((char*)x, y)

//////////////   Generic Byte Swapping Algorithm ///////////////////////////////
void ByteSwap(char * b, int n){ 
 // This function swithes between big endianness to little endianness and vise 
 // versa taken from http://www.codeproject.com/KB/cpp/endianness.aspx  
  register int i = 0;
  register int j = n-1;
  while (i<j){
    std::swap(b[i], b[j]);
    i++, j--;
  }
}
//////////////   Mapping RQ Type name to Rq type /////////////////////////
Rq* TypeMap(string rqtype) {
  if      (!rqtype.compare("string") or !rqtype.compare("char")
           or !rqtype.compare("int8") or !rqtype.compare("uint8"))
    return (Rq*)(new RqC()); 
  else if (!rqtype.compare("bool") or !rqtype.compare("logical"))
    return (Rq*)(new RqO()); 
  else if (!rqtype.compare("short") or !rqtype.compare("int16"))
    return (Rq*)(new RqS()); 
  else if (!rqtype.compare("unsigned short") or !rqtype.compare("uint16"))
    return (Rq*)(new RqUS()); 
  else if (!rqtype.compare("int") or !rqtype.compare("int32"))
    return (Rq*)(new RqI()); 
  else if (!rqtype.compare("unsigned int") or !rqtype.compare("uint32"))
    return (Rq*)(new RqUI());  
  else if (!rqtype.compare("long long") or !rqtype.compare("int64"))
    return (Rq*)(new RqLL());  
  else if (!rqtype.compare("unsigned long long") or !rqtype.compare("uint64"))
    return (Rq*)(new RqUL());  
  else if (!rqtype.compare("float") or !rqtype.compare("single"))
    return (Rq*)(new RqF()); 
  else if (!rqtype.compare("double"))
    return (Rq*)(new RqD());     
  return NULL;
}
//////////////   Convert RQ Dimension name to real dimensions //////////////////
vector<size_t> Dimensions(string rqdim) {
  vector<size_t> dim;
  size_t start=0, end=0;
  while (start < rqdim.size() and end < rqdim.size()) {
    end = rqdim.find(',', start);
    dim.push_back( atoi(rqdim.substr(start, end-start).c_str()) );
    start = end+1;
  }
  return dim;
}
//////////////   Reading a RQ File class   /////////////////////////////////////
RQFileIO::RQFileIO() :
  header_size(0),
  event_size(0),
  livetime_size(0),
  current_header(0),
  current_event(0),
  current_livetime(0),
  first_header(1),
  first_event(1),
  first_livetime(1)
{ }
//______________________________________________________________________________
RQFileIO::~RQFileIO() {
  raw_header.clear();
  raw_events.clear();
  raw_livetime.clear();
  xml.clear();  
}
//______________________________________________________________________________
void RQFileIO::AddHeader(RQFileIO *from, size_t sequence) {
	raw_header.push_back(from->raw_header[sequence]);
	header.nsequences++;
	if (first_header) {
		GetHeader(header.nsequences-1);
		first_header = false;
	}
}
//______________________________________________________________________________
void RQFileIO::AddEvents(RQFileIO *from, size_t sequence) {
	raw_events.push_back(from->raw_events[sequence]);
	events.nsequences++;
	if (first_event) {
		GetEvent(events.nsequences-1);
		first_event = false;
	}
}
//______________________________________________________________________________
void RQFileIO::AddLivetime(RQFileIO *from, size_t sequence) {
	raw_livetime.push_back(from->raw_livetime[sequence]);
	livetime.nsequences++;
	if (first_livetime) {
		GetLivetime(livetime.nsequences-1);
		first_livetime = false;
	}
}
//______________________________________________________________________________
void RQFileIO::Copy(RQFileIO *from) {
	CopyStructure(from);
	raw_header = from->raw_header;
  raw_events = from->raw_events;
  raw_livetime = from->raw_livetime;
	CopyStructure(from);
	header.nsequences = from->header.nsequences;
	events.nsequences = from->events.nsequences;
	livetime.nsequences = from->livetime.nsequences;
	first_header = first_event = first_livetime = true;
	if (header.nsequences) {
		GetHeader(0);
		first_header = false;
	}
	if (events.nsequences) {
		GetEvent(0);
		first_event = false;
	}
	if (livetime.nsequences) {
		GetLivetime(0);
		first_livetime = false;
	}
}
//______________________________________________________________________________
void RQFileIO::CopyStructure(RQFileIO *from) {
	// Create all the rqs for the header string 
  vector<string> str = ParseString(from->header.header_string);
  header_size = 0;
  for (size_t i=0; i<str.size(); i+=3) {
    header_size += header.AddRQ(str[i+0], str[i+1], str[i+2]);
  }
	// Create all the rqs for the events
  str = ParseString(from->events.header_string);
  event_size = 0;
  for (size_t i=0; i<str.size(); i+=3) {
    event_size += events.AddRQ(str[i+0], str[i+1], str[i+2]);
  }
	// Create all the rqs for the livetime
  str = ParseString(from->livetime.header_string);
  livetime_size = 0;
  for (size_t i=0; i<str.size(); i+=3) {
    livetime_size += livetime.AddRQ(str[i+0], str[i+1], str[i+2]);
  }
}
//______________________________________________________________________________
void RQFileIO::CopyXML(RQFileIO* from) {
	xml = from->xml;
}
//______________________________________________________________________________
void RQFileIO::Clear() {
  raw_header.clear();
  raw_events.clear();
  raw_livetime.clear();
  xml.clear();
  header.Clear();
  events.Clear();
  livetime.Clear();
  header_size = event_size = livetime_size = 0;
  current_header = current_event = current_livetime = 0;
  first_header = first_event = first_livetime = true;
}
//______________________________________________________________________________
bool RQFileIO::GetHeader(size_t sequence) {
  // check if sequence exist
  if (sequence>=(size_t)header.nsequences) return false;
  
  // check if the nsequences changed
  size_t size = header.BuildHeader();
  size_t start, bytes, dim, offset = 0;
  Rq* tmp;
  if ((size_t)header.nsequences != raw_header.size()) {
  	size_t previous = raw_header.size();
  	raw_header.resize(header.nsequences, vector<char>(size,0));
  	
  	for (size_t i=previous; i<raw_header.size(); i++) {
      raw_header[i].resize(size, 0);
      
      // from the previous rqs to last rqs 
      offset = 0;
      for (size_t k=0; k<header.rqs.size(); k++) {
        tmp = header.rqs[k];
        bytes = tmp->GetBytes();
        dim = tmp->GetLinearDim();
        
        for (size_t j=0; j<dim; j++) {   
          tmp->GetRaw(&raw_header[i][offset + j*bytes], j);
        }
        offset += bytes*dim;
      }
    }
  }
  
  // check if any new rqs were added 
  for (size_t i=0; i<raw_header.size() and header_size<size; i++) {
    raw_header[i].resize(size, 0);
    
    // from the previous rqs to last rqs 
    offset = header_size;
    start = header.GetHeaderIndex(header_size);
    for (size_t k=start; k<header.rqs.size(); k++) {
      tmp = header.rqs[k];
      bytes = tmp->GetBytes();
      dim = tmp->GetLinearDim();
      
      for (size_t j=0; j<dim; j++) {   
        tmp->GetRaw(&raw_header[i][offset + j*bytes], j);
      }
      offset += bytes*dim;
    }
  }
  header_size = size;
  
  // store old sequence back into the raw data. 
  offset = 0;
  for (size_t i=0; i<header.rqs.size() and !first_header; i++) {
    tmp = header.rqs[i];
    bytes = tmp->GetBytes();
    dim = tmp->GetLinearDim();
    
    for (size_t j=0; j<dim; j++) {   
      tmp->GetRaw(&raw_header[current_header][offset + j*bytes], j);
    }
    offset += bytes*dim;
  }  
  
  // start reading in the file 1 rq variable at a time 
  offset = 0;
  for (size_t i=0; i<header.rqs.size(); i++) {
    tmp = header.rqs[i];
    bytes = tmp->GetBytes();
    dim = tmp->GetLinearDim();
    
    for (size_t j=0; j<dim; j++) {   
      tmp->SetRaw(&raw_header[sequence][offset + j*bytes], j);
    }
    offset += bytes*dim;   
  }
  
  // update the current_header
  current_header = sequence;
  first_header = false;
  return true;
}
//______________________________________________________________________________
bool RQFileIO::GetEvent(size_t sequence) {
  // check 
  if (sequence>=(size_t)events.nsequences) return false;
  
  // check if the nsequences changed
  size_t size = events.BuildHeader();
  size_t start, bytes, dim, offset = 0;
  Rq* tmp;
  if ((size_t)events.nsequences != raw_events.size()) {
  	size_t previous = raw_events.size();
  	raw_events.resize(events.nsequences, vector<char>(size,0));
  	
  	for (size_t i=previous; i<raw_events.size(); i++) {
      raw_events[i].resize(size, 0);
      
      // from the previous rqs to last rqs 
      offset = 0;
      for (size_t k=0; k<events.rqs.size(); k++) {
        tmp = events.rqs[k];
        bytes = tmp->GetBytes();
        dim = tmp->GetLinearDim();
        
        for (size_t j=0; j<dim; j++) {   
          tmp->GetRaw(&raw_events[i][offset + j*bytes], j);
        }
        offset += bytes*dim;
      }
    }
  }
  
  // check if any new rqs were added 
  for (size_t i=0; i<raw_events.size() and event_size<size; i++) {
    raw_events[i].resize(size, 0);
    
    // from the previous rqs to last rqs 
    offset = event_size;
    start = events.GetHeaderIndex(event_size);
    for (size_t k=start; k<events.rqs.size(); k++) {
      tmp = events.rqs[k];
      bytes = tmp->GetBytes();
      dim = tmp->GetLinearDim();
      
      for (size_t j=0; j<dim; j++) {   
        tmp->GetRaw(&raw_events[i][offset + j*bytes], j);
      }
      offset += bytes*dim;
    }
  }
  event_size = size;
  
  // store old sequence back into the raw data. 
  offset = 0;
  for (size_t i=0; i<events.rqs.size() and !first_event; i++) {
    tmp = events.rqs[i];
    bytes = tmp->GetBytes();
    dim = tmp->GetLinearDim();
    
    for (size_t j=0; j<dim; j++) {   
      tmp->GetRaw(&raw_events[current_event][offset + j*bytes], j);
    }
    offset += bytes*dim;
  }
  
  // start reading in the file 1 rq variable at a time 
  offset = 0;
  for (size_t i=0; i<events.rqs.size(); i++) {
    tmp = events.rqs[i];
    bytes = tmp->GetBytes();
    dim = tmp->GetLinearDim();
    
    for (size_t j=0; j<dim; j++) {   
      tmp->SetRaw(&raw_events[sequence][offset + j*bytes], j);
    }
    offset += bytes*dim;   
  }
  
  // update the current_event
  current_event = sequence;
  first_event = false;
  return true;
}
//______________________________________________________________________________
bool RQFileIO::GetLivetime(size_t sequence) {
  // check 
  if (sequence>=(size_t)livetime.nsequences) return false;
  
  // check if the nsequences changed
  size_t size = livetime.BuildHeader();
  size_t start, bytes, dim, offset = 0;
  Rq* tmp;
  if ((size_t)livetime.nsequences != raw_livetime.size()) {
  	size_t previous = raw_livetime.size();
  	raw_livetime.resize(livetime.nsequences, vector<char>(size,0));
  	
  	for (size_t i=previous; i<raw_livetime.size(); i++) {
      raw_livetime[i].resize(size, 0);
      
      // from the previous rqs to last rqs 
      offset = 0;
      for (size_t k=0; k<livetime.rqs.size(); k++) {
        tmp = livetime.rqs[k];
        bytes = tmp->GetBytes();
        dim = tmp->GetLinearDim();
        
        for (size_t j=0; j<dim; j++) {   
          tmp->GetRaw(&raw_livetime[i][offset + j*bytes], j);
        }
        offset += bytes*dim;
      }
    }
  }
  
  // check if any new rqs were added 
  for (size_t i=0; i<raw_livetime.size() and livetime_size<size; i++) {
    raw_livetime[i].resize(size, 0);
    
    // from the previous rqs to last rqs 
    offset = livetime_size;
    start = livetime.GetHeaderIndex(livetime_size);
    for (size_t k=start; k<livetime.rqs.size(); k++) {
      tmp = livetime.rqs[k];
      bytes = tmp->GetBytes();
      dim = tmp->GetLinearDim();
      
      for (size_t j=0; j<dim; j++) {   
        tmp->GetRaw(&raw_livetime[i][offset + j*bytes], j);
      }
      offset += bytes*dim;
    }
  }
  livetime_size = size;
  
  // store old sequence back into the raw data. 
  offset = 0;
  for (size_t i=0; i<livetime.rqs.size() and !first_livetime; i++) {
    tmp = livetime.rqs[i];
    bytes = tmp->GetBytes();
    dim = tmp->GetLinearDim();
    
    for (size_t j=0; j<dim; j++) {   
      tmp->GetRaw(&raw_livetime[current_livetime][offset + j*bytes], j);
    }
    offset += bytes*dim;
  }
  
  // start reading in the file 1 rq variable at a time 
  offset = 0;
  for (size_t i=0; i<livetime.rqs.size(); i++) {
    tmp = livetime.rqs[i];
    bytes = tmp->GetBytes();
    dim = tmp->GetLinearDim();
    
    for (size_t j=0; j<dim; j++) {   
      tmp->SetRaw(&raw_livetime[sequence][offset + j*bytes], j);
    }
    offset += bytes*dim;   
  }
  
  // update the current_livetime
  current_livetime = sequence;
  first_livetime = false;
  return true;
}
//______________________________________________________________________________
bool RQFileIO::ReadFile(string filename) {  
  // remove all previous stuff
  Clear();
  
  // open file 
  ifstream file (filename.c_str(), ifstream::binary);
    
  // determine byteswap
  bool byteswap = false;
  unsigned int endianess = 0;
  file.readbinary(endianess);
  if (endianess == 0x04030201) byteswap = true;
  
  // read in xml string 
  unsigned int xml_length = 0;
  file.readbinary(xml_length);
  if (byteswap) EndianByteSwap(xml_length);
  xml.resize(xml_length);
  file.readstring(xml.data(), xml_length);  
  
  // read the header RQBlock header 
  file.readbinary(header.header_length);
  if (byteswap) EndianByteSwap(header.header_length);
  header.header_string.resize(header.header_length);
  file.readstring(header.header_string.data(), header.header_length);
  file.readbinary(header.nsequences);
  if (byteswap) EndianByteSwap(header.nsequences);
  
  // Create all the rqs for the header string 
  vector<string> str = ParseString(header.header_string);
  if (str.size() % 3 != 0) {
    cerr << "Error: 1st RQBlock header string has the incorrect format.\n";
    return false;
  }
  header_size = 0;
  for (size_t i=0; i<str.size(); i+=3) {
    header_size += header.AddRQ(str[i+0], str[i+1], str[i+2]);
  }  
  
  // read all the header binary information 
  size_t bytes, dim, offset; // in case of byteswap is necessary
  raw_header.resize(header.nsequences, vector<char> (header_size, 0));
  for (size_t i=0; i<raw_header.size() and !byteswap; i++) // no byteswap
    file.readstring(&raw_header[i][0], header_size);
  for (size_t i=0; i<raw_header.size() and byteswap; i++) { // yes byteswap
    offset = 0;
    for (size_t j=0; j<header.rqs.size(); j++) {
      bytes = header.rqs[j]->GetBytes();
      dim = header.rqs[j]->GetLinearDim();      
      for (size_t k=0; k<dim; k++) {
        file.readstring(&raw_header[i][offset+j*bytes], bytes);
        ByteSwap(&raw_header[i][offset+j*bytes], bytes);
      }
      offset += dim*bytes;
    }
  }
  
  // load the first header block 
  if (header.nsequences != 0) GetHeader(0);
  
  // read the events RQBlock header string
  file.readbinary(events.header_length);
  if (byteswap) EndianByteSwap(events.header_length);
  events.header_string.resize(events.header_length);
  file.readstring(events.header_string.data(), events.header_length);
  file.readbinary(events.nsequences);
  if (byteswap) EndianByteSwap(events.nsequences);
  
  // Create all the rqs for the events
  str = ParseString(events.header_string);
  if (str.size() % 3 != 0) {
    cerr << "Error: 2nd RQBlock header string has the incorrect format.\n";
    return false;
  }
  event_size = 0;
  for (size_t i=0; i<str.size(); i+=3) {
    event_size += events.AddRQ(str[i+0], str[i+1], str[i+2]);
  }  
  
  // read in the whole events 
  raw_events.resize(events.nsequences, vector<char> (event_size, 0));
  for (size_t i=0; i<raw_events.size() and !byteswap; i++) 
    file.readstring(&raw_events[i][0], event_size);
  for (size_t i=0; i<raw_events.size() and byteswap; i++) { // yes byteswap
    offset = 0;
    for (size_t j=0; j<events.rqs.size(); j++) {
      bytes = events.rqs[j]->GetBytes();
      dim = events.rqs[j]->GetLinearDim();      
      for (size_t k=0; k<dim; k++) {
        file.readstring(&raw_events[i][offset+j*bytes], bytes);
        ByteSwap(&raw_events[i][offset+j*bytes], bytes);
      }
      offset += dim*bytes;
    }
  }
  
  // load the first event 
  if (events.nsequences != 0) GetEvent(0);
  
  // read the livetime RQBlock header string
  file.readbinary(livetime.header_length);
  if (byteswap) EndianByteSwap(livetime.header_length);
  livetime.header_string.resize(livetime.header_length);
  file.readstring(livetime.header_string.data(), livetime.header_length);
  file.readbinary(livetime.nsequences);
  if (byteswap) EndianByteSwap(livetime.nsequences);
  
  // Create all the rqs for the livetime
  str = ParseString(livetime.header_string);
  if (str.size() % 3 != 0) {
    cerr << "Error: 3rd RQBlock header string has the incorrect format.\n";
    return false;
  }
  livetime_size = 0;
  for (size_t i=0; i<str.size(); i+=3) {
    livetime_size += livetime.AddRQ(str[i+0], str[i+1], str[i+2]);
  }
  
  // read in the livetimes 
  raw_livetime.resize(livetime.nsequences, vector<char> (livetime_size, 0));
  for (size_t i=0; i<raw_livetime.size(); i++) 
    file.readstring(&raw_livetime[i][0], livetime_size);
  for (size_t i=0; i<raw_livetime.size() and byteswap; i++) { // yes byteswap
    offset = 0;
    for (size_t j=0; j<livetime.rqs.size(); j++) {
      bytes = livetime.rqs[j]->GetBytes();
      dim = livetime.rqs[j]->GetLinearDim();      
      for (size_t k=0; k<dim; k++) {
        file.readstring(&raw_livetime[i][offset+j*bytes], bytes);
        ByteSwap(&raw_livetime[i][offset+j*bytes], bytes);
      }
      offset += dim*bytes;
    }
  }
  
  // load the first event 
  if (livetime.nsequences != 0) GetLivetime(0);
  
  file.close();
  return true;
}
//______________________________________________________________________________
bool RQFileIO::WriteFile(string filename) {
  
  // store the previous values
  GetHeader(current_header);
  GetEvent(current_event);
  GetLivetime(current_livetime);
  
  // open file 
  ofstream file(filename.c_str(), ofstream::binary);
  
  // write the byteswap 
  unsigned int endianess = 0x01020304;
  file.writebinary(endianess);
  
  // write the xml 
  unsigned int xml_length = xml.size();
  file.writebinary(xml_length);
  file.writestring(xml.data(), xml_length);
  
  // build all the header string information for the RQBlocks
  header_size = header.BuildHeader();
  event_size = events.BuildHeader();
  livetime_size = livetime.BuildHeader();
  
  // write the header RQBlock header 
  file.writebinary(header.header_length);
  file.writestring(header.header_string.data(), header.header_length);
  file.writebinary(header.nsequences);
  
  // write the header RQBlock data 
  for (size_t i=0; i<raw_header.size(); i++) 
    file.writestring(raw_header[i].data(), header_size);
  
  // write the events RQBlock 
  file.writebinary(events.header_length);
  file.writestring(events.header_string.data(), events.header_length);
  file.writebinary(events.nsequences);
  
  // write the events RQBlock data 
  for (size_t i=0; i<raw_events.size(); i++) 
    file.writestring(raw_events[i].data(), event_size);
  
  // write the header RQBlock 
  file.writebinary(livetime.header_length);
  file.writestring(livetime.header_string.data(), livetime.header_length);
  file.writebinary(livetime.nsequences);
  
  // write the livetime RQBlock data 
  for (size_t i=0; i<raw_livetime.size(); i++) 
    file.writestring(raw_livetime[i].data(), livetime_size);
  
  file.close();
  return true;
}
//______________________________________________________________________________
vector<string> RQFileIO::ParseString(string header) {
  vector<string> str;
  size_t start=0, end=0;
  while (start < header.size() and end < header.size()) {
    end = header.find(';', start);
    str.push_back( header.substr(start, end-start).c_str() );
    start = end+1;
  }
  return str;
}
////////////////////   Manages the Access to RQ quantities /////////////////////
RQBlock::RQBlock() :
  header_length(0),
  nsequences(0)
{ }
//______________________________________________________________________________
RQBlock::~RQBlock() {
  Clear();
}
//______________________________________________________________________________
void RQBlock::Clear() {
  for (size_t i=0; i<rqs.size(); i++) {
    delete rqs[i]; 
    rqs[i] = NULL;
  }
  rqs.clear();
  header_string.clear();
  header_length = nsequences = 0;
}
//______________________________________________________________________________
size_t RQBlock::AddRQ(string rqname, string rqtype, string rqdim, double initialize) {
  // check if rqname already exist
  for (size_t i=0; i<rqs.size(); i++) {
    if (!rqname.compare(rqs[i]->rqname)) {
      cerr << "Warning: RQ name already exist\n";
      return 0;
    }    
  }  
  // get correct variable type 
  Rq* tmp = TypeMap(rqtype);
  if (!tmp) {
    cerr << "Error: Cannot add new rq. Unknown rq type found \""
         << rqtype << "\"\n";
    return 0;
  }  
  tmp->Build(rqname, rqtype, rqdim, initialize);
  rqs.push_back(tmp);  
  return tmp->GetTotalBytes();
}
//______________________________________________________________________________
size_t RQBlock::BuildHeader() {
  size_t size = 0;  
  header_string.clear();
  Rq* tmp = NULL;
  for (size_t i=0; i<rqs.size(); i++) {
    tmp = rqs[i];
    header_string += tmp->rqname;
    header_string += ';';
    header_string += tmp->rqtype;
    header_string += ';';
    header_string += tmp->rqdim;
    header_string += ';';
    size += tmp->GetTotalBytes();
  }
  header_length = header_string.size();
  return size;
}
//______________________________________________________________________________
size_t RQBlock::GetHeaderIndex(size_t bytes) {
  size_t size = 0;
  size_t index = 0;
  Rq* tmp = NULL;
  for (size_t i=0; i<rqs.size(); i++) {
    if (size >= bytes) break;
    tmp = rqs[i];
    size += tmp->GetTotalBytes();
    index++;
  }
  return index;
}
//______________________________________________________________________________
Rq* RQBlock::Get(string rqname) {
  for (size_t i=0; i<rqs.size(); i++)
    if (!rqname.compare(rqs[i]->rqname))
      return rqs[i];
  return NULL;
}
//////////////////////////   Generic RQ variable ///////////////////////////////
void Rq::Build(string n, string t, string d, double value) {
  rqname=n; rqtype=t; rqdim=d;
  dimensions = Dimensions(rqdim);
  linear_dim = 1;
  for (size_t i=0; i<dimensions.size(); i++) linear_dim *= dimensions[i];
  if (dimensions.size() == 0) linear_dim = 0;
  total_bytes = bytes*linear_dim;
  initial_value = value;
  Clear();
}
//______________________________________________________________________________
int Rq::RowToColumn(size_t i, size_t j) {
  if (dimensions.size() == 1) if (i<dimensions[0]) return i;
  if (dimensions.size() == 2) 
    if (i<dimensions[0] and j<dimensions[1]) return i+dimensions[0]*j;
  return -1;
}
//______________________________________________________________________________
void Rq::Print(ostream &ss/*=cout*/) {
  ss << rqname << ";" << rqtype << ";" << rqdim << ";\n";
  if (type == 1)
    ss << GetString() << "\n";
  else if (dimensions.size() == 2) {
    for (size_t i=0; i<dimensions[0]; i++) {
      for (size_t j=0; j<dimensions[1]; j++) {
        if (type == 2) ss << GetInt(i,j) << " ";
        else if (type == 3) ss << GetDouble(i,j) << " ";
      }
      ss << "\n";
    }
  } else if (dimensions.size() == 1) {
    for (size_t i=0; i<dimensions[0]; i++) {
      if (type == 2) ss << GetInt(i) << " ";
      else if (type == 3) ss << GetDouble(i) << " ";
    } 
    ss << "\n";
  }
}
//////////////   String (char) RQ type /////////////////////////////////////////
string RqC::GetString() {
  string tmp;
  for (size_t i=0; i<data.size(); i++) tmp += data[i];
  return tmp;
}
//______________________________________________________________________________
long long RqC::GetInt(size_t i, size_t j) {
  int index = RowToColumn(i, j);
  if (index == -1) return 0; // invalid section 
  return data[index];
}
//______________________________________________________________________________
double RqC::GetDouble(size_t i, size_t j) {
  int index = RowToColumn(i, j);
  if (index == -1) return 0; // invalid section 
  return data[index];
}
//______________________________________________________________________________
void RqC::SetString(string s) {
  Clear();
  for (size_t i=0; i<s.size() && i<linear_dim; i++) data[i] = s[i];
}
//______________________________________________________________________________
void RqC::SetInt(long long x, size_t i, size_t j) {
  int index = RowToColumn(i, j);
  if (index == -1) return; // invalid section 
  data[index] = (char)x;
}
//______________________________________________________________________________
void RqC::SetDouble(double x, size_t i, size_t j) {
  int index = RowToColumn(i, j);
  if (index == -1) return; // invalid section 
  data[index] = (char)x;
}
//////////////   Boolean (bool) RQ type ////////////////////////////////////////
string RqO::GetString() {
  string tmp;
  for (size_t i=0; i<data.size(); i++) tmp += char(data[i]);
  return tmp;
}
//______________________________________________________________________________
long long RqO::GetInt(size_t i, size_t j) {
  int index = RowToColumn(i, j);
  if (index == -1) return 0; // invalid section 
  return data[index];
}
//______________________________________________________________________________
double RqO::GetDouble(size_t i, size_t j) {
  int index = RowToColumn(i, j);
  if (index == -1) return 0; // invalid section 
  return data[index];
}
//______________________________________________________________________________
void RqO::SetString(string s) {
  Clear();
  for (size_t i=0; i<s.size() && i<linear_dim; i++) data[i] = bool(s[i]);
}
//______________________________________________________________________________
void RqO::SetInt(long long x, size_t i, size_t j) {
  int index = RowToColumn(i, j);
  if (index == -1) return; // invalid section 
  data[index] = bool(x);
}
//______________________________________________________________________________
void RqO::SetDouble(double x, size_t i, size_t j) {
  int index = RowToColumn(i, j);
  if (index == -1) return; // invalid section 
  data[index] = bool(x);
}
//////////////   Short (16bit integer) RQ type /////////////////////////////////
string RqS::GetString() {
  string tmp;
  for (size_t i=0; i<data.size(); i++) tmp += char(data[i]);
  return tmp;
}
//______________________________________________________________________________
long long RqS::GetInt(size_t i, size_t j) {
  int index = RowToColumn(i, j);
  if (index == -1) return 0; // invalid section 
  return data[index];
}
//______________________________________________________________________________
double RqS::GetDouble(size_t i, size_t j) {
  int index = RowToColumn(i, j);
  if (index == -1) return 0; // invalid section 
  return data[index];
}
//______________________________________________________________________________
void RqS::SetString(string s) {
  Clear();
  for (size_t i=0; i<s.size() && i<linear_dim; i++) data[i] = short(s[i]);
}
//______________________________________________________________________________
void RqS::SetInt(long long x, size_t i, size_t j) {
  int index = RowToColumn(i, j);
  if (index == -1) return; // invalid section 
  data[index] = short(x);
}
//______________________________________________________________________________
void RqS::SetDouble(double x, size_t i, size_t j) {
  int index = RowToColumn(i, j);
  if (index == -1) return; // invalid section 
  data[index] = short(x);
}
//////////////   Unsigned short (16bit integer) RQ type ////////////////////////
string RqUS::GetString() {
  string tmp;
  for (size_t i=0; i<data.size(); i++) tmp += char(data[i]);
  return tmp;
}
//______________________________________________________________________________
long long RqUS::GetInt(size_t i, size_t j) {
  int index = RowToColumn(i, j);
  if (index == -1) return 0; // invalid section 
  return data[index];
}
//______________________________________________________________________________
double RqUS::GetDouble(size_t i, size_t j) {
  int index = RowToColumn(i, j);
  if (index == -1) return 0; // invalid section 
  return data[index];
}
//______________________________________________________________________________
void RqUS::SetString(string s) {
  Clear();
  for (size_t i=0; i<s.size() && i<linear_dim; i++) data[i] = short(s[i]);
}
//______________________________________________________________________________
void RqUS::SetInt(long long x, size_t i, size_t j) {
  int index = RowToColumn(i, j);
  if (index == -1) return; // invalid section 
  data[index] = short(x);
}
//______________________________________________________________________________
void RqUS::SetDouble(double x, size_t i, size_t j) {
  int index = RowToColumn(i, j);
  if (index == -1) return; // invalid section 
  data[index] = short(x);
}
//////////////   Integer (32bit integer) RQ type ///////////////////////////////
string RqI::GetString() {
  string tmp;
  for (size_t i=0; i<data.size(); i++) tmp += char(data[i]);
  return tmp;
}
//______________________________________________________________________________
long long RqI::GetInt(size_t i, size_t j) {
  int index = RowToColumn(i, j);
  if (index == -1) return 0; // invalid section 
  return data[index];
}
//______________________________________________________________________________
double RqI::GetDouble(size_t i, size_t j) {
  int index = RowToColumn(i, j);
  if (index == -1) return 0; // invalid section 
  return data[index];
}
//______________________________________________________________________________
void RqI::SetString(string s) {
  Clear();
  for (size_t i=0; i<s.size() && i<linear_dim; i++) data[i] = int(s[i]);
}
//______________________________________________________________________________
void RqI::SetInt(long long x, size_t i, size_t j) {
  int index = RowToColumn(i, j);
  if (index == -1) return; // invalid section 
  data[index] = int(x);
}
//______________________________________________________________________________
void RqI::SetDouble(double x, size_t i, size_t j) {
  int index = RowToColumn(i, j);
  if (index == -1) return; // invalid section 
  data[index] = int(x);
}
//////////////   Unsigned integer (32bit integer) RQ type //////////////////////
string RqUI::GetString() {
  string tmp;
  for (size_t i=0; i<data.size(); i++) tmp += char(data[i]);
  return tmp;
}
//______________________________________________________________________________
long long RqUI::GetInt(size_t i, size_t j) {
  int index = RowToColumn(i, j);
  if (index == -1) return 0; // invalid section 
  return data[index];
}
//______________________________________________________________________________
double RqUI::GetDouble(size_t i, size_t j) {
  int index = RowToColumn(i, j);
  if (index == -1) return 0; // invalid section 
  return data[index];
}
//______________________________________________________________________________
void RqUI::SetString(string s) {
  Clear();
  for (size_t i=0; i<s.size() && i<linear_dim; i++) data[i] = int(s[i]);
}
//______________________________________________________________________________
void RqUI::SetInt(long long x, size_t i, size_t j) {
  int index = RowToColumn(i, j);
  if (index == -1) return; // invalid section 
  data[index] = int(x);
}
//______________________________________________________________________________
void RqUI::SetDouble(double x, size_t i, size_t j) {
  int index = RowToColumn(i, j);
  if (index == -1) return; // invalid section 
  data[index] = int(x);
}
//////////////   Long interger (64bit integer) RQ type /////////////////////////
string RqLL::GetString() {
  string tmp;
  for (size_t i=0; i<data.size(); i++) tmp += char(data[i]);
  return tmp;
}
//______________________________________________________________________________
long long RqLL::GetInt(size_t i, size_t j) {
  int index = RowToColumn(i, j);
  if (index == -1) return 0; // invalid section 
  return data[index];
}
//______________________________________________________________________________
double RqLL::GetDouble(size_t i, size_t j) {
  int index = RowToColumn(i, j);
  if (index == -1) return 0; // invalid section 
  return data[index];
}
//______________________________________________________________________________
void RqLL::SetString(string s) {
  Clear();
  for (size_t i=0; i<s.size() && i<linear_dim; i++) data[i] = (long long)(s[i]);
}
//______________________________________________________________________________
void RqLL::SetInt(long long x, size_t i, size_t j) {
  int index = RowToColumn(i, j);
  if (index == -1) return; // invalid section 
  data[index] = x;
}
//______________________________________________________________________________
void RqLL::SetDouble(double x, size_t i, size_t j) {
  int index = RowToColumn(i, j);
  if (index == -1) return; // invalid section 
  data[index] = x;
}
//////////////   Unsigned long integer (64bit integer) RQ type /////////////////
string RqUL::GetString() {
  string tmp;
  for (size_t i=0; i<data.size(); i++) tmp += char(data[i]);
  return tmp;
}
//______________________________________________________________________________
long long RqUL::GetInt(size_t i, size_t j) {
  int index = RowToColumn(i, j);
  if (index == -1) return 0; // invalid section 
  return data[index];
}
//______________________________________________________________________________
double RqUL::GetDouble(size_t i, size_t j) {
  int index = RowToColumn(i, j);
  if (index == -1) return 0; // invalid section 
  return data[index];
}
//______________________________________________________________________________
void RqUL::SetString(string s) {
  Clear();
  for (size_t i=0; i<s.size() && i<linear_dim; i++) data[i] = (long long)(s[i]);
}
//______________________________________________________________________________
void RqUL::SetInt(long long x, size_t i, size_t j) {
  int index = RowToColumn(i, j);
  if (index == -1) return; // invalid section 
  data[index] = x;
}
//______________________________________________________________________________
void RqUL::SetDouble(double x, size_t i, size_t j) {
  int index = RowToColumn(i, j);
  if (index == -1) return; // invalid section 
  data[index] = x;
}
//////////////   Floating (32bit) RQ type //////////////////////////////////////
string RqF::GetString() {
  string tmp;
  for (size_t i=0; i<data.size(); i++) tmp += char(data[i]);
  return tmp;
}
//______________________________________________________________________________
long long RqF::GetInt(size_t i, size_t j) {
  int index = RowToColumn(i, j);
  if (index == -1) return 0; // invalid section 
  return data[index];
}
//______________________________________________________________________________
double RqF::GetDouble(size_t i, size_t j) {
  int index = RowToColumn(i, j);
  if (index == -1) return 0; // invalid section 
  return data[index];
}
//______________________________________________________________________________
void RqF::SetString(string s) {
  Clear();
  for (size_t i=0; i<s.size() && i<linear_dim; i++) data[i] = float(s[i]);
}
//______________________________________________________________________________
void RqF::SetInt(long long x, size_t i, size_t j) {
  int index = RowToColumn(i, j);
  if (index == -1) return; // invalid section 
  data[index] = x;
}
//______________________________________________________________________________
void RqF::SetDouble(double x, size_t i, size_t j) {
  int index = RowToColumn(i, j);
  if (index == -1) return; // invalid section 
  data[index] = x;
}
//////////////   Floating (64bit) RQ type //////////////////////////////////////
string RqD::GetString() {
  string tmp;
  for (size_t i=0; i<data.size(); i++) tmp += char(data[i]);
  return tmp;
}
//______________________________________________________________________________
long long RqD::GetInt(size_t i, size_t j) {
  int index = RowToColumn(i, j);
  if (index == -1) return 0; // invalid section 
  return data[index];
}
//______________________________________________________________________________
double RqD::GetDouble(size_t i, size_t j) {
  int index = RowToColumn(i, j);
  if (index == -1) return 0; // invalid section 
  return data[index];
}
//______________________________________________________________________________
void RqD::SetString(string s) {
  Clear();
  for (size_t i=0; i<s.size() && i<linear_dim; i++) data[i] = double(s[i]);
}
//______________________________________________________________________________
void RqD::SetInt(long long x, size_t i, size_t j) {
  int index = RowToColumn(i, j);
  if (index == -1) return; // invalid section 
  data[index] = x;
}
//______________________________________________________________________________
void RqD::SetDouble(double x, size_t i, size_t j) {
  int index = RowToColumn(i, j);
  if (index == -1) return; // invalid section 
  data[index] = x;
}
