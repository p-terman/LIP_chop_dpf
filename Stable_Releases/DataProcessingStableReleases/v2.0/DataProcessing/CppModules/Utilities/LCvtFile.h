#ifndef LCVTFILE_H
#define LCVTFILE_H
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include "LEvtFile.h"
#include "LCvtEvent.h"
#include "LEvtUtilities.h"


class LCvtFile {
  public:
    ////////////////////////////
    //////// Members ///////////
    ////////////////////////////
    std::ifstream input;
    
    unsigned int endianness;
    bool byteswap;
    
    unsigned int settings_string_len;
    char* settings_string;
    boost::property_tree::ptree settings_pt;

    // File Header
    unsigned int fh_endianness;
    unsigned int fh_date_time;
    unsigned int fh_location;
    unsigned int fh_number_events;
    std::vector<unsigned int> fh_gid_event_num;
    std::vector<unsigned int> fh_gid_byte_loc;

    // Livetime Header
    unsigned short lt_num_sequences;
    std::vector<unsigned long long> lh_timestamp_latch;
    std::vector<unsigned long long> lh_timestamp_end;

    // Event GID
    std::vector<LCvtEvent*> cvt_events;

    //  The file loaded into the object
    std::string input_filename;
    
    ////////////////////////////
    //////// Methods ///////////
    ////////////////////////////

    LCvtFile();
    LCvtFile(std::string path, std::string cvt_file);
    LCvtFile(LEvtFile *evt);
    ~LCvtFile();
    
    void GetSettings();
    void GetFileHeader();
    void GetLivetime();
    void GetEvents();
    void WriteCvtFile(std::string path, std::string cvt_file);
    void UpdateFileLocations();
    void Convert(std::vector<double> pmt_gains);
};

#endif
