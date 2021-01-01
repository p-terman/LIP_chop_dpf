#include "LEvtUtilities.h"

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
void DelimitStringToVec(std::string str, std::vector<std::string>* blocks, 
                        std::string separator) {
  blocks->clear();
  std::string headerstring = str;

  size_t pos1 = 0, pos2 = 0;
  size_t length = headerstring.length();
  while ( pos1 < length && pos2 < headerstring.npos ){
      pos2 = headerstring.find_first_of(separator, pos1);
      if( pos2 == std::string::npos) {
        blocks->push_back(headerstring.substr(pos1, std::string::npos));
        return;
      }
      if (pos2 != pos1){
        std::string sub = headerstring.substr(pos1, pos2-pos1);
        blocks->push_back( sub );
      }
      pos1 = pos2+1;
  }

  return ;
}