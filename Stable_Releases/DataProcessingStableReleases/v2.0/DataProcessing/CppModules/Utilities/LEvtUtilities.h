#ifndef LEVTUTILITIES_H
#define LEVTUTILITIES_H
#include <vector>
#include <string>

// A shortcut to the endian byteswapping
#define EndianByteSwap(x) ByteSwap((unsigned char *) &x,sizeof(x))

// The following two help read in binary data without cluttering the screen.
#define readval(x) read((char*)&x,sizeof(x));if(byteswap)EndianByteSwap(x)
#define readstr(x,y) read((char*)x,y);if(byteswap)EndianByteSwap(x)
#define readvec(x) read(reinterpret_cast<char*>(&(x.front())),x.size()*sizeof(x.front()));if(byteswap)EndianByteSwap(x)
//#define readvec2d(x) read((char*)&(x.front()),x.size()*(x.front().size())*sizeof(x.front().front()));

void DelimitStringToVec(std::string str, std::vector<std::string>* blocks, std::string separator=";");

void ByteSwap(unsigned char * b, int n);

#endif