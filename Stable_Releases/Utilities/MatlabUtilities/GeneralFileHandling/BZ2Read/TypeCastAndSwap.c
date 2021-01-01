#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* This c code is intended for use as a mex file to perform typecasting and
 * byte swapping -- it may be replaced by the typecast and swapbytes functions
 * in matlab 7.1 (R14SP3) and later.
 *
 * Usage is:
 * data = TypeCastAndSwap(inputdata, format, byteswap_input, byteswap_output)
 * inputdata:  Matlab vector (numeric or char)
 * format:     desired type of output data ([u]int(8|16|32|64), double, or single
 *                                             also understands c-style type naming)
 *                     if the input format is blank or not understood, the format is not
 *                     changed.
 * byteswap_input:   boolean, true to reverse byte order based on the input element size
 * byteswap_output:  boolean, true to reverse byte order based on the output element size
 * data:       output Matlab vector (same singleton dimension as input -- if input is scalar,
 *                                   output will be a column vector)
 *
 * The function creates a matlab array of the desired type and copies the data 
 * to the new array -- if the input is the wrong length to give an integer number
 * of output elements, the data is padded with zeros and a warning is printed.
 * After filling the new array, the bytes in each element are reversed according to the
 * byteswap arguments.  If both byteswap input and output are true, the input byteswapping
 * is done first (it's left as an excercise to the reader to show that this doesn't matter).
 * If the code is instructed to byteswap single-byte data, a warning is issued.
 *
 * CED, 06/14/06
 * Release Version 1.0
 * */

void FlipBytes(void* input_buffer, int elementsize, int count) {
  char temp;
  char* buffer = NULL;
  int i;
  int j;
  int halfelementsize = (elementsize/2);
  buffer = (char*) input_buffer;
  for (i=0; i<((count-1)*elementsize+1); i += elementsize)
    for (j = 0; j < halfelementsize; j++) {
      temp = buffer[i+j];
      buffer[i+j] = buffer[i+elementsize-j-1];
      buffer[i+elementsize-j-1] = temp;
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    char newtype[16];
    int numbytes;
    int numelements;
    int elementsize;
    int new_numbytes;
    int new_numelements;
    int new_elementsize;
    mxClassID datatype;
    int dims[] = {1,1};
    
    if (nlhs!=1 || nrhs!=4 || !mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || mxGetNumberOfDimensions(prhs[0])>2 || (mxGetM(prhs[0])>1 && mxGetN(prhs[0])>1) || !mxIsChar(prhs[1]) || !mxIsLogicalScalar(prhs[2]) || !mxIsLogicalScalar(prhs[3])) {
        mexWarnMsgTxt("Improper input to TypeCastAndSwap -- first argument should be a real vector, second argument should be a string determining new type, third and fourth arguments should be boolean, true to swap bytes on input and output respectively.  One output argument.");
    }

    numelements = mxGetNumberOfElements(prhs[0]);
    elementsize = mxGetElementSize(prhs[0]);
    numbytes = numelements * elementsize;
    mxGetString(prhs[1],newtype,16);
    
    if (strcmp(newtype,"uint8")==0 || strcmp(newtype,"unsigned char")==0) {
        new_elementsize = 1;
        datatype = mxUINT8_CLASS;
    }
    else if (strcmp(newtype,"int8")==0 || strcmp(newtype,"signed char")==0) {
        new_elementsize = 1;
        datatype = mxINT8_CLASS;
    }
    else if (strcmp(newtype,"uint16")==0 || strcmp(newtype,"unsigned short")==0) {
        new_elementsize = 2;
        datatype = mxUINT16_CLASS;
    }
    else if (strcmp(newtype,"int16")==0 || strcmp(newtype,"signed short")==0 || strcmp(newtype,"short")==0) {
        new_elementsize = 2;
        datatype = mxINT16_CLASS;
    }
    else if (strcmp(newtype,"uint32")==0 || strcmp(newtype,"unsigned long")==0 || strcmp(newtype,"unsigned int")==0) {
        new_elementsize = 4;
        datatype = mxUINT32_CLASS;
    }
    else if (strcmp(newtype,"int32")==0 || (strcmp(newtype,"signed long")==0) || strcmp(newtype,"long")==0 || strcmp(newtype,"signed int")==0 || strcmp(newtype,"int")==0) {
        new_elementsize = 4;
        datatype = mxINT32_CLASS;
    }
    else if (strcmp(newtype,"uint64")==0) {
        new_elementsize = 8;
        datatype = mxUINT64_CLASS;
    }
    else if (strcmp(newtype,"int64")==0) {
        new_elementsize = 8;
        datatype = mxINT64_CLASS;
    }
    else if (strcmp(newtype,"single")==0 || strcmp(newtype,"float")==0) {
        new_elementsize = 4;
        datatype = mxSINGLE_CLASS;
    }
    else if (strcmp(newtype,"double")==0) {
        new_elementsize = 8;
        datatype = mxDOUBLE_CLASS;
    }
    else {
        new_elementsize = mxGetElementSize(prhs[0]);
        datatype = mxGetClassID(prhs[0]);
        if (strcmp(newtype,"")!=0)
            mexWarnMsgTxt("Invalid class identifier, keeping current class.  Valid classes are [u]int(8|16|32|64), [[un]signed] (char|short|int|long), float, single, and double");
    }
    
    new_numbytes = numbytes - (numbytes % new_elementsize);
    if (new_numbytes != numbytes) {
        mexWarnMsgTxt("Input data does not divide evenly into desired output format -- input padded with zeros.");
        new_numbytes += new_elementsize;
    }
    new_numelements = new_numbytes / new_elementsize;
    
    if (mxGetN(prhs[0])==1)
        dims[0] = new_numelements;
    else
        dims[1] = new_numelements;
    
    plhs[0] = mxCreateNumericArray(2,dims,datatype,mxREAL);
    memcpy(mxGetData(plhs[0]),mxGetData(prhs[0]),numbytes);
    
    if (mxIsLogicalScalarTrue(prhs[2])) {
        if (elementsize>1)
            FlipBytes(mxGetData(plhs[0]),elementsize,numelements);
        else
            mexWarnMsgTxt("Swapbytes_in was true for single-byte input");
    }
    
    if (mxIsLogicalScalarTrue(prhs[3]))  {
        if (new_elementsize>1)
            FlipBytes(mxGetData(plhs[0]),new_elementsize,new_numelements);
        else
            mexWarnMsgTxt("Swapbytes_out was true for a single-byte output");
    }
}

