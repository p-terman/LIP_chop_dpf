#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* This c code does what MEM_fread.m was supposed to do, but was too slow
 * for multiple calls to effectively handle.  Have a blast with it. */

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

/* declare variables */
    mxArray* MEM_fid;
    mxArray* count;
    mxArray* precision;
    mxArray* gbl_MEM_FID_INFO;
    int i_lhs;
    int fid;
    char merror[16];
    bool valid_count;
    bool valid_fid;
    double numtoread_d[2];
    int numtoread[2];
    char input_precision[8];
    char output_precision[8];
    char precision_string[32];
    int precision_string_length;
    int first_asterisk;
    int first_equals;
    int bytes_available;
    int input_elementsize;
    int output_elementsize;
    int count_read;
    int bytes_read;
    mxClassID input_datatype;
    mxClassID output_datatype;
    unsigned char* data;
    unsigned char* source_data;
    unsigned char* output_data;
    int c;
    
/* get global variables */
    gbl_MEM_FID_INFO = (mxArray*) mexGetVariablePtr("global","gbl_MEM_FID_INFO");
    
    if (nlhs!=3) {
        for (i_lhs=0;i_lhs<nlhs;i_lhs++)
            plhs[i_lhs] = mxCreateDoubleMatrix(0,0,mxREAL);
        return;
    }
    
    switch (nrhs) {
/*        case 1:
            MEM_fid = prhs[0];
            break; */
        case 2:
            MEM_fid = prhs[0];
            count = prhs[1];
            precision = mxCreateString("uint8=>char");
            break;
        case 3:
            MEM_fid = prhs[0];
            count = prhs[1];
            precision = prhs[2];
            break;
/*        case 4:
            MEM_fid = prhs[0];
            count = prhs[1];
            precision = prhs[2];
            skip = prhs[3];
            break;
        case 5:
            MEM_fid = prhs[0];
            count = prhs[1];
            precision = prhs[2];
            skip = prhs[3];
            machineformat = prhs[4];
            break; */
        default:
            for (i_lhs=0;i_lhs<nlhs;i_lhs++)
                plhs[i_lhs] = mxCreateDoubleMatrix(0,0,mxREAL);
            return;
    }
    
/* Get fid */    
    if (mxIsNumeric(MEM_fid) && !mxIsEmpty(MEM_fid))
        fid = (unsigned int) mxGetScalar(MEM_fid);
    else {
        for (i_lhs=0;i_lhs<nlhs;i_lhs++)
            plhs[i_lhs] = mxCreateDoubleMatrix(0,0,mxREAL);
        return;
    }
    
/* check validity of count input */    
    valid_count = false;
    if (mxIsDouble(count) && !mxIsEmpty(count)) {
        numtoread_d[0] = mxGetScalar(count);
        if (mxGetNumberOfElements(count)==1)
            numtoread_d[1] = 1;
        else
            numtoread_d[1] = (mxGetPr(count))[1];
        numtoread[0] = (int) numtoread_d[0];
        numtoread[1] = (int) numtoread_d[1];
        valid_count = (numtoread[0] > 0 && numtoread[1] > 0 && !mxIsInf(numtoread_d[0]) && !mxIsInf(numtoread_d[1]));
    }
    
    if (!valid_count) {
        for (i_lhs=0;i_lhs<nlhs;i_lhs++)
            plhs[i_lhs] = mxCreateDoubleMatrix(0,0,mxREAL);
        return;
    }
    
/* check validity of fid */    
    valid_fid = false;
    if (mxIsStruct(gbl_MEM_FID_INFO))
        if (mxGetNumberOfElements(gbl_MEM_FID_INFO) >= fid)
            if (~mxIsEmpty(mxGetField(gbl_MEM_FID_INFO,fid-1,"maxblocksize"))) {
                mxGetString(mxGetField(gbl_MEM_FID_INFO,fid-1,"merror"),merror,15);
                if (strcmp(merror,"NOT_AT_EOF")!=0)
                    valid_fid = true;
            }
    
    if (!valid_fid) {
        for (i_lhs=0;i_lhs<nlhs;i_lhs++)
            plhs[i_lhs] = mxCreateDoubleMatrix(0,0,mxREAL);
        return;
    }
    
/* parse precision string */
    if (mxIsChar(precision)) {
        mxGetString(precision,precision_string,31);
        precision_string_length = strlen(precision_string);
        first_asterisk = strcspn(precision_string,"*");
        first_equals = strcspn(precision_string,"=");
        if (first_asterisk==0 && first_equals==precision_string_length) {
            strncpy(input_precision,&(precision_string[1]),precision_string_length);
            strncpy(output_precision,&(precision_string[1]),precision_string_length);
        }
        else if (first_asterisk==precision_string_length && first_equals>0 && first_equals<(precision_string_length-2)) {
            strncpy(input_precision,precision_string,first_equals);
            input_precision[first_equals] = '\0';
            strncpy(output_precision,&(precision_string[first_equals+2]),precision_string_length-first_equals-2);
            output_precision[precision_string_length-first_equals-2] = '\0';
        }
        else if (first_asterisk==precision_string_length && first_equals==precision_string_length) {
            strcpy(input_precision,precision_string);
            strcpy(output_precision,"double");
        }
        else {
            for (i_lhs=0;i_lhs<nlhs;i_lhs++)
                plhs[i_lhs] = mxCreateDoubleMatrix(0,0,mxREAL);
            return;
        }
    }
    else {
        for (i_lhs=0;i_lhs<nlhs;i_lhs++)
            plhs[i_lhs] = mxCreateDoubleMatrix(0,0,mxREAL);
        return;
    }

/* set input_datatypes */
    if (strcmp(input_precision,"uint8")==0 || strcmp(input_precision,"unsigned char")==0) {
        input_elementsize = 1;
        input_datatype = mxUINT8_CLASS;
    }
    else if (strcmp(input_precision,"int8")==0 || strcmp(input_precision,"signed char")==0) {
        input_elementsize = 1;
        input_datatype = mxINT8_CLASS;
    }
    else if (strcmp(input_precision,"uint16")==0 || strcmp(input_precision,"unsigned short")==0) {
        input_elementsize = 2;
        input_datatype = mxUINT16_CLASS;
    }
    else if (strcmp(input_precision,"int16")==0 || strcmp(input_precision,"signed short")==0 || strcmp(input_precision,"short")==0) {
        input_elementsize = 2;
        input_datatype = mxINT16_CLASS;
    }
    else if (strcmp(input_precision,"uint32")==0 || strcmp(input_precision,"unsigned long")==0 || strcmp(input_precision,"unsigned int")==0) {
        input_elementsize = 4;
        input_datatype = mxUINT32_CLASS;
    }
    else if (strcmp(input_precision,"int32")==0 || (strcmp(input_precision,"signed long")==0) || strcmp(input_precision,"long")==0 || strcmp(input_precision,"signed int")==0 || strcmp(input_precision,"int")==0) {
        input_elementsize = 4;
        input_datatype = mxINT32_CLASS;
    }
    else if (strcmp(input_precision,"single")==0 || strcmp(input_precision,"float")==0) {
        input_elementsize = 4;
        input_datatype = mxSINGLE_CLASS;
    }
    else if (strcmp(input_precision,"double")==0) {
        input_elementsize = 8;
        input_datatype = mxDOUBLE_CLASS;
    }
    else {
        for (i_lhs=0;i_lhs<nlhs;i_lhs++)
            plhs[i_lhs] = mxCreateDoubleMatrix(0,0,mxREAL);
        return;
    }

/* set output_datatypes */
    if (strcmp(output_precision,"uint8")==0 || strcmp(output_precision,"unsigned char")==0) {
        output_elementsize = 1;
        output_datatype = mxUINT8_CLASS;
    }
    else if (strcmp(output_precision,"int8")==0 || strcmp(output_precision,"signed char")==0) {
        output_elementsize = 1;
        output_datatype = mxINT8_CLASS;
    }
    else if (strcmp(output_precision,"uint16")==0 || strcmp(output_precision,"unsigned short")==0) {
        output_elementsize = 2;
        output_datatype = mxUINT16_CLASS;
    }
    else if (strcmp(output_precision,"int16")==0 || strcmp(output_precision,"signed short")==0 || strcmp(output_precision,"short")==0) {
        output_elementsize = 2;
        output_datatype = mxINT16_CLASS;
    }
    else if (strcmp(output_precision,"uint32")==0 || strcmp(output_precision,"unsigned long")==0 || strcmp(output_precision,"unsigned int")==0) {
        output_elementsize = 4;
        output_datatype = mxUINT32_CLASS;
    }
    else if (strcmp(output_precision,"int32")==0 || (strcmp(output_precision,"signed long")==0) || strcmp(output_precision,"long")==0 || strcmp(output_precision,"signed int")==0 || strcmp(output_precision,"int")==0) {
        output_elementsize = 4;
        output_datatype = mxINT32_CLASS;
    }
    else if (strcmp(output_precision,"single")==0 || strcmp(output_precision,"float")==0) {
        output_elementsize = 4;
        output_datatype = mxSINGLE_CLASS;
    }
    else if (strcmp(output_precision,"double")==0) {
        output_elementsize = 8;
        output_datatype = mxDOUBLE_CLASS;
    }
    else {
        for (i_lhs=0;i_lhs<nlhs;i_lhs++)
            plhs[i_lhs] = mxCreateDoubleMatrix(0,0,mxREAL);
        return;
    }

/* get data */
    bytes_available = (int) (mxGetScalar(mxGetField(gbl_MEM_FID_INFO,fid-1,"currentblocksize")) + mxGetScalar(mxGetField(gbl_MEM_FID_INFO,fid-1,"offset")) - mxGetScalar(mxGetField(gbl_MEM_FID_INFO,fid-1,"filepos")));
    if (bytes_available<(numtoread[0]*numtoread[1]*input_elementsize))
        count_read = (bytes_available - (bytes_available%input_elementsize));
    else
        count_read = numtoread[0]*numtoread[1];

    bytes_read = count_read * input_elementsize;
    plhs[1] = mxCreateDoubleScalar((double) count_read);
    plhs[2] = mxCreateDoubleScalar((double) bytes_read);
    
    if (numtoread[1]==1)
        numtoread[0] = count_read;
    else
        numtoread[1] = (int) ceil(((double) count_read) / ((double) numtoread[0]));
    plhs[0] = mxCreateNumericArray(2,numtoread,output_datatype,mxREAL);

    if (strcmp(input_precision,output_precision)==0)
        data = (unsigned char*) mxGetData(plhs[0]);
    else
        data = (unsigned char*) malloc((size_t) bytes_read);
    
    source_data = (unsigned char*) mxGetData(mxGetField(gbl_MEM_FID_INFO,fid-1,"data"));
    source_data = &(source_data[((int) (mxGetScalar(mxGetField(gbl_MEM_FID_INFO,fid-1,"filepos")))) - ((int) (mxGetScalar(mxGetField(gbl_MEM_FID_INFO,fid-1,"offset"))))]);
    memcpy(data,source_data,bytes_read);
    
    if (mxIsLogicalScalarTrue(mxGetField(gbl_MEM_FID_INFO,fid-1,"swapbytes")))
        FlipBytes(data,input_elementsize,count_read);
    
    if (strcmp(input_precision,output_precision)!=0) {
        output_data = (unsigned char*) mxGetData(plhs[0]);
        switch (input_datatype) {
            case mxUINT8_CLASS:
                switch (output_datatype) {
                    case mxUINT8_CLASS:
                        for (c=0;c<count_read;c++)
                            *((unsigned char*) &(output_data[c*output_elementsize])) = (unsigned char) *((unsigned char*) &(data[c*input_elementsize]));
                        break;
                    case mxINT8_CLASS:
                        for (c=0;c<count_read;c++)
                            *((signed char*) &(output_data[c*output_elementsize])) = (signed char) *((unsigned char*) &(data[c*input_elementsize]));
                        break;
                    case mxUINT16_CLASS:
                        for (c=0;c<count_read;c++)
                            *((unsigned short*) &(output_data[c*output_elementsize])) = (unsigned short) *((unsigned char*) &(data[c*input_elementsize]));
                        break;
                    case mxINT16_CLASS:
                        for (c=0;c<count_read;c++)
                            *((signed short*) &(output_data[c*output_elementsize])) = (signed short) *((unsigned char*) &(data[c*input_elementsize]));
                        break;
                    case mxUINT32_CLASS:
                        for (c=0;c<count_read;c++)
                            *((unsigned long*) &(output_data[c*output_elementsize])) = (unsigned long) *((unsigned char*) &(data[c*input_elementsize]));
                        break;
                    case mxINT32_CLASS:
                        for (c=0;c<count_read;c++)
                            *((signed long*) &(output_data[c*output_elementsize])) = (signed long) *((unsigned char*) &(data[c*input_elementsize]));
                        break;
                    case mxSINGLE_CLASS:
                        for (c=0;c<count_read;c++)
                            *((float*) &(output_data[c*output_elementsize])) = (float) *((unsigned char*) &(data[c*input_elementsize]));
                        break;
                    case mxDOUBLE_CLASS:
                        for (c=0;c<count_read;c++)
                            *((double*) &(output_data[c*output_elementsize])) = (double) *((unsigned char*) &(data[c*input_elementsize]));
                        break;
                }
                break;
            case mxINT8_CLASS:
                switch (output_datatype) {
                    case mxUINT8_CLASS:
                        for (c=0;c<count_read;c++)
                            *((unsigned char*) &(output_data[c*output_elementsize])) = (unsigned char) *((signed char*) &(data[c*input_elementsize]));
                        break;
                    case mxINT8_CLASS:
                        for (c=0;c<count_read;c++)
                            *((signed char*) &(output_data[c*output_elementsize])) = (signed char) *((signed char*) &(data[c*input_elementsize]));
                        break;
                    case mxUINT16_CLASS:
                        for (c=0;c<count_read;c++)
                            *((unsigned short*) &(output_data[c*output_elementsize])) = (unsigned short) *((signed char*) &(data[c*input_elementsize]));
                        break;
                    case mxINT16_CLASS:
                        for (c=0;c<count_read;c++)
                            *((signed short*) &(output_data[c*output_elementsize])) = (signed short) *((signed char*) &(data[c*input_elementsize]));
                        break;
                    case mxUINT32_CLASS:
                        for (c=0;c<count_read;c++)
                            *((unsigned long*) &(output_data[c*output_elementsize])) = (unsigned long) *((signed char*) &(data[c*input_elementsize]));
                        break;
                    case mxINT32_CLASS:
                        for (c=0;c<count_read;c++)
                            *((signed long*) &(output_data[c*output_elementsize])) = (signed long) *((signed char*) &(data[c*input_elementsize]));
                        break;
                    case mxSINGLE_CLASS:
                        for (c=0;c<count_read;c++)
                            *((float*) &(output_data[c*output_elementsize])) = (float) *((signed char*) &(data[c*input_elementsize]));
                        break;
                    case mxDOUBLE_CLASS:
                        for (c=0;c<count_read;c++)
                            *((double*) &(output_data[c*output_elementsize])) = (double) *((signed char*) &(data[c*input_elementsize]));
                        break;
                }
                break;
            case mxUINT16_CLASS:
                switch (output_datatype) {
                    case mxUINT8_CLASS:
                        for (c=0;c<count_read;c++)
                            *((unsigned char*) &(output_data[c*output_elementsize])) = (unsigned char) *((unsigned short*) &(data[c*input_elementsize]));
                        break;
                    case mxINT8_CLASS:
                        for (c=0;c<count_read;c++)
                            *((signed char*) &(output_data[c*output_elementsize])) = (signed char) *((unsigned short*) &(data[c*input_elementsize]));
                        break;
                    case mxUINT16_CLASS:
                        for (c=0;c<count_read;c++)
                            *((unsigned short*) &(output_data[c*output_elementsize])) = (unsigned short) *((unsigned short*) &(data[c*input_elementsize]));
                        break;
                    case mxINT16_CLASS:
                        for (c=0;c<count_read;c++)
                            *((signed short*) &(output_data[c*output_elementsize])) = (signed short) *((unsigned short*) &(data[c*input_elementsize]));
                        break;
                    case mxUINT32_CLASS:
                        for (c=0;c<count_read;c++)
                            *((unsigned long*) &(output_data[c*output_elementsize])) = (unsigned long) *((unsigned short*) &(data[c*input_elementsize]));
                        break;
                    case mxINT32_CLASS:
                        for (c=0;c<count_read;c++)
                            *((signed long*) &(output_data[c*output_elementsize])) = (signed long) *((unsigned short*) &(data[c*input_elementsize]));
                        break;
                    case mxSINGLE_CLASS:
                        for (c=0;c<count_read;c++)
                            *((float*) &(output_data[c*output_elementsize])) = (float) *((unsigned short*) &(data[c*input_elementsize]));
                        break;
                    case mxDOUBLE_CLASS:
                        for (c=0;c<count_read;c++)
                            *((double*) &(output_data[c*output_elementsize])) = (double) *((unsigned short*) &(data[c*input_elementsize]));
                        break;
                }
                break;
            case mxINT16_CLASS:
                switch (output_datatype) {
                    case mxUINT8_CLASS:
                        for (c=0;c<count_read;c++)
                            *((unsigned char*) &(output_data[c*output_elementsize])) = (unsigned char) *((signed short*) &(data[c*input_elementsize]));
                        break;
                    case mxINT8_CLASS:
                        for (c=0;c<count_read;c++)
                            *((signed char*) &(output_data[c*output_elementsize])) = (signed char) *((signed short*) &(data[c*input_elementsize]));
                        break;
                    case mxUINT16_CLASS:
                        for (c=0;c<count_read;c++)
                            *((unsigned short*) &(output_data[c*output_elementsize])) = (unsigned short) *((signed short*) &(data[c*input_elementsize]));
                        break;
                    case mxINT16_CLASS:
                        for (c=0;c<count_read;c++)
                            *((signed short*) &(output_data[c*output_elementsize])) = (signed short) *((signed short*) &(data[c*input_elementsize]));
                        break;
                    case mxUINT32_CLASS:
                        for (c=0;c<count_read;c++)
                            *((unsigned long*) &(output_data[c*output_elementsize])) = (unsigned long) *((signed short*) &(data[c*input_elementsize]));
                        break;
                    case mxINT32_CLASS:
                        for (c=0;c<count_read;c++)
                            *((signed long*) &(output_data[c*output_elementsize])) = (signed long) *((signed short*) &(data[c*input_elementsize]));
                        break;
                    case mxSINGLE_CLASS:
                        for (c=0;c<count_read;c++)
                            *((float*) &(output_data[c*output_elementsize])) = (float) *((signed short*) &(data[c*input_elementsize]));
                        break;
                    case mxDOUBLE_CLASS:
                        for (c=0;c<count_read;c++)
                            *((double*) &(output_data[c*output_elementsize])) = (double) *((signed short*) &(data[c*input_elementsize]));
                        break;
                }
                break;
            case mxUINT32_CLASS:
                switch (output_datatype) {
                    case mxUINT8_CLASS:
                        for (c=0;c<count_read;c++)
                            *((unsigned char*) &(output_data[c*output_elementsize])) = (unsigned char) *((unsigned long*) &(data[c*input_elementsize]));
                        break;
                    case mxINT8_CLASS:
                        for (c=0;c<count_read;c++)
                            *((signed char*) &(output_data[c*output_elementsize])) = (signed char) *((unsigned long*) &(data[c*input_elementsize]));
                        break;
                    case mxUINT16_CLASS:
                        for (c=0;c<count_read;c++)
                            *((unsigned short*) &(output_data[c*output_elementsize])) = (unsigned short) *((unsigned long*) &(data[c*input_elementsize]));
                        break;
                    case mxINT16_CLASS:
                        for (c=0;c<count_read;c++)
                            *((signed short*) &(output_data[c*output_elementsize])) = (signed short) *((unsigned long*) &(data[c*input_elementsize]));
                        break;
                    case mxUINT32_CLASS:
                        for (c=0;c<count_read;c++)
                            *((unsigned long*) &(output_data[c*output_elementsize])) = (unsigned long) *((unsigned long*) &(data[c*input_elementsize]));
                        break;
                    case mxINT32_CLASS:
                        for (c=0;c<count_read;c++)
                            *((signed long*) &(output_data[c*output_elementsize])) = (signed long) *((unsigned long*) &(data[c*input_elementsize]));
                        break;
                    case mxSINGLE_CLASS:
                        for (c=0;c<count_read;c++)
                            *((float*) &(output_data[c*output_elementsize])) = (float) *((unsigned long*) &(data[c*input_elementsize]));
                        break;
                    case mxDOUBLE_CLASS:
                        for (c=0;c<count_read;c++)
                            *((double*) &(output_data[c*output_elementsize])) = (double) *((unsigned long*) &(data[c*input_elementsize]));
                        break;
                }
                break;
            case mxINT32_CLASS:
                switch (output_datatype) {
                    case mxUINT8_CLASS:
                        for (c=0;c<count_read;c++)
                            *((unsigned char*) &(output_data[c*output_elementsize])) = (unsigned char) *((signed long*) &(data[c*input_elementsize]));
                        break;
                    case mxINT8_CLASS:
                        for (c=0;c<count_read;c++)
                            *((signed char*) &(output_data[c*output_elementsize])) = (signed char) *((signed long*) &(data[c*input_elementsize]));
                        break;
                    case mxUINT16_CLASS:
                        for (c=0;c<count_read;c++)
                            *((unsigned short*) &(output_data[c*output_elementsize])) = (unsigned short) *((signed long*) &(data[c*input_elementsize]));
                        break;
                    case mxINT16_CLASS:
                        for (c=0;c<count_read;c++)
                            *((signed short*) &(output_data[c*output_elementsize])) = (signed short) *((signed long*) &(data[c*input_elementsize]));
                        break;
                    case mxUINT32_CLASS:
                        for (c=0;c<count_read;c++)
                            *((unsigned long*) &(output_data[c*output_elementsize])) = (unsigned long) *((signed long*) &(data[c*input_elementsize]));
                        break;
                    case mxINT32_CLASS:
                        for (c=0;c<count_read;c++)
                            *((signed long*) &(output_data[c*output_elementsize])) = (signed long) *((signed long*) &(data[c*input_elementsize]));
                        break;
                    case mxSINGLE_CLASS:
                        for (c=0;c<count_read;c++)
                            *((float*) &(output_data[c*output_elementsize])) = (float) *((signed long*) &(data[c*input_elementsize]));
                        break;
                    case mxDOUBLE_CLASS:
                        for (c=0;c<count_read;c++)
                            *((double*) &(output_data[c*output_elementsize])) = (double) *((signed long*) &(data[c*input_elementsize]));
                        break;
                }
                break;
            case mxSINGLE_CLASS:
                switch (output_datatype) {
                    case mxUINT8_CLASS:
                        for (c=0;c<count_read;c++)
                            *((unsigned char*) &(output_data[c*output_elementsize])) = (unsigned char) *((float*) &(data[c*input_elementsize]));
                        break;
                    case mxINT8_CLASS:
                        for (c=0;c<count_read;c++)
                            *((signed char*) &(output_data[c*output_elementsize])) = (signed char) *((float*) &(data[c*input_elementsize]));
                        break;
                    case mxUINT16_CLASS:
                        for (c=0;c<count_read;c++)
                            *((unsigned short*) &(output_data[c*output_elementsize])) = (unsigned short) *((float*) &(data[c*input_elementsize]));
                        break;
                    case mxINT16_CLASS:
                        for (c=0;c<count_read;c++)
                            *((signed short*) &(output_data[c*output_elementsize])) = (signed short) *((float*) &(data[c*input_elementsize]));
                        break;
                    case mxUINT32_CLASS:
                        for (c=0;c<count_read;c++)
                            *((unsigned long*) &(output_data[c*output_elementsize])) = (unsigned long) *((float*) &(data[c*input_elementsize]));
                        break;
                    case mxINT32_CLASS:
                        for (c=0;c<count_read;c++)
                            *((signed long*) &(output_data[c*output_elementsize])) = (signed long) *((float*) &(data[c*input_elementsize]));
                        break;
                    case mxSINGLE_CLASS:
                        for (c=0;c<count_read;c++)
                            *((float*) &(output_data[c*output_elementsize])) = (float) *((float*) &(data[c*input_elementsize]));
                        break;
                    case mxDOUBLE_CLASS:
                        for (c=0;c<count_read;c++)
                            *((double*) &(output_data[c*output_elementsize])) = (double) *((float*) &(data[c*input_elementsize]));
                        break;
                }
                break;
            case mxDOUBLE_CLASS:
                switch (output_datatype) {
                    case mxUINT8_CLASS:
                        for (c=0;c<count_read;c++)
                            *((unsigned char*) &(output_data[c*output_elementsize])) = (unsigned char) *((double*) &(data[c*input_elementsize]));
                        break;
                    case mxINT8_CLASS:
                        for (c=0;c<count_read;c++)
                            *((signed char*) &(output_data[c*output_elementsize])) = (signed char) *((double*) &(data[c*input_elementsize]));
                        break;
                    case mxUINT16_CLASS:
                        for (c=0;c<count_read;c++)
                            *((unsigned short*) &(output_data[c*output_elementsize])) = (unsigned short) *((double*) &(data[c*input_elementsize]));
                        break;
                    case mxINT16_CLASS:
                        for (c=0;c<count_read;c++)
                            *((signed short*) &(output_data[c*output_elementsize])) = (signed short) *((double*) &(data[c*input_elementsize]));
                        break;
                    case mxUINT32_CLASS:
                        for (c=0;c<count_read;c++)
                            *((unsigned long*) &(output_data[c*output_elementsize])) = (unsigned long) *((double*) &(data[c*input_elementsize]));
                        break;
                    case mxINT32_CLASS:
                        for (c=0;c<count_read;c++)
                            *((signed long*) &(output_data[c*output_elementsize])) = (signed long) *((double*) &(data[c*input_elementsize]));
                        break;
                    case mxSINGLE_CLASS:
                        for (c=0;c<count_read;c++)
                            *((float*) &(output_data[c*output_elementsize])) = (float) *((double*) &(data[c*input_elementsize]));
                        break;
                    case mxDOUBLE_CLASS:
                        for (c=0;c<count_read;c++)
                            *((double*) &(output_data[c*output_elementsize])) = (double) *((double*) &(data[c*input_elementsize]));
                        break;
                }
                break;
        }
        free(data);
    }
    
    return;
}

