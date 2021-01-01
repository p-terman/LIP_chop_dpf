#ifdef _WIN32
#define _WIN32_USEDLL_NOT
#endif
/* This code does not seem to work with dynamic linking -- use static
 * libraries.  If you want to try with dll's in windows, delete _NOT from
 * line two.  */

/* scan down past DLL junk for documentation */

#ifdef _WIN32_USEDLL
#define BZ_IMPORT
#endif

#include "mex.h"
#include "zlib.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* code for windows only */

#ifdef _WIN32_USEDLL
#include <io.h>
#endif


#ifdef _WIN32_USEDLL

#define BZ2_LIBNAME "libbz2-1.0.2.DLL" 

#include <windows.h>
static int BZ2DLLLoaded = 0;
static HINSTANCE BZ2DLLhLib;
int BZ2DLLLoadLibrary(void)
{
   HINSTANCE hLib;

   if(BZ2DLLLoaded==1){return 0;}
   hLib=LoadLibrary(BZ2_LIBNAME);
   if(hLib == NULL){
      fprintf(stderr,"Can't load %s\n",BZ2_LIBNAME);
      return -1;
   }
   BZ2_bzReadOpen= GetProcAddress(hLib,"BZ2_bzReadOpen");
   BZ2_bzReadClose= GetProcAddress(hLib,"BZ2_bzReadClose");
   BZ2_bzRead= GetProcAddress(hLib,"BZ2_bzRead");

   if (!BZ2_bzReadOpen || !BZ2_bzReadClose || !BZ2_bzRead ) {
      fprintf(stderr,"GetProcAddress failed.\n");
      return -1;
   }
   BZ2DLLLoaded=1;
   BZ2DLLhLib=hLib;
   return 0;

}
int BZ2DLLFreeLibrary(void)
{
   if(BZ2DLLLoaded==0){return 0;}
   FreeLibrary(BZ2DLLhLib);
   BZ2DLLLoaded=0;
   return 0;
}
#endif /* WIN32 */


/* This c code is intended for use as a mex file to open, navigate, read, and
 * close .bz2 files.  This is done by passing an fid array back and forth between
 * this code and Matlab every time this function is called.  These elements are stored
 * as uint32 (unsigned long), and are as follows:
 * gzfid[0] - pointer to the FILE object
 * gzfid[1] - pointer to the BZFILE object
 * gzfid[2] - latest bzerror code (should be converted to int32 to read,
 *            codes for different values are listed in bzlib.h
 * gzfid[3] - current position in file (first byte is zero)
 * gzfid[4] - max buffer to read when scanning through file -- used for setting
 *            position only, does not limit the size of buffer read out in 'read'
 *            mode.
 *
 * The code should be used as follows:
 * To open file:
 *    gzfid = GZRead('open', filename, uint32(max_buffer_size))
 * To set position in file:
 *    gzfid = GZRead('setpos', gzfid, uint32(new_position), [suppress_eof_warning])
 *                          (suppress_eof_warning is optional bool argument, default is false)
 * To read compressed data from current position:
 *    [gzfid data] = GZRead('read', gzfid, uint32(num_bytes_to_read), 'datatype', swapbytes)
 * To close file
 *    gzfid = GZRead('close', gzfid)
 *
 * When reading data, the data is returned as a column-vector of desired datatype.
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

bool mexGZFile_open(mxArray*** p_plhs, const mxArray *prhs[]) {
    /* This function opens the desired file, creating the FILE and BZFILE objects. 
     * The filename is stored in prhs[1], and the pointers to the FILE and BZFILE
     * objects created are stored in plhs[0] (along with current position, which is
     * zero) to be passed to this function in future calls.  */

    /* declare variables */
    char filename[256];
    /* FILE* rawfile; */
    gzFile* gzippedfile;
    /* int* gzerror; */
    unsigned long gzfid[5];
    const int dims[] = {1,5};
    
    /* open file */
    mxGetString(prhs[1],filename,256);
    /* rawfile = fopen(filename,"rb");
    if (rawfile==NULL) {
        mexWarnMsgTxt("cold not open file in GZRead, no action taken");
        return false;
    } */

    /* open file as gzipped stream */
    /* gzerror = (int*) mxMalloc(sizeof(int));
    *gzerror = 0; */
    gzippedfile = malloc(sizeof(gzFile));
    /* mexMakeMemoryPersistent(gzippedfile); */
    gzippedfile = gzopen(filename,"rb");

    /* save file info for future calls */
    gzfid[0] = (unsigned long) 0;
    gzfid[1] = (unsigned long) gzippedfile;
    gzfid[2] = (unsigned long) 0;
    gzfid[3] = 0;
    gzfid[4] = *((unsigned long*) mxGetData(prhs[2]));
    (*p_plhs)[0] = mxCreateNumericArray(2,dims,mxUINT32_CLASS,mxREAL);
    memcpy(mxGetData((*p_plhs)[0]),gzfid,5*sizeof(unsigned long));
    
    /* free memory */
    /* mxFree(bzerror); */
    return true;
}

bool mexGZFile_setpos(mxArray*** p_plhs, const mxArray *prhs[], bool suppress_eof_warning) {
    /* This function reads through the BZFILE to the desired point.  Pointers to
     * the relevant FILE and BZFILE instances are stored in prhs[1], along with the
     * current file position, and the new file position is specified in prhs[2].
     * plhs[0] returns pointers to the FILE and BZFILE instances (which will change,
     * if the file must be rewound) along with the new file position, and plhs[1] returns
     * the requested data as a uint8 vector.*/
    
    /* declare variables */
    unsigned long* gzfid;
    /* FILE* rawfile; */
    gzFile* gzippedfile;
    /* int* bzerror; */
    unsigned long currentpos;
    unsigned long maxbuffer_size;
    unsigned long newpos;
    unsigned char* dummy_buffer;
    unsigned long bytes_read;
    const int dims[] = {1,5};
    
    /* get file info */
    /* bzerror = (int*) mxMalloc(sizeof(int)); */
    gzfid = (unsigned long*) mxGetData(prhs[1]);
    /* rawfile = (FILE*) gzfid[0]; */
    gzippedfile = (gzFile*) gzfid[1];
    /* *bzerror = (int) gzfid[2]; */
    currentpos = gzfid[3];
    maxbuffer_size = gzfid[4];
    
    /* get desired position */
    newpos = *((unsigned long*) mxGetData(prhs[2]));
    
    /* gzseek(gzippedfile,newpos,0);
    currentpos = (unsigned long) gztell(gzippedfile);*/ /* gzseek and gztell don't seem to work */
    
    /* rewind file if necessary */
    if (currentpos > newpos) {
        gzrewind(gzippedfile);
        currentpos = 0;
    }

    /* read through file to desired point, in maximum steps of maxbuffer */
    dummy_buffer = (unsigned char*) mxMalloc(maxbuffer_size*sizeof(unsigned char));
    
    while ((newpos - currentpos) >= maxbuffer_size) {
        bytes_read = gzread(gzippedfile, dummy_buffer, maxbuffer_size);
        currentpos = currentpos + bytes_read;
        if (bytes_read != maxbuffer_size) {
            if (!suppress_eof_warning)
                mexWarnMsgTxt("failed to reach desired position in bzippedfile, could be end of file.");
            newpos = currentpos;
        }
    }
    
    if (newpos > currentpos) {
        bytes_read = gzread(gzippedfile, dummy_buffer, newpos-currentpos);
        currentpos = currentpos + bytes_read;
        if (newpos != currentpos)
            if (!suppress_eof_warning)
                mexWarnMsgTxt("failed to reach desired position in bzippedfile, could be end of file.");
    }
    
    /* save file info for future calls */
    gzfid = (unsigned long*) mxMalloc(5*sizeof(unsigned long));
    gzfid[0] = 0;
    gzfid[1] = (unsigned long) gzippedfile;
    gzfid[2] = 0;
    gzfid[3] = currentpos;
    gzfid[4] = maxbuffer_size;
    (*p_plhs)[0] = mxCreateNumericArray(2,dims,mxUINT32_CLASS,mxREAL);
    memcpy(mxGetData((*p_plhs)[0]),gzfid,5*sizeof(unsigned long));

    /* free memory */
    /* mxFree(bzerror); */
    mxFree(gzfid);
    mxFree(dummy_buffer);
    return true;
}

bool mexGZFile_read(mxArray*** p_plhs, const mxArray *prhs[]) {
    /* This function reads the desired number of bytes from the BZFILE.  Pointers
     * to the relevant FILE and BZFILE are stored in prhs[1], along with the current
     * file position, and the number of desired bytes is specified in prhs[2].
     * plhs[0] returns pointers to the FILE and BZFILE instances along with the new
     * file position.  */

    /* declare variables */
    unsigned long* gzfid;
    /* FILE* rawfile; */
    gzFile* gzippedfile;
    /* int* bzerror; */
    unsigned long currentpos;
    unsigned long maxbuffer_size;
    char newtype[16];
    unsigned long numbytes;
    unsigned long numelements;
    unsigned long elementsize;
    unsigned long bytes_read;
    mxClassID datatype;
    const int dims[] = {1,5};
    int bdims[] = {1,1};
    
    /* get file info */
    /* bzerror = (int*) mxMalloc(sizeof(int)); */
    gzfid = (unsigned long*) mxGetData(prhs[1]);
    /* rawfile = (FILE*) gzfid[0]; */
    gzippedfile = (gzFile*) gzfid[1];
    /* *bzerror = (int) gzfid[2]; */
    currentpos = gzfid[3];
    maxbuffer_size = gzfid[4];
    
    /* get numbytes to read */
    numelements = *((unsigned long*) mxGetData(prhs[2]));
    mxGetString(prhs[3],newtype,16);
    if (strcmp(newtype,"uint8")==0 || strcmp(newtype,"unsigned char")==0) {
        elementsize = 1;
        datatype = mxUINT8_CLASS;
    }
    else if (strcmp(newtype,"int8")==0 || strcmp(newtype,"signed char")==0) {
        elementsize = 1;
        datatype = mxINT8_CLASS;
    }
    else if (strcmp(newtype,"uint16")==0 || strcmp(newtype,"unsigned short")==0) {
        elementsize = 2;
        datatype = mxUINT16_CLASS;
    }
    else if (strcmp(newtype,"int16")==0 || strcmp(newtype,"signed short")==0 || strcmp(newtype,"short")==0) {
        elementsize = 2;
        datatype = mxINT16_CLASS;
    }
    else if (strcmp(newtype,"uint32")==0 || strcmp(newtype,"unsigned long")==0 || strcmp(newtype,"unsigned int")==0) {
        elementsize = 4;
        datatype = mxUINT32_CLASS;
    }
    else if (strcmp(newtype,"int32")==0 || (strcmp(newtype,"signed long")==0) || strcmp(newtype,"long")==0 || strcmp(newtype,"signed int")==0 || strcmp(newtype,"int")==0) {
        elementsize = 4;
        datatype = mxINT32_CLASS;
    }
    else if (strcmp(newtype,"uint64")==0) {
        elementsize = 8;
        datatype = mxUINT64_CLASS;
    }
    else if (strcmp(newtype,"int64")==0) {
        elementsize = 8;
        datatype = mxINT64_CLASS;
    }
    else if (strcmp(newtype,"single")==0 || strcmp(newtype,"float")==0) {
        elementsize = 4;
        datatype = mxSINGLE_CLASS;
    }
    else if (strcmp(newtype,"double")==0) {
        elementsize = 8;
        datatype = mxDOUBLE_CLASS;
    }
    else {
        elementsize = 1;
        datatype = mxUINT8_CLASS;
        if (strcmp(newtype,"")!=0)
            mexWarnMsgTxt("Invalid class identifier, using uint8.  Valid classes are [u]int(8|16|32|64), [[un]signed] (char|short|int|long), float, single, and double");
    }
    
    numbytes = elementsize * numelements;
    
    /* read buffer */
    bdims[0] = numelements;
    (*p_plhs)[1] = mxCreateNumericArray(2,bdims,datatype,mxREAL);
    bytes_read = gzread(gzippedfile, mxGetData((*p_plhs)[1]), numbytes);
    currentpos = currentpos + bytes_read;
    /*if (bytes_read != numbytes)
        mexWarnMsgTxt("failed to reach desired position in bzippedfile, could be end of file.");
    */
    /* flip bytes if desired */
    if (mxIsLogicalScalarTrue(prhs[4]) && elementsize>1)
        FlipBytes(mxGetData((*p_plhs)[1]),elementsize,numelements);
    
    /* save file info for future calls */
    gzfid = (unsigned long*) mxMalloc(5*sizeof(unsigned long));
    gzfid[0] = 0;
    gzfid[1] = (unsigned long) gzippedfile;
    gzfid[2] = 0;
    gzfid[3] = currentpos;
    gzfid[4] = maxbuffer_size;
    (*p_plhs)[0] = mxCreateNumericArray(2,dims,mxUINT32_CLASS,mxREAL);
    memcpy(mxGetData((*p_plhs)[0]),gzfid,5*sizeof(unsigned long));

    /* free memory */
    /* mxFree(bzerror); */
    mxFree(gzfid);
    return true;
}

bool mexGZFile_close(mxArray*** p_plhs, const mxArray *prhs[]) {
    /* This function closes the BZFILE and FILE, and returns zeros in plhs[0].  */
    
    /* declare variables */
    unsigned long* gzfid;
    /* FILE* rawfile; */
    gzFile* gzippedfile;
    /* int* bzerror; */
    unsigned long currentpos;
    unsigned long maxbuffer_size;
    const int dims[] = {1,5};
    
    /* get file info */
    /* bzerror = (int*) mxMalloc(sizeof(int)); */
    gzfid = (unsigned long*) mxGetData(prhs[1]);
    /* rawfile = (FILE*) gzfid[0]; */
    gzippedfile = (gzFile*) gzfid[1];
    /* *bzerror = (int) gzfid[2]; */
    currentpos = gzfid[3];
    maxbuffer_size = gzfid[4];

    /* close BZFILE stream */
    gzclose(gzippedfile);
    /* free(gzippedfile); */
    /* apparently gzclose also frees gzippedfile */
    
    /* close FILE */
    /* fclose(rawfile); */
    
    /* save file info for future calls */
    gzfid = (unsigned long*) mxMalloc(5*sizeof(unsigned long));
    gzfid[0] = 0;
    gzfid[1] = (unsigned long) gzippedfile;
    gzfid[2] = 0;
    gzfid[3] = currentpos;
    gzfid[4] = maxbuffer_size;
    (*p_plhs)[0] = mxCreateNumericArray(2,dims,mxUINT32_CLASS,mxREAL);
    memcpy(mxGetData((*p_plhs)[0]),gzfid,5*sizeof(unsigned long));

    /* free memory */
    /* mxFree(bzerror); */
    mxFree(gzfid);
#ifdef _WIN32_USEDLL    
    BZ2DLLFreeLibrary();
#endif
    return true;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    /* This function checks the validity of the inputs, and calls the proper function
     * to perform the desired operation.  */

    /* declare variables */
    char op_mode[8];
    double errnum = -1;
    int n=0;
    bool suppress_eof_warning;
    
    /* windows dll use only (does not work -- use static libraries) */
#ifdef _WIN32_USEDLL
    if (BZ2DLLLoaded == 0) {
        if(BZ2DLLLoadLibrary()<0) {
            mexWarnMsgTxt("Loading of dll failed.  Giving up.\n");
            return;
        }
    }
#endif    

    /* check first argument is string */
    if (nrhs<1 || !mxIsChar(prhs[0]))
        mexWarnMsgTxt("Requires string for first argument -- valid values are 'open', 'setpos', 'read', and 'close', no action taken");
    else {

        /* check proper number of return arguments */
        mxGetString(prhs[0],op_mode,7);
        if (!((nlhs==1 && strcmp(op_mode,"read")!=0) || (nlhs==2 && strcmp(op_mode,"read")==0)))
            mexWarnMsgTxt("GZRead must return exactly one or two arguments, namely the information on the currently opened file, and the requested data in 'read' mode, no action taken");
        else {
            if (strcmp(op_mode,"open")==0) {
                if (nrhs<3 || !mxIsChar(prhs[1]) || !mxIsUint32(prhs[2]) || mxGetNumberOfElements(prhs[2])<1)
                    mexWarnMsgTxt("'open' mode requires second argument to be filename (string), and third argument to be maxbuffer_size (uint32), no action taken");
                else if (mexGZFile_open(&plhs,prhs))
                    return;
            }
            else if (strcmp(op_mode,"setpos")==0) {
                if (nrhs<4)
                    suppress_eof_warning = false;
                else
                    suppress_eof_warning = mxIsLogicalScalarTrue(prhs[3]);
                if (nrhs<3 || !mxIsUint32(prhs[1]) || mxGetNumberOfElements(prhs[1])<5 || !mxIsUint32(prhs[2]) || mxGetNumberOfElements(prhs[2])<1)
                    mexWarnMsgTxt("'setpos' mode requires second argument to be bz_fileID array, third argument to be file_position (uint32), and fourth argument to be boolean (print warnings or not) -- no action taken");
                else if (mexGZFile_setpos(&plhs,prhs,suppress_eof_warning))
                    return;
            }
            else if (strcmp(op_mode,"read")==0) {
                if (nrhs<5 || !mxIsUint32(prhs[1]) || mxGetNumberOfElements(prhs[1])<5 || !mxIsUint32(prhs[2]) || mxGetNumberOfElements(prhs[2])<1 || !mxIsChar(prhs[3]) || !mxIsLogicalScalar(prhs[4]))
                    mexWarnMsgTxt("'read' mode requires second argument to be bz_fileID array, third argument to be number of elements to be read (uint32), fourth argument to be desired datatype, and fifth argument to be boolean for switching endianness -- no action taken");
                else if (mexGZFile_read(&plhs,prhs))
                    return;
            }
            else if (strcmp(op_mode,"close")==0) {
                if (nrhs<2 || !mxIsUint32(prhs[1]) || mxGetNumberOfElements(prhs[1])<5)
                    mexWarnMsgTxt("'close' mode requires second argument to be bz_fileID array, no action taken");
                else if (mexGZFile_close(&plhs,prhs))
                    return;
            }
            else
                mexWarnMsgTxt("invalid first argument -- valid values are 'open', 'setpos', 'read', and 'close'");
        }
    }
    /* if the code reaches this far, there was an error.  
     * Set first return argument to second input argument (if it exists) and the rest to -1 */
    if (nrhs>1) {
        plhs[0] = mxDuplicateArray(prhs[1]);
        n++;
    }

    for (n=n;n<nlhs;n++) {
        plhs[n] = mxCreateDoubleMatrix(1,1,mxREAL);
        *mxGetPr(plhs[n])=errnum;
    }
}

