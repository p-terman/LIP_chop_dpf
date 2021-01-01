
/*
Development of Mexed version of LUX01LoadRQ1s
2009-05-08 JJC
*/

#include "math.h"
#include "mex.h"
#include "string.h"
#include <stdio.h>
#include <stdlib.h>

/* use flip bytes because daquiri is in little endian, but gsk-10 is in big endian */
void FlipBytes(void* input_buffer, int elementsize, int count) { /* lovely code from eric dahl. call after fread if switching endianness */
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



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
int idebug = 0;


mxArray *mxdata ;
double *data_holder ;
mxArray **mxdata_holder ;
mxArray **mxdata_holder2 ;
mxArray *mxnl ;
int nl ;
mxArray *mxswitch_endian ;
int switch_endian ;
mxArray *mxvarsize_sum ;
int varsize_sum ;
mxArray *mxbegin_read ;
int begin_read ;
mxArray *mxlinewords ;
int linewords ;
mxArray *mxmax_nl ;
int max_nl ;
mxArray *mxfilename ;
mxArray *mxvartype ;
mxArray *mxentrytype ;
int *entrytype ;
mxArray *mxvarsize ;
mxArray *mxvartype1 ;
mxArray *mxvarsize1 ;
mxArray *mxname1 ;
char **vartype ;
int line_offset ;
int *varsize ;
FILE * fp_rq1 ;
char *filename ;
char **name ;
char **dname ;
mxArray *mxname ;
mxArray *mxvarbytes ;
mxArray *varbytes1 ;
int *varbytes ;
int nw ;
int ii ;
int dims[2] ;
int dims2[3] ;
int event_offset ;
int nli ;

if (idebug>=1) {mexPrintf("Entered Mex Function\n");}
/* Get all variables needed from the caller function (LUX01LoadRQ1s.m) */
mxdata = mexGetVariable("caller", "data_mex");
mxnl = mexGetVariable("caller", "nl");
nl = (int)mxGetPr(mxnl)[0];
mxvarsize_sum = mexGetVariable("caller", "varsize_sum");
if (idebug>=2) {mexPrintf("varsize_sum = %d\n", (int)mxGetPr(mxvarsize_sum)[0]);}
varsize_sum = (int)mxGetPr(mxvarsize_sum)[0];
if (idebug>=2) {mexPrintf("varsize_sum = %d\n", varsize_sum);}
mxbegin_read = mexGetVariable("caller", "begin_read");
begin_read = (int)mxGetPr(mxbegin_read)[0];
if (idebug>=2) {mexPrintf("begin_read = %d\n", begin_read);}
mxlinewords = mexGetVariable("caller", "linewords");
linewords = (int)mxGetPr(mxlinewords)[0];
if (idebug>=2) {mexPrintf("linewords = %d\n", linewords);}
mxmax_nl = mexGetVariable("caller", "max_nl");
max_nl = (int)mxGetPr(mxmax_nl)[0];
if (idebug>=2) {mexPrintf("max_nl = %d\n", max_nl);}
mxnl = mexGetVariable("caller", "nl");
nl = (int)mxGetPr(mxnl)[0];
if (idebug>=2) {mexPrintf("nl = %d\n", nl);}
mxswitch_endian = mexGetVariable("caller", "switch_endian");
switch_endian = (int)mxGetPr(mxswitch_endian)[0];
if (idebug>=2) {mexPrintf("switch_endian = %d\n", switch_endian);}
mxfilename = mexGetVariable("caller", "filename");
filename = mxCalloc((mxGetN(mxfilename)+1),sizeof(char));
mxGetString(mxfilename,filename,(mxGetN(mxfilename)+1));
if (idebug>=2) {mexPrintf("filename = %s\n", filename);}
mxentrytype = mexGetVariable("caller", "entrytype");
entrytype = (int*)mxGetPr(mxentrytype);

mxvartype = mexGetVariable("caller","vartype");
mxvarsize = mexGetVariable("caller","varsize");
mxname = mexGetVariable("caller","name");
mxvarbytes = mexGetVariable("caller","varbytes");
/* Get all variables needed from the caller function (LUX01LoadRQ1s.m) */


varbytes = mxCalloc((mxGetN(mxvarbytes)+1),sizeof(int));

varsize = mxCalloc(linewords,sizeof(int));
vartype = mxCalloc(linewords,sizeof(*vartype));
name = mxCalloc(linewords,sizeof(*name));
dname = mxCalloc(linewords,sizeof(*dname));


dims[0] = 1;
dims[1] = 20;

mxdata_holder = mxCalloc(linewords,sizeof(*mxdata_holder));
mxdata_holder2 = mxCalloc(linewords,sizeof(*mxdata_holder2));

for (nw=0; nw<linewords; nw++)
{
	mxname1 = mxGetCell(mxname, nw);
	name[nw] = mxArrayToString(mxname1);
	mxvartype1 = mxGetCell(mxvartype, nw);
	vartype[nw] = mxArrayToString(mxvartype1);
	mxvarsize1 = mxGetCell(mxvarsize, nw);
	varbytes[nw] = (int)mxGetPr(mxvarbytes)[nw];
	varsize[nw]=1;
	for (ii=0; ii<mxGetN(mxvarsize1); ii++)
	{varsize[nw] = varsize[nw]*((int)mxGetPr(mxvarsize1)[ii]);}
	if (idebug>=2) {mexPrintf("[%d]: %s    %s    %d    %d\n",nw, name[nw], vartype[nw], varsize[nw], varbytes[nw]);}
	
	dims[1] = varsize[nw];
	
	mxdata_holder[nw] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS,mxREAL);
}
dims2[2]=(max_nl-nl);
dims2[0]=1;
fp_rq1 = fopen(filename, "rb"); /* open file to read */
if (fp_rq1 == NULL)
{mexPrintf("file open failed!!\n\n");
return;}
if (idebug>=2) {mexPrintf("file opened\n");}
for (nw=0; nw<linewords; nw++)
{
	
	if (idebug>=2) {mexPrintf("dims2[0] %d\n",nw);}
	dims2[1]=varsize[nw];
	if (idebug>=2) {mexPrintf("dims2[0] %d\n",nw);}
	
	mxdata_holder2[nw] = mxCreateNumericArray(3,dims2,mxDOUBLE_CLASS,mxREAL);
	if (idebug>=2) {mexPrintf("dims2[0] %d\n",nw);}
}
if (idebug>=2) {mexPrintf("data holder created\n");}

/*
nl = max_nl-1;
linewords=1;
*/
nli=0;
while (nl < max_nl) {

	nl = nl+1;
	
	line_offset = varsize_sum * (nl-1) ;

	
	fseek(fp_rq1, begin_read+line_offset, SEEK_SET); 
   
	

	for (nw=0; nw<linewords; nw++) 
	{

		
	

		fread(mxGetData(mxdata_holder[nw]), varbytes[nw], varsize[nw], fp_rq1);

        if (switch_endian==1)
        {if (idebug>=2) {mexPrintf("switching endianness in mex file... hold on to your hats!\n");}
         FlipBytes(mxGetData(mxdata_holder[nw]), varbytes[nw], varsize[nw]);}
		if (idebug>=4) {mexPrintf("---- %d ----\n", nw);}
		
		event_offset = nli*varsize[nw];
	
		for (ii=0; ii<varsize[nw]; ii++) {
           
			 (double)mxGetPr(mxdata_holder2[nw])[ii+event_offset] = (double)mxGetPr(mxdata_holder[nw])[ii] ;
		}
		
		
		
	}

/*if (idebug>=3) {mexPrintf("==== %d ====\n", nl);}
*/




nli=nli+1;
}



for (nw=0; nw<linewords; nw++)
{
	mxSetFieldByNumber(mxdata, 0, nw, mxdata_holder2[nw]);
}








fclose(fp_rq1);
/*
mxFree(*vartype);
mxFree(filename);
mxFree(*mxdata_holder);
mxFree(*mxdata_holder2);
*/

mexPutVariable("caller","data_mex",mxdata);
/*mxFree(mxdata);*/
if (idebug>=1) {mexPrintf("Leaving Mex Function\n");}

}






/*  A bit of fun testing out the method. */ /*
#include "math.h"
#include "mex.h"
#include "string.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{


const mxArray * barb ;
const mxArray * tarr ;


barb = mexGetVariable("caller","varb");
if (barb==NULL)
  mexErrMsgTxt("Variable 'varb' not in workspace.");

tarr = mexGetVariable("caller","garr");
if (tarr==NULL)
  mexErrMsgTxt("Variable 'garr' not in workspace.");
  


 
mexPrintf("barb = %d\n", (int)mxGetPr(barb)[0]);
mexPrintf("tarr = %d\n", (int)mxGetPr(tarr)[0]);

mxGetPr(barb)[0] = mxGetPr(barb)[0]*2;

mexPutVariable("caller", "varb", barb);

mexPrintf("hello?\n");
}
*/














