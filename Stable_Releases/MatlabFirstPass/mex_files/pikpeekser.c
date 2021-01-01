/*
mex routine to find features (peaks) in data
usage: [peeks] = pikpeeks(trc,thr);

trc is the summed data trace, a vector of length [nb_samples]
as such, pikpeeks is designed to operate on traces which have been summed across all channels

thr is the % of max that the routine will use to determine when a pulse ends
thr = 0.10 is recommended

pfs 060527
ldv 090220 - moved all variable definitions to the top of the program
*/

#include "mex.h"
#include "matrix.h"

void mexFunction(	int nlhs, 
					mxArray *plhs[],
					int nrhs,
					const mxArray *prhs[]
					 ) {
					 
	int pts	, ch, evt;
	const int nouts = 7; /* how many outputs in the output array */
	int ii, jj, kk, nn, nd;
	const int *d;
	double *_in;
	double *_out;
	const int i = 0;
	double *_thr;
	double thr;
	/* initialize requisite variables to zero */
	double m = 0;
	double m1 = 0;
	int c = 0;
	int c1 = 0;
    /* others */
	double thr10 = 0;
	double thr50 = 0;
	int timedill = 0;
	int timemin = 0;
    int timedillcap = 0;
	int cmin = 0;
	int c0l = 0;
	int c10l = 0;
	int c50l = 0;
    int tmp = 0;
	int timedilr = 0;
	int cmax = 0;
	int c0r = 0;
	int c10r = 0;
	int c50r = 0;
	int c10r_flag_keep_looking = 0;

    
    
    
	 /* Find the dimension of the data */
    d = mxGetDimensions(prhs[i]);
    nd = mxGetNumberOfDimensions(prhs[i]);
	pts = d[0];
	ch = d[1];
    evt = d[2];
	/* mexPrintf("pts = %d\n",d[0]); */
	/* mexPrintf("mex!\n"); */
	
	/* point to the input data */
    _in = mxGetPr(prhs[i]);
    _thr = mxGetPr(prhs[1]);
    thr = _thr[0];  /* threshold fraction of peakmax */
	/* mexPrintf("threshold = %f\n",thr); */
	
	/* Create an mxArray for the output */    
		plhs[0] = mxCreateDoubleMatrix( 1 , nouts , mxREAL ); /* rows , cols */

	/* Create an array for the output's data */
		_out = (double *) mxMalloc( nouts * sizeof(double) );
	
	/****************** find peaks ! ***********************/
     
	for (ii = 0; ii < pts; ii++) { /* find the maximum in the trc */
		if ( _in[ii] > m1 ) {
			m1 = _in[ii];
			c1 = ii;
		}
	}
	/* use temp variables m,c to compare subsequent maxima */
	m = m1; /* Max */
	c = c1; /* Index of max */
	thr10 = m1*thr; /* convert '%' input to an actual number */
	thr50 = m1*0.5; /* convert '%' input to an actual number */
	
	/* Go to left
	dt = Inf
	for index from time_max -> time_max-dt
	- 50% point is a span with [c+1]>50% and [c]<=50%, 
	      also reset dt = max(timedil*(time_max - c),100 ns)
	- 10% is first with [c]<=10%, can't stop until we find this 
	*/
	/* Currently at c */
	timedill = 6.; /* How many times 50% spread we will search */
	timemin = 15; /* Min time we search to for a pulse end */
    timedillcap = 50; /* max samples between c50l and c0l (assuming c10l is found) */
	cmin = 0; /* min time we will go to */
	c0l = 0;
	c10l = 0;
	c50l = 0;
	
	while ((c > 0) && ((c > cmin) || (c10l==0))) { /* find half-max */
			c--;
			/* mexPrintf("c = %d\n",c); */
			if ((_in[c+1]>thr50) && (_in[c]<=thr50)) {
				c50l=c; 
				/*mexPrintf("cmin = %d\n",cmin); */
				/*cmin = c1 - max(timemin,timedil*(c1-c));*/
				tmp = timedill*(c1-c);
                if (tmp > (c1-c+timedillcap))
                    tmp = c1-c+timedillcap;
				if (tmp > timemin) {
					cmin = c1 - tmp;
				}else {
					cmin = c1 - timemin;
				}
			}
			/* mexPrintf("c50l = %d\n",c50l); */
			if ((c10l==0) && (_in[c]<=thr10)) {
				c10l=c; 
			}
	}
	c0l = c; 
	
	c=c1; /* re-initialize c */
	/* Right-most point in peak */
	timedilr = 8.; /* How many times 50% spread we will search */
	cmax = pts-1; /* max time we will go to */
	c0r = pts-1;
	c10r = pts-1;
	c50r = pts-1;
	c10r_flag_keep_looking = 1; /* Set to 0 when below 10% -  possible condition  for finishing */
		while ((c < pts-1) && ((c < cmax) || (c10r_flag_keep_looking))) { /* find half-max */
			c++;
			if ((_in[c-1]>thr50) && (_in[c]<=thr50)) {
				c50r=c; 
				/*cmax = c1 + max(timemin,timedil*(c-c1));*/
				tmp = timedilr*(c-c1);
                if (tmp > (c-c1+timedillcap))
                    tmp = c-c1+timedillcap;
				if (tmp > timemin) {
					cmax = c1 + tmp;
				}else {
					cmax = c1 + timemin;
				}
			}
			
			/* If trace rises to 2 * thr then reset the 10% marker */
			if ((~c10r_flag_keep_looking) && (_in[c]>=2*thr10)) {
				c10r_flag_keep_looking = 1; 
			}
			if ((c10r_flag_keep_looking) && (_in[c]<=thr10)) {
				c10r=c; 
				c10r_flag_keep_looking = 0;
			}
		}
		c0r = c; /* Left most point in peak */
	
	/* must add 1 to each value bc C starts w 0, MATLAB w 1 */ 
	_out[0] = c0l + 1;
	_out[1] = c10l + 1;
	_out[2] = c50l + 1;
	_out[3] = c1 + 1;
	_out[4] = c50r + 1;
	_out[5] = c10r + 1;
	_out[6] = c0r + 1;
	
	
	/* Assign the data array to the output array */
		mxSetPr(plhs[0], _out);

}
