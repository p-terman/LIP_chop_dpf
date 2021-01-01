/* perusePeeks.c

all documentation is contained in perusePeeks.m

121211 pfs
140206 ss - removed timing calculations now in separate module
140620 ss - modification to split possible multiple scatters 

*/

#include "mex.h"
#include "matrix.h"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	/* variable declarations */					 
	int ii=0;
	int jj=0;
	int ch=0;
	int c2=0;
	int tt=0;
	int nSamplesToSlideFilter;
	double *eventSamples, *tmpSamples;
	double *eventRecord, *tmpRecord;
	double *parameters;
	double *fullBoxAreasRolling;
	double *fullBoxCsum;
    double *average;
	double m2=0;
	double tmpBoxArea=0; /* instantaneous, incremental value */
	const int *d;
	/* outputs: */
	double *afTiming;     int dims0 = 2;
        
	double *fBoxArea;     int dims2 = 3;
        
	
	
	 /* get the dimension of the input data */
		d = mxGetDimensions(prhs[0]);
		
	/* point to the input data */
		tmpSamples = mxGetPr(prhs[0]);
		tmpRecord = mxGetPr(prhs[1]);
		parameters = mxGetPr(prhs[2]);

	/* tmpSamples is a pod concatenation and may contain gaps, which we must fill in */
		int gapSamples;
		int recordLength = tmpSamples[d[0]-1]-tmpSamples[0] + 1 ;
		/* mexPrintf("\n%1.0f - %1.0f = %1.0f\n",tmpSamples[d[0]-1],tmpSamples[0],recordLength); */
		
		eventSamples = mxCalloc(recordLength,sizeof(double));
		eventRecord  = mxCalloc(recordLength,sizeof(double));

		for (ii=0;ii<recordLength;ii++) { /* initialize */
			eventSamples[ii] = tmpSamples[0]+ii;
			eventRecord[ii] = 0;
		}

		int nn=0;
		for (ii=0;ii<(d[0]-1);ii++) {
			gapSamples = tmpSamples[ii+1]-tmpSamples[ii];
			if (gapSamples==1) {
				eventRecord[nn] = tmpRecord[ii];
				/*mexPrintf("%5.3f\n",eventRecord[nn]);*/
				nn=nn+1;
			}
			else if (gapSamples>1) {
				eventRecord[nn] = tmpRecord[ii]; /* first point */
				/*mexPrintf("%5.3f\n",eventRecord[nn]);*/
				nn=nn+1;
				for (jj=0;jj<(gapSamples-1);jj++) { /* fill in zeros */
					eventRecord[nn] = 0; 
					/*mexPrintf("%5.3f\n",eventRecord[nn]);*/
					nn=nn+1;
				}
			}
			else {
				mexPrintf("ii=%d jj=%d nn=%d gapSamples=%d\n",ii,jj,nn,gapSamples);
				mexPrintf("major issue with indexing, bailing out (perusePeeks.c)."); return;
			}
		}
		
		
		int preBoxSamples        = (int)parameters[0]; /* better be an integer number of samples! */
		int fullBoxSamples   	 = (int)parameters[1]; 
		int postBoxSamples       = (int)parameters[2];
		double edgeFraction  	 = parameters[3];
		double txFraction  	     = parameters[4];

		int skinnyBoxSamples     = (int)parameters[5];
		int maximumGap           = parameters[6]; 
		int nLookAhead           = (int)parameters[7];
		int nLookBehind          = (int)parameters[8];
		double noiseThreshold    = parameters[9];
/*
                mexPrintf("Pulse finding in prod version with following input parameters:\n");
                mexPrintf("  preBoxSamples %d\n", preBoxSamples);
                mexPrintf("  fullBoxSamples %d\n", fullBoxSamples);
                mexPrintf("  postBoxSamples %d\n", postBoxSamples);
                mexPrintf("  edgeFraction %0.3f\n", edgeFraction);
                mexPrintf("  skinnyBoxSamples %d\n", skinnyBoxSamples);
                mexPrintf("  maximumGap %d\n", maximumGap);
                mexPrintf("  nLookAhead %d\n", nLookAhead);
                mexPrintf("  nLookBehind %d\n", nLookBehind);
                mexPrintf("  noiseThre %0.3f\n", noiseThreshold);
*/               
	/* allocate some memory */
		afTiming    = mxCalloc( dims0,sizeof(double));
		fBoxArea    = mxCalloc( dims2,sizeof(double));
        
		
		fullBoxAreasRolling = mxCalloc(recordLength,sizeof(double));
		fullBoxCsum         = mxCalloc(recordLength,sizeof(double));
        average             = mxCalloc(recordLength,sizeof(double));
    
		
	/* initialize arrays to zero (is this necessary, or paranoia? -- necessary 120703) */
		for (ii=0;ii<dims0;ii++) {
			afTiming[ii]=0;
		}

		for (ii=0;ii<dims2;ii++) {
			fBoxArea[ii]=0;
		}
        for (ii=0;ii<recordLength;ii++) {
        average[ii]=0;
        }


	/* all set, algorithm follows -------------------------------------------------- */
		
	/* error trapping */
		if (fullBoxSamples>recordLength) {
			/* mexPrintf("%d > %d\n",fullBoxSamples,recordLength); */
			/*flags[0] = 1;*/
			fullBoxSamples = recordLength; /* so we don't run off the end of the record */
			nSamplesToSlideFilter = 1;	
		}
		else {
			nSamplesToSlideFilter = (recordLength-fullBoxSamples-1);	
		}
		 /*	mexPrintf("recordLength = [%d %d]\n",d[0],d[1]); */
		/* 	mexPrintf("fullBoxSamples = %d\n",fullBoxSamples); */

	/* run the fullBoxFilter for each starting time tt */
		for (tt=0; tt<nSamplesToSlideFilter; tt++) {
			fullBoxAreasRolling[tt]=0; tmpBoxArea=0;
			for (ii=tt; ii<(tt+fullBoxSamples); ii++) { /* recall indexing column-wise */
				tmpBoxArea = tmpBoxArea + eventRecord[ii];
			}	
			fullBoxAreasRolling[tt] = tmpBoxArea;
		}	
		/* which fullBoxAreasRolling is largest? */
		m2=0; c2=0;
		for (tt=0; tt<nSamplesToSlideFilter; tt++) {
			if (fullBoxAreasRolling[tt]>=m2) {
				m2 = fullBoxAreasRolling[tt];
				c2 = tt;
			}
			/*mexPrintf("tt %d :: c2 %d :: m2 %2.1f (%2.1f)\n",tt,c2,m2,fullBoxAreasRolling[tt]);  */
		}


	if (fullBoxSamples>0){		
		
        double maxValue=0;
        int maxIndex=0;
        for (ii=c2;ii<(c2+fullBoxSamples);ii++) {
            if (eventRecord[ii] > maxValue) {
                maxValue = eventRecord[ii];
                maxIndex = ii;
                
            }
        }
        /* mexPrintf("maxValue = %3.2f\n",maxValue);mexPrintf("maxIndex = %d\n",maxIndex);*/

/* ********* perusePeeks2 mods */
			/* Should we modify the box start and box width? to determine this, */
			/* check for pulse end using a maximum gap cut on lightly filtered signal. */  
			int t0_max_gap = 0;
			int n_gaps = 0;
			for (ii=maxIndex;ii>=0;ii--) {
				/* look forward and backwards wrt scan dir and calculate average */
				int n_sum=0;
				double sum=0.;
				int ii_sum;
				for(ii_sum=(ii+nLookBehind); ii_sum>=(ii-nLookAhead); ii_sum--){
					if(ii_sum>recordLength || ii_sum<0) continue;   
					sum+=eventRecord[ii_sum];
					n_sum+=1;
				}
				sum = sum/n_sum;

				/* count sequential samples with average below baseline threshold */
				if(sum < noiseThreshold) n_gaps += 1; 
				else n_gaps = 0;

				/* flag as end of pulse if have quiet period greater than maximumGap */
				/* or we run out of samples during a quiet period */
				if (n_gaps>maximumGap || (n_gaps > 0 && ii == (0+nLookBehind))) {
					t0_max_gap = ii+n_gaps-1-nLookBehind; /* pulse start defined as start of quiet period */
					break;
				}
				if (t0_max_gap<0){ /* do not run off the eventRecord */
					t0_max_gap=0;
				}
				/* mexPrintf("%d,%1.0f,%d\n",ii,eventSamples[ii],n_gaps); */
			}

			/* as before but now scan forwards in time along the falling edge of pulse */ 
			n_gaps = 0;
			int t2_max_gap = recordLength-1;
			for (ii=maxIndex;ii<recordLength;ii++) {
				int n_sum=0;
				double sum=0.;
				int ii_sum;
				for(ii_sum=(ii-nLookBehind); ii_sum<=(ii+nLookAhead); ii_sum++){
					if(ii_sum>recordLength || ii_sum<0) continue;   
					sum+=eventRecord[ii_sum];
					n_sum+=1;
				}
				sum = sum/n_sum;

				if(sum < noiseThreshold) n_gaps += 1; 
				else n_gaps = 0;

				if (n_gaps > maximumGap || (n_gaps > 0 && ii == (recordLength-1-nLookAhead))){
					t2_max_gap = ii-n_gaps+1+nLookAhead; /* pulse start defined as start of quiet period */
					break;
				}
				/* mexPrintf("%d,%1.0f,%d\n",ii,eventSamples[ii],n_gaps); */

			}
	
			/* variable box width definition*/
			/* now use interval found using maximum gap method as pulse region */
			fullBoxSamples = (t2_max_gap-t0_max_gap+1);
			c2 = t0_max_gap;
/*			if (fullBoxSamples<skinnyBoxSamples) {
				fullBoxSamples=skinnyBoxSamples; 
			} */
			/*	mexPrintf("fullBoxSamples %d ... c2 %d\n",fullBoxSamples,c2); */

                        /* double check we haven't run out of waveform */
                        int maxwidth = (recordLength-1)-(c2-1);
                        if(fullBoxSamples > maxwidth){
                            fullBoxSamples = maxwidth;
                        }
			
        
        
        
/* ********* perusePeeks2 mods */
	
	
		/* assign afTiming */ 
			afTiming[0] = c2; /* index at which max area full box begins */
			afTiming[1] = (c2+fullBoxSamples-1); /* box ends */
			/* mexPrintf("%1.1f\n",afTiming[0]); */
			/* mexPrintf("%1.1f\n",afTiming[6]); */
			/* mexPrintf("%d\n",recordLength); */
        
        /* experimental - is pulse a double scatter? first look right of peak */
        
        double maxValue2=0;
        int maxIndex2=0;
        for (ii=c2;ii<(c2+fullBoxSamples);ii++) {
            if (eventRecord[ii] > maxValue2) {
                maxValue2 = eventRecord[ii];
                maxIndex2 = ii;
            }
        }
        
        int flag_10r = 0;
        int flag_10l = 0;
        int flag_10rr = 0;
        int flag_10ll = 0;
        int flag_20r = 0;
        int flag_20l = 0;
        int index_10r = 0;
        int index_10l = 0;
        int index_20r = 0;
        int index_20l = 0;
        int min_l;
        int min_r;
        int index_min_l = 0;
        int index_min_r = 0;
        int len_10r = 0;
        int len_10l = 0;
        int len_l = 0;
        int len_r = 0;
        int fflagl = 0;
        int fflagr = 0;
        
        
        if (fullBoxSamples > 600) {
            int ll;
            for (ll=maxIndex2;ll<(c2+fullBoxSamples-(nLookAhead+24));ll++) {
                double avg = 0;
                int n_avg = 0;
                int ll_sum;
                for(ll_sum=(ll-(nLookBehind+24)); ll_sum<=(ll+(nLookAhead+24)); ll_sum++){
					if(ll_sum>recordLength || ll_sum<0) continue;
					avg+=eventRecord[ll_sum];
					n_avg+=1;
				}
                average[ll] = avg/n_avg;
                
                
                if (average[ll] < 0.4*eventRecord[maxIndex2] && flag_10r == 0) {
                    flag_10r = 1;  /* record pulse has fallen */
                    index_10r = ll;
                }
            
                if (average[ll] > 0.1*eventRecord[maxIndex2] && average[ll] > average[ll-1] && flag_10r == 1 && flag_20r == 0) {
                    flag_20r = 1; /* record pulse has rerisen */
                    index_20r = ll;
                }
                
                if (flag_20r==1 && average[ll] >= 0.1*eventRecord[maxIndex2] && average[ll] > average[ll-1]) {
                    len_r += 1; /* how long is potential pulse rising */
                }
                if (flag_20r ==1 && average[ll] < 0.1*eventRecord[maxIndex2]) break;
            }
                if (len_r > 50) {
                    min_r = eventRecord[index_10r];
                    int mm;
                    for (mm=index_10r; mm<=index_20r; mm++) {
                        if (eventRecord[mm] <= min_r && eventRecord[mm]>0) {
                            min_r = eventRecord[mm];
                            index_min_r = mm;
                            fflagr = 1; /* to prevent afTiming being set to zero if not satisfied */
                        }
                    }
                    if (fflagr == 1) {
                    afTiming[1] = index_min_r; /* reset end point to lowest between two pulses */
                    }
                }
            
			for(ll=maxIndex2; ll>c2+(nLookBehind+24); ll--) { /* now look left of peak */
                double avg = 0;
                int n_avg = 0;
                int ll_sum;
                for(ll_sum=(ll-(nLookBehind+24)); ll_sum<=(ll+(nLookAhead+24)); ll_sum++){
					if(ll_sum>recordLength || ll_sum<0) continue;
					avg+=eventRecord[ll_sum];
					n_avg+=1;
				}
                average[ll] = avg/n_avg;
                if (average[ll] < 0.4*eventRecord[maxIndex2] && flag_10l ==0) {
                    flag_10l = 1;  /* record pulse has fallen */
                    index_10l = ll;
                }
                if (average[ll] > 0.1*eventRecord[maxIndex2] && average[ll] > average[ll+1] && flag_10l == 1 && flag_20l == 0) {
                    flag_20l = 1; /* record pulse has rerisen */
                    index_20l = ll;
                }
                
                if (flag_20l==1 && average[ll] >= 0.1*eventRecord[maxIndex2] && average[ll] > average[ll+1]) {
                    len_l += 1; /* how long is potential pulse rising */
                }
                if (flag_20l ==1 && average[ll] < 0.1*eventRecord[maxIndex2]) break; /* break if pulse starts to fall */
            
            }
            if (len_l > 50) {
                min_l = eventRecord[index_10l];
                int mm;
                for (mm=index_10l; mm>index_20l; mm--) {
                    if (eventRecord[mm] <= min_l && eventRecord[mm]>0) {
                        min_l = eventRecord[mm];
                        index_min_l = mm;
                        fflagl = 1; /* to prevent afTiming being set to zero if not satisfied */
                    }
                    
                }
                if (fflagl == 1) {
                    afTiming[0] = index_min_l; /* reset start point to lowest between two pulses */
                }
            }
        }
        
        
	
			/* need the cumulative sum */
				fullBoxCsum[(int)afTiming[0]] = eventRecord[(int)afTiming[0]]; /* the first point */
				for (tt=afTiming[0]; tt<=afTiming[1]-1; tt++) { /* -1 here, +1 line below */
					fullBoxCsum[tt+1] = fullBoxCsum[tt] + eventRecord[tt+1];
					/* mexPrintf("%d ::%4.1f :: %4.1f\n",tt+1,eventRecord[tt+1],fullBoxCsum[tt+1]); */
				}		
			/* assign boxArea. but first, get the pre- and post- box "wing" areas */
	
			/* pre-fullBox */
			tmpBoxArea=0;
			for (ii=afTiming[0]-1; ii>=(afTiming[0]-preBoxSamples); ii--) {
				if (ii<0) { /* error-trap */
					break; /* flag= ?? */
				}
				tmpBoxArea = tmpBoxArea + eventRecord[ii];
			}
			fBoxArea[0] = tmpBoxArea;
	
			/* fullBox */
			tmpBoxArea=0;
			for (ii=afTiming[0]; ii<=(afTiming[1]); ii++) {
				if (ii>=recordLength) { /* error-trap */
					break; /* flag= ?? */
				}
				tmpBoxArea = tmpBoxArea + eventRecord[ii];
			}
			fBoxArea[1] = tmpBoxArea; 
	
			/* post-fullBox */
			tmpBoxArea=0;
			for (ii=afTiming[1]+1; ii<=(afTiming[1]+postBoxSamples); ii++) {
				if (ii>=recordLength) { /* error-trap */
					break; /* flag= ?? */
				}
				tmpBoxArea = tmpBoxArea + eventRecord[ii];
				/* mexPrintf("%d :: %4.1f\n",ii,eventRecord[ii]); */
			}
			fBoxArea[2] = tmpBoxArea;
		
	
	
		/* ONCE EVERYTHING ELSE IS DONE. return sample value, not sample index */
			for (ii=0;ii<dims0;ii++) { 
				tt = (int)afTiming[ii];
				afTiming[ii]=eventSamples[tt]; /* don't want the sample index, but the actual sample (time) value */
				/* mexPrintf("%d %1.0f\n",tt,eventSamples[tt]); */
			}	
			/* mexPrintf("afTiming: [%1.0f %1.0f %1.0f %1.0f %1.0f]\n",afTiming[0],afTiming[1],afTiming[2],afTiming[3],afTiming[4]); */
			/* mexPrintf("fBoxArea: [%1.1f %1.1f %1.1f]\n",fBoxArea[0],fBoxArea[1],fBoxArea[2]); */
	
	}	
	/* Create an mxArray for the output */    
		/*  now using dynamic memory allocation 121110 pfs */
		plhs[0] = mxCreateNumericMatrix(1,dims0,mxDOUBLE_CLASS,mxREAL);
		plhs[1] = mxCreateNumericMatrix(1,dims2,mxDOUBLE_CLASS,mxREAL);


	/* Assign the data array to the output array */
		mxSetPr(plhs[0], afTiming);
		mxSetPr(plhs[1], fBoxArea);
	

		
}

