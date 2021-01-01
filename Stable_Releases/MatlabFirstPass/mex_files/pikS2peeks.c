/* [pulse_times pulse_areas pulse_diffs trace_diffs] = ...
 *     pikS2peeks(data, s2window, s1window, threshold, rel_area_thresh, x50_cutpoint, num_to_find)
 */

#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>


/* la la la */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    /* inputs */
    double* data;
    int num_samples;
    int s1window;
    int s2window;
    double threshold;
    double rel_area_thresh;
    double x50_cutpoint;
    int num_to_find;
    /* outputs */
    double* pulse_times;
    double* pulse_areas;
    double* pulse_diffs;
    double* trace_diffs;
    /* indices */
    int i_lhs;
    int i_smp;
    int i_smpd;
    int i_smpw;
    int i_t;
    int i_pulse;
    int ii_pulse;
    int n_workingpulse;
    /* other working variables */
    double s2area;
    double s1area;
    double max_s1area;
    double* local_s1areas;
    int num_local_s1s;
    double area_diff;
    bool in_pulse;
    bool even_in_flattop;
    double local_threshold;
    double round_fixer;
    double area50_threshold;
    int cutpoint;
    int buffer;
    int* pulse_order;
    int target_pulse;
    int source_pulse;
    
    /* set default inputs */
    num_to_find = 3;
    x50_cutpoint = 2.5;
    rel_area_thresh = 0;
    threshold = .5;
    s1window = 20;
    s2window = 200;
    data = NULL;
    num_samples = 0;
    
    /* load available inputs */
    if (nrhs>=7)
        if (!mxIsEmpty(prhs[6]))
            num_to_find = mxGetScalar(prhs[6]);
    
    if (nrhs>=6)
        if (!mxIsEmpty(prhs[5]))
            x50_cutpoint = mxGetScalar(prhs[5]);
    
    if (nrhs>=5)
        if (!mxIsEmpty(prhs[4]))
            rel_area_thresh = mxGetScalar(prhs[4]);
    
    if (nrhs>=4)
        if (!mxIsEmpty(prhs[3]))
            threshold = mxGetScalar(prhs[3]);
    
    if (nrhs>=3)
        if (!mxIsEmpty(prhs[2]))
            s1window = (int) mxGetScalar(prhs[2]);
    
    if (nrhs>=2)
        if (!mxIsEmpty(prhs[1]))
            s2window = (int) mxGetScalar(prhs[1]);
    
    if (nrhs >=1)
        if (!mxIsEmpty(prhs[0])) {
            data = (double*) mxGetPr(prhs[0]);
            num_samples = mxGetNumberOfElements(prhs[0]);
        }
    
    /* check validity of inputs */
    if (s1window%2==0)
        s1window++;
    if (s2window%2==0)
        s2window++;

    if ((s1window >= s2window) || (s2window >= num_samples)) {
        mexWarnMsgTxt("s1 window longer than s2 window or s2 window is longer than trace length in pikS2peaks.mex");
        for (i_lhs=0;i_lhs<nlhs;i_lhs++)
            plhs[i_lhs] = mxCreateDoubleMatrix(0,0,mxREAL);
        return;
    }
    
    /* set outputs */
    for (i_lhs=4;i_lhs<nlhs;i_lhs++)
        plhs[i_lhs] = mxCreateDoubleMatrix(0,0,mxREAL);

    if (nlhs>=1) {
        plhs[0] = mxCreateDoubleMatrix(7,num_to_find+1,mxREAL);
        pulse_times = (double*) mxGetPr(plhs[0]);
    }
    else
        pulse_times = (double*) calloc((num_to_find+1)*7,sizeof(double));
    
    if (nlhs>=2) {
        plhs[1] = mxCreateDoubleMatrix(1,num_to_find+1,mxREAL);
        pulse_areas = (double*) mxGetPr(plhs[1]);
    }
    else
        pulse_areas = (double*) calloc(num_to_find+1,sizeof(double));
    
    if (nlhs>=3) {
        plhs[2] = mxCreateDoubleMatrix(1,num_to_find+1,mxREAL);
        pulse_diffs = (double*) mxGetPr(plhs[2]);
    }
    else
        pulse_diffs = (double*) calloc(num_to_find+1,sizeof(double));
    
    if (nlhs>=4) {
        plhs[3] = mxCreateDoubleMatrix(num_samples,3,mxREAL);
        trace_diffs = (double*) mxGetPr(plhs[3]);
    }
    else
        trace_diffs = (double*) calloc(3*num_samples,sizeof(double));
    

    /* initialize some variables */
    n_workingpulse = 0;
    num_local_s1s = s2window-s1window+1;
    local_s1areas = (double*) calloc(num_local_s1s,sizeof(double));
    pulse_order = (int*) calloc(2*(num_to_find+1),sizeof(int));
    buffer = ((s2window+1)/2)-s1window;
    round_fixer = .01;
    for (i_smp=0;i_smp<(num_local_s1s);i_smp++)
        local_s1areas[i_smp] = 0;
    
    s2area = 0;
    s1area = 0;
    for (i_smp=0;i_smp<s2window;i_smp++)
        s2area += data[i_smp];
    for (i_smp=0;i_smp<s1window;i_smp++)
        s1area += data[i_smp];
    local_s1areas[0] = s1area;
    max_s1area = s1area;
    for (i_smp=s1window;i_smp<s2window;i_smp++) {
        s1area += (data[i_smp] - data[i_smp-s1window]);
        local_s1areas[i_smp-s1window+1] = s1area;
        if (s1area>max_s1area)
            max_s1area = s1area;
    }
    area_diff = s2area - max_s1area;
    trace_diffs[(s2window-1)/2] = area_diff;
    trace_diffs[((s2window-1)/2)+num_samples] = s2area;
    trace_diffs[((s2window-1)/2)+2*num_samples] = max_s1area;
    
    if (area_diff > threshold) {
        in_pulse = true;
        for (i_t=4;i_t<7;i_t++)
            pulse_times[n_workingpulse*7 + i_t] = 0;
        for (i_t=0;i_t<3;i_t++)
            pulse_times[n_workingpulse*7 + i_t] = i_t+1;        
        pulse_times[n_workingpulse*7 + 3] = (s2window+1)/2;
        pulse_areas[n_workingpulse] = s2area;
        pulse_diffs[n_workingpulse] = area_diff;
        local_threshold = area_diff*rel_area_thresh;
        area50_threshold = area_diff*.5*s2window/(s2window-s1window);
        even_in_flattop = false;
    }
       
    /* find S2 pulses */
    for (i_smp=(s2window+1)/2;i_smp<(num_samples-((s2window-1)/2));i_smp++) {
        i_smpd = i_smp + ((s2window-1)/2);
        s2area += (data[i_smpd] - data[i_smpd - s2window]);
        s1area += (data[i_smpd] - data[i_smpd - s1window]);
        if (s1area>max_s1area) {
            max_s1area = s1area;
            local_s1areas[(i_smpd-s1window+1)%num_local_s1s] = s1area;
        }
        else if (max_s1area == local_s1areas[(i_smpd-s1window+1)%num_local_s1s]) {
            local_s1areas[(i_smpd-s1window+1)%num_local_s1s] = s1area;
            max_s1area = local_s1areas[0];
            for (i_smpw=1;i_smpw<num_local_s1s;i_smpw++)
                if (max_s1area<local_s1areas[i_smpw])
                    max_s1area = local_s1areas[i_smpw];
        }
        else
            local_s1areas[(i_smpd-s1window+1)%num_local_s1s] = s1area;
        area_diff = s2area - max_s1area;
        
        trace_diffs[i_smp] = area_diff;
        trace_diffs[i_smp+num_samples] = s2area;
        trace_diffs[i_smp+2*num_samples] = max_s1area;
        
        if ((!in_pulse && (area_diff>threshold)) || (in_pulse && (pulse_times[n_workingpulse*7+4]==0)) || ((area_diff > threshold) && (s2area > local_threshold) && (i_smp<cutpoint))) {
            if (in_pulse) {
                if (area_diff > pulse_diffs[n_workingpulse]) {
                    pulse_times[n_workingpulse*7+3] = i_smp+1;
                    pulse_times[n_workingpulse*7+4] = 0;
                    pulse_times[n_workingpulse*7+5] = 0;
                    pulse_times[n_workingpulse*7+6] = 0;
                    pulse_diffs[n_workingpulse] = area_diff;
                    pulse_areas[n_workingpulse] = s2area;
                    local_threshold = area_diff*rel_area_thresh + round_fixer;
                    area50_threshold = area_diff*.5*s2window/(s2window-s1window);
                    even_in_flattop = false;
                }
                else if (area_diff == pulse_diffs[n_workingpulse]) {
                    pulse_times[n_workingpulse*7+4] = 0;
                    pulse_times[n_workingpulse*7+5] = 0;
                    pulse_times[n_workingpulse*7+6] = 0;
                    if (even_in_flattop) {
                        pulse_times[n_workingpulse*7+3]++;
                        even_in_flattop = false;
                    }
                    else
                        even_in_flattop = true;
                }
                else if ((s2area < area50_threshold) && (pulse_times[n_workingpulse*7+4]==0)) {
                    pulse_times[n_workingpulse*7+4] = i_smp;
                    cutpoint = pulse_times[n_workingpulse*7+3] + x50_cutpoint*(pulse_times[n_workingpulse*7+4]-pulse_times[n_workingpulse*7+3]);
                }
            }
            else {
                in_pulse = true;
                pulse_times[n_workingpulse*7] = i_smp - buffer + 1;
                pulse_times[n_workingpulse*7+1] = i_smp + 1;
                pulse_times[n_workingpulse*7+2] = i_smp + 1;
                pulse_times[n_workingpulse*7+3] = i_smp + 1;
                pulse_times[n_workingpulse*7+4] = 0;
                pulse_times[n_workingpulse*7+5] = 0;
                pulse_times[n_workingpulse*7+6] = 0;
                pulse_areas[n_workingpulse] = s2area;
                pulse_diffs[n_workingpulse] = area_diff;
                local_threshold = area_diff*rel_area_thresh + round_fixer;
                area50_threshold = area_diff*.5*s2window/(s2window-s1window);
                even_in_flattop = false;
            }
        }
        else if (in_pulse) {
            if (s2area<=local_threshold) {
                if ((pulse_times[n_workingpulse*7+6]==0) || (pulse_times[n_workingpulse*7+6]>(i_smp-buffer+1)))
                    pulse_times[n_workingpulse*7+6] = i_smp-buffer+1;
                if ((pulse_times[n_workingpulse*7+5]==0) || (pulse_times[n_workingpulse*7+5]>(i_smp-buffer)))
                    pulse_times[n_workingpulse*7+5] = i_smp-buffer;
            }
            else if (i_smp>=cutpoint) {
                if ((pulse_times[n_workingpulse*7+6]==0) || (pulse_times[n_workingpulse*7+6]>(i_smp+1)))
                    pulse_times[n_workingpulse*7+6] = i_smp+1;
                if ((pulse_times[n_workingpulse*7+5]==0) || (pulse_times[n_workingpulse*7+5]>i_smp))
                    pulse_times[n_workingpulse*7+5] = i_smp;
            }
            else if (pulse_times[n_workingpulse*7+5]==0) {
                pulse_times[n_workingpulse*7+5] = i_smp;
                pulse_times[n_workingpulse*7+6] = i_smp + buffer;
            }
            /* require correct order with respect to pulse peak... this can be a problem
             * when the look-back pulse-end is triggered by a negative-going swing */
            if (pulse_times[n_workingpulse*7+5] <= (pulse_times[n_workingpulse*7+3]+1))
                pulse_times[n_workingpulse*7+5] = pulse_times[n_workingpulse*7+3]+2;
            if (pulse_times[n_workingpulse*7+6] <= pulse_times[n_workingpulse*7+5])
                pulse_times[n_workingpulse*7+6] = pulse_times[n_workingpulse*7+5]+1;

            if ((pulse_times[n_workingpulse*7+6]!=0) && (i_smp>=(pulse_times[n_workingpulse*7+6] + buffer - 1))) {
                in_pulse = false;
                local_threshold = threshold;
                cutpoint = 0;
                n_workingpulse = 0;
                for (i_pulse=1;i_pulse<(num_to_find+1);i_pulse++)
                    if (pulse_diffs[i_pulse]<pulse_diffs[n_workingpulse])
                        n_workingpulse = i_pulse;
            }
        }
    }
    
    /* wrap up pulse end */
    if (in_pulse) {
        if (pulse_times[n_workingpulse*7+4]==0)
            pulse_times[n_workingpulse*7+4] = num_samples;
        if (pulse_times[n_workingpulse*7+5]==0) {
            pulse_times[n_workingpulse*7+5] = num_samples;
            pulse_times[n_workingpulse*7+6] = num_samples;
        }
    }
    
    /* find leading times on left hand sides */
    for (i_pulse=0;i_pulse<(num_to_find+1);i_pulse++) {
        if (pulse_diffs[i_pulse]>0) {
            area50_threshold = pulse_diffs[i_pulse]*.5*s2window/(s2window-s1window);
            for (i_smp=(pulse_times[i_pulse*7+3]-1);trace_diffs[i_smp+num_samples]>area50_threshold;i_smp--);
            pulse_times[i_pulse*7+2] = i_smp + 2;
            cutpoint = pulse_times[i_pulse*7+3] - x50_cutpoint*(pulse_times[i_pulse*7+3]-pulse_times[i_pulse*7+2]);
            if (cutpoint<pulse_times[i_pulse*7])
                cutpoint = pulse_times[i_pulse*7];
            for (i_smp=i_smp;i_smp>cutpoint;i_smp--)
                if (trace_diffs[i_smp+num_samples]<=(pulse_diffs[i_pulse]*rel_area_thresh+round_fixer)) {
                    if (pulse_times[i_pulse*7+1]<(i_smp+buffer+2))
                        pulse_times[i_pulse*7+1] = i_smp + buffer + 2;
                    pulse_times[i_pulse*7] = i_smp + buffer + 1;
                    break;
                }
            if (pulse_times[i_pulse*7]<(cutpoint+1))
                pulse_times[i_pulse*7] = cutpoint+1;
            if (pulse_times[i_pulse*7+1]<(cutpoint+2))
                pulse_times[i_pulse*7+1] = (cutpoint+2);
            
            /* require correct order with respect to pulse peak... this can be a problem
             * when the look-back pulse-end is triggered by a negative-going swing */
            if (pulse_times[i_pulse*7+1] >= (pulse_times[i_pulse*7+3]-1))
                pulse_times[i_pulse*7+1] = pulse_times[i_pulse*7+3]-2;
            if (pulse_times[i_pulse*7] >= pulse_times[i_pulse*7+1])
                pulse_times[i_pulse*7] = pulse_times[i_pulse*7+1]-1;

            if (pulse_times[i_pulse*7+2]<=pulse_times[i_pulse*7+1])
                pulse_times[i_pulse*7+2]=pulse_times[i_pulse*7+1]+1;
            if (pulse_times[i_pulse*7+4]>=pulse_times[i_pulse*7+5])
                pulse_times[i_pulse*7+4]=pulse_times[i_pulse*7+5]-1;
        }
    }
    
    /* sort pulses */
    pulse_order[0] = 0; /* first column is position of nth pulse found, second column is which pulse is nth biggest */
    for (i_pulse=1;i_pulse<(num_to_find+1);i_pulse++) {
        pulse_order[i_pulse] = i_pulse;
        for (ii_pulse=0;ii_pulse<i_pulse;ii_pulse++)
            if (pulse_diffs[i_pulse]>pulse_diffs[ii_pulse]) {
                pulse_order[ii_pulse]++;
                pulse_order[i_pulse]--;
            }
    }
    for (i_pulse=0;i_pulse<(num_to_find+1);i_pulse++)
        pulse_order[pulse_order[i_pulse] + num_to_find + 1] = i_pulse;
    
    source_pulse = num_to_find;
    target_pulse = pulse_order[source_pulse + num_to_find+1];
    if (source_pulse != target_pulse) {
        pulse_diffs[target_pulse] = pulse_diffs[source_pulse];
        pulse_areas[target_pulse] = pulse_areas[source_pulse];
        memcpy(&(pulse_times[target_pulse*7]),&(pulse_times[source_pulse*7]),7*sizeof(double));
        pulse_order[target_pulse] = pulse_order[source_pulse];
        pulse_order[pulse_order[target_pulse]+num_to_find+1] = target_pulse;
    }
    
    for (i_pulse=0;i_pulse<num_to_find-1;i_pulse++) {
        target_pulse = i_pulse;
        source_pulse = pulse_order[target_pulse+num_to_find+1];
        if (target_pulse != source_pulse) {
            pulse_diffs[num_to_find] = pulse_diffs[target_pulse];
            pulse_areas[num_to_find] = pulse_areas[target_pulse];
            memcpy(&(pulse_times[num_to_find*7]),&(pulse_times[target_pulse*7]),7*sizeof(double));
            pulse_diffs[target_pulse] = pulse_diffs[source_pulse];
            pulse_areas[target_pulse] = pulse_areas[source_pulse];
            memcpy(&(pulse_times[target_pulse*7]),&(pulse_times[source_pulse*7]),7*sizeof(double));
            pulse_diffs[source_pulse] = pulse_diffs[num_to_find];
            pulse_areas[source_pulse] = pulse_areas[num_to_find];
            memcpy(&(pulse_times[source_pulse*7]),&(pulse_times[num_to_find*7]),7*sizeof(double));
            
            pulse_order[source_pulse] = pulse_order[target_pulse];
            pulse_order[pulse_order[target_pulse]+num_to_find+1] = source_pulse;
            pulse_order[target_pulse]=target_pulse;
            pulse_order[target_pulse+num_to_find+1]=target_pulse;
        }
    }   
    
    /* free dynamically allocated variables */
    free(local_s1areas);
    free(pulse_order);
    if (nlhs<1)
        free(pulse_times);
    if (nlhs < 2)
        free(pulse_areas);
    if (nlhs < 3)
        free(pulse_diffs);
    if (nlhs<4)
        free(trace_diffs);

}