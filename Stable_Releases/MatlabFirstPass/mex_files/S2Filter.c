/* [area_diff s2area max_s1area] = ...
 *     S2Filter(data, s2window, s1window, threshold)
 */

#include "mex.h"
#include <stdlib.h>
#include <math.h>

/* main function -- called by matlab */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    /*** first declare all variables (some C compilers require this) ***/
     
    /* inputs (data from prhs goes here) */
    double* data;
    int num_samples;
    int s1window;
    int s2window;
    double threshold;

    /* outputs (set plhs to point to these) */
    double* area_diff;
    double* s2area;
    double* max_s1area;

    /* indices (declare all indices used later on) */
    int i_rhs;
    int i_lhs;
    int i_smp;
    int i_smpw;

    /* other working variables */
    double s1area;
    double* local_s1areas;
    int num_local_s1s;
    int half_s2window;
    int half_s1window;
    
    
    /*** now parse prhs and plhs *///
     
    /* set default inputs */
    data = NULL;
    num_samples = 0;
    s2window = 200;
    s1window = 20;
    threshold = 0;
    
    /* load available inputs */
    i_rhs = 3;
    if (nrhs>i_rhs)
        if (!mxIsEmpty(prhs[i_rhs]))
            threshold = (double) mxGetScalar(prhs[i_rhs]);
    
    i_rhs = 2;
    if (nrhs>i_rhs)
        if (!mxIsEmpty(prhs[i_rhs]))
            s1window = (int) mxGetScalar(prhs[i_rhs]);
    
    i_rhs = 1;
    if (nrhs>i_rhs)
        if (!mxIsEmpty(prhs[i_rhs]))
            s2window = (int) mxGetScalar(prhs[i_rhs]);
    
    i_rhs = 0;
    if (nrhs>i_rhs)
        if (!mxIsEmpty(prhs[i_rhs])) {
            data = (double*) mxGetPr(prhs[i_rhs]);
            num_samples = mxGetNumberOfElements(prhs[i_rhs]);
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
    for (i_lhs=3;i_lhs<nlhs;i_lhs++)
        plhs[i_lhs] = mxCreateDoubleMatrix(0,0,mxREAL);

    i_lhs = 0;
    if (nlhs>i_lhs) {
        plhs[i_lhs] = mxCreateDoubleMatrix(num_samples,1,mxREAL);
        area_diff = (double*) mxGetPr(plhs[i_lhs]);
    }
    else
        area_diff = (double*) mxCalloc(num_samples,sizeof(double));
    
    i_lhs = 1;
    if (nlhs>i_lhs) {
        plhs[i_lhs] = mxCreateDoubleMatrix(num_samples,1,mxREAL);
        s2area = (double*) mxGetPr(plhs[i_lhs]);
    }
    else
        s2area = (double*) mxCalloc(num_samples,sizeof(double));
    
    i_lhs = 2;
    if (nlhs>i_lhs) {
        plhs[i_lhs] = mxCreateDoubleMatrix(num_samples,1,mxREAL);
        max_s1area = (double*) mxGetPr(plhs[i_lhs]);
    }
    else
        max_s1area = (double*) mxCalloc(num_samples,sizeof(double));
    

    /*** Now begin the real code ***/
    /* initialize some variables */
    num_local_s1s = s2window-s1window+1;
    local_s1areas = (double*) mxCalloc(num_local_s1s,sizeof(double));
    half_s2window = (s2window+1)/2;
    half_s1window = (s1window+1)/2;
    
    /**** initialize filter ****/
    /** initialize s2area **/
    s2area[0] = 0;
    for (i_smp=0; i_smp<half_s2window; i_smp++)
        s2area[0] += data[i_smp];
    for (i_smp=1; i_smp<half_s2window; i_smp++)
        s2area[i_smp] = s2area[i_smp-1] + data[i_smp+half_s2window-1];
    
    /** initialize max_s1area **/
    /* initialize local_s1areas */
    s1area = 0;
    for (i_smpw=0;i_smpw<half_s1window;i_smpw++)
        s1area += data[i_smpw];
    local_s1areas[0] = s1area;
    for (i_smpw=1; i_smpw<half_s1window; i_smpw++) {
        s1area += data[i_smpw+half_s1window-1];
        local_s1areas[i_smpw] = s1area;
    }
    for (i_smpw=half_s1window; i_smpw<num_local_s1s; i_smpw++) {
        s1area += (data[i_smpw+half_s1window-1] - data[i_smpw-half_s1window]);
        local_s1areas[i_smpw] = s1area;
    }

    /* initialize first half_s2window-half_s1window max_s1area */
    max_s1area[0] = local_s1areas[0];
    for (i_smpw=0; i_smpw<(half_s2window-half_s1window+1); i_smpw++)
        if (local_s1areas[i_smpw]>max_s1area[0])
            max_s1area[0] = local_s1areas[i_smpw];
    for (i_smp=1; i_smp<(half_s2window-half_s1window+1); i_smp++) {
        i_smpw = i_smp+half_s2window-half_s1window;
        if (local_s1areas[i_smpw]>max_s1area[i_smp-1])
            max_s1area[i_smp] = local_s1areas[i_smpw];
        else
            max_s1area[i_smp] = max_s1area[i_smp-1];
    }
    
    /* initialize rest of first half_s2window max_s1area */
    for (i_smp=(half_s2window-half_s1window+1); i_smp<half_s2window; i_smp++) {
        i_smpw = i_smp+half_s2window-half_s1window;
        s1area += (data[i_smp+half_s2window-1] - data[i_smp+half_s2window-1-s1window]);
        if (s1area>max_s1area[i_smp-1]) {
            max_s1area[i_smp] = s1area;
            local_s1areas[i_smpw%num_local_s1s] = s1area;
        }
        else if (max_s1area[i_smp-1] == local_s1areas[i_smpw%num_local_s1s]) {
            local_s1areas[i_smpw%num_local_s1s] = s1area;
            max_s1area[i_smp] = local_s1areas[0];
            for (i_smpw=1;i_smpw<num_local_s1s;i_smpw++)
                if (max_s1area[i_smp]<local_s1areas[i_smpw])
                    max_s1area[i_smp] = local_s1areas[i_smpw];
        }
        else {
            local_s1areas[i_smpw%num_local_s1s] = s1area;
            max_s1area[i_smp] = max_s1area[i_smp-1];
        }
    }
    
    /* initialize area_diff */
    for (i_smp = 0; i_smp<half_s2window; i_smp++)
        area_diff[i_smp] = s2area[i_smp] - max_s1area[i_smp];
    
    /*** pass filter over bulk of waveform ****/
    for (i_smp=half_s2window; i_smp<(num_samples-half_s2window+1); i_smp++) {
        i_smpw = i_smp+half_s2window-half_s1window;
        s2area[i_smp] = s2area[i_smp-1] + data[i_smp+half_s2window-1] - data[i_smp-half_s2window];
        s1area += (data[i_smp+half_s2window-1] - data[i_smp+half_s2window-1-s1window]);
        if (s1area>max_s1area[i_smp-1]) {
            max_s1area[i_smp] = s1area;
            local_s1areas[i_smpw%num_local_s1s] = s1area;
        }
        else if (max_s1area[i_smp-1] == local_s1areas[i_smpw%num_local_s1s]) {
            local_s1areas[i_smpw%num_local_s1s] = s1area;
            max_s1area[i_smp] = local_s1areas[0];
            for (i_smpw=1;i_smpw<num_local_s1s;i_smpw++)
                if (max_s1area[i_smp]<local_s1areas[i_smpw])
                    max_s1area[i_smp] = local_s1areas[i_smpw];
        }
        else {
            local_s1areas[i_smpw%num_local_s1s] = s1area;
            max_s1area[i_smp] = max_s1area[i_smp-1];
        }
        area_diff[i_smp] = s2area[i_smp] - max_s1area[i_smp];
    }
    
    /*** filter end of waveform ***/
    /* find max_s1area up to half_s1window from end */
    for (i_smp=(num_samples-half_s2window+1); i_smp<(num_samples-half_s1window+1); i_smp++) {
        i_smpw = i_smp+half_s2window-half_s1window;
        s1area += (data[i_smp+half_s2window-1] - data[i_smp+half_s2window-1-s1window]);
        if (s1area>max_s1area[i_smp-1]) {
            max_s1area[i_smp] = s1area;
            local_s1areas[i_smpw%num_local_s1s] = s1area;
        }
        else if (max_s1area[i_smp-1] == local_s1areas[i_smpw%num_local_s1s]) {
            local_s1areas[i_smpw%num_local_s1s] = s1area;
            max_s1area[i_smp] = local_s1areas[0];
            for (i_smpw=1;i_smpw<num_local_s1s;i_smpw++)
                if (max_s1area[i_smp]<local_s1areas[i_smpw])
                    max_s1area[i_smp] = local_s1areas[i_smpw];
        }
        else {
            local_s1areas[i_smpw%num_local_s1s] = s1area;
            max_s1area[i_smp] = max_s1area[i_smp-1];
        }
    }
    
    /* find max_s1area to the end */
    for (i_smp=(num_samples-half_s1window+1); i_smp<num_samples; i_smp++) {
        i_smpw = i_smp+half_s2window-half_s1window;
        s1area -= data[i_smp+half_s2window-1-s1window];
        if (s1area>max_s1area[i_smp-1]) {
            max_s1area[i_smp] = s1area;
            local_s1areas[i_smpw%num_local_s1s] = s1area;
        }
        else if (max_s1area[i_smp-1] == local_s1areas[i_smpw%num_local_s1s]) {
            local_s1areas[i_smpw%num_local_s1s] = s1area;
            max_s1area[i_smp] = local_s1areas[0];
            for (i_smpw=1;i_smpw<num_local_s1s;i_smpw++)
                if (max_s1area[i_smp]<local_s1areas[i_smpw])
                    max_s1area[i_smp] = local_s1areas[i_smpw];
        }
        else {
            local_s1areas[i_smpw%num_local_s1s] = s1area;
            max_s1area[i_smp] = max_s1area[i_smp-1];
        }
    }
    
    /* find s2area and area_diff to the end */
    for (i_smp=(num_samples-half_s2window+1); i_smp<num_samples; i_smp++) {
        s2area[i_smp] = s2area[i_smp-1] - data[i_smp-half_s2window];
        area_diff[i_smp] = s2area[i_smp] - max_s1area[i_smp];
    }
    
    
    /*** apply s2area threshold to area_diff ***/
    i_smp = 0;
    while (i_smp<num_samples) {
        if (s2area[i_smp] <= threshold) {
            /* create window of zeros in area_diff around this point */
            if (i_smp<half_s2window)
                for (i_smpw=0;i_smpw<(i_smp+half_s2window);i_smpw++)
                    area_diff[i_smpw] = 0;
            else if (i_smp>(num_samples-half_s2window))
                for (i_smpw=(i_smp-half_s2window+1);i_smpw<num_samples;i_smpw++)
                    area_diff[i_smpw] = 0;
            else
                for (i_smpw=(i_smp-half_s2window+1);i_smpw<(i_smp+half_s2window);i_smpw++)
                    area_diff[i_smpw] = 0;
            i_smp++;
            
            /* extend window as long as s2area is below threshold */
            while ((i_smp<=(num_samples-half_s2window)) && (s2area[i_smp] <= threshold)) {
                area_diff[i_smp+half_s2window-1] = 0;
                i_smp++;
            }
            
            /* if we've already gone to the end, no need to go further */
            if (i_smp>(num_samples-half_s2window))
                i_smp = num_samples;
        }
        else
            i_smp++;
    }



    /*** free dynamically allocated variables ***/
    mxFree(local_s1areas);
    if (nlhs<1)
        mxFree(area_diff);
    if (nlhs < 2)
        mxFree(s2area);
    if (nlhs < 3)
        mxFree(max_s1area);
}
