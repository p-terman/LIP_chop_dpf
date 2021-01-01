% LUXFirstPass_CreateRQStruct
%
% DO NOT CALL THIS MANUALLY, IT DOESN'T MAKE SENSE
%
% To be called from LUXFirstPass
% Assign variables to dp structure before writing binary file
%
% 2010-02-02 JJC, DCM
%

%clear dp variable before assigning!
clear dp

%dp{1}.set = filename_prefix;
dp{1}.t0 = t0+pretrigger;
dp{1}.t10l = t10l+pretrigger;
dp{1}.t50l = t50l+pretrigger;
dp{1}.t1 = t1+pretrigger;
dp{1}.t50r = t50r+pretrigger;
dp{1}.t10r = t10r+pretrigger;
dp{1}.t2 = t2+pretrigger;
dp{1}.ft10l = ft10l+pretrigger;
dp{1}.ft50l = ft50l+pretrigger;
dp{1}.ft1 = ft1+pretrigger;
dp{1}.ft50r = ft50r+pretrigger;
dp{1}.ft10r = ft10r+pretrigger;

dp{1}.peak_height = peak_height;
dp{1}.whichpeak = whichpeak;
dp{1}.whichpulse = whichpulse; % added 2009-10-23 JJC
dp{1}.s2candidate_areadiffs=s2candidate_areadiffs;
%dp{1}.a0 = a0;
dp{1}.evt_sum_per_ch = evt_sum_per_ch;
%dp{1}.p_p = p_p;
dp{1}.evt_mean = evt_mean;
dp{1}.evt_std = evt_std;
dp{1}.npts_above_thresh = npts_above_thresh;
dp{1}.head = head;
dp{1}.tail = tail;
%dp{1}.a1 = a1;
%dp{1}.a2 = a2;
dp{1}.t_mean = t_mean+pretrigger;
dp{1}.t_std = t_std+pretrigger;
dp{1}.peak_area_per_ch = peak_area_per_ch;
dp{1}.peak_height_per_ch = peak_height_per_ch; % added back 090804 pfs
dp{1}.peak_height_mV_per_ch = peak_height_mV_per_ch;
dp{1}.prompt_fraction = prompt_fraction; % added 090804 pfs - renamed from fp 20090915 jjc
dp{1}.preS1 = preS1; % added 090804 pfs

%         dp{1}.qqq = qqq;
%       dp{1}.xxx = 'hah!  you thought this was a real RQ!!';
%       dp{1}.xx = 'sold out.';
%dp{1}.coin = coin;
%dp{1}.outlier_bl = outlier_bl;
%dp{1}.summed_bl = summed_bl;
dp{1}.baseline_mV = baseline_mV;
dp{1}.sat = sat;
%dp{1}.bs = bs;
%dp{1}.flattened = flattened;
%dp{1}.digi = digi; % commented out 090804 pfs

%dp{1}.tLOAD = tLOAD;
%dp{1}.tANA = tANA;
%dp{1}.tTime = datestr(now); % record time stamp that job finished
dp{1}.evt_list = evt_list;
%dp{1}.verzsh = verzsh;
%dp{1}.par = par;
dp{1}.clock_time_sec = event_areas.clock_time_sec;
dp{1}.timestamp = timestamp;
