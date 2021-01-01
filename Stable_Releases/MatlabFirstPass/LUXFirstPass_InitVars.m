% LUXFirstPass_InitVars
%
% DO NOT CALL THIS MANUALLY, IT DOESN'T MAKE SENSE
%
% To be called from LUXFirstPass
% Create and initialize variables for use in LUXFirstPass
%
% 2010-01-27 JJC, DCM
%


%% Useful shortcuts
pts = s.daq_settings.sis3301.global.posttrigger + 1 + s.daq_settings.sis3301.global.pretrigger;
chs = length(s.daq_settings.sis3301.global.data_chs);

max_nb_s2s = s.ana_settings.max_nb_s2s;% how many potential S2-like pulses we want to find
max_nb_pulses = s.ana_settings.max_nb_pulses;% how many peaks total we will look for
pretrigger = s.daq_settings.sis3301.global.pretrigger;


%% Some other settings

evt_num_limit = 1;  % Number of chunks the ana_settingslysis will be divided into


%% Initialize variables to store things in

t0 =  -(pretrigger+1).*ones(max_nb_pulses,nb_evts_in_file); % t0 = 1
t10l = -(pretrigger+1).*ones(max_nb_pulses,nb_evts_in_file);
t50l = -(pretrigger+1).*ones(max_nb_pulses,nb_evts_in_file);
t10r = -(pretrigger+1).*ones(max_nb_pulses,nb_evts_in_file);
t50r = -(pretrigger+1).*ones(max_nb_pulses,nb_evts_in_file);
t1 = -(pretrigger+1).*ones(max_nb_pulses,nb_evts_in_file);
t2 = -(pretrigger+1).*ones(max_nb_pulses,nb_evts_in_file); % t2 = pts (the end of the trace)
ft10l = -(pretrigger+1).*ones(max_nb_s2s,nb_evts_in_file);
ft50l = -(pretrigger+1).*ones(max_nb_s2s,nb_evts_in_file);
ft10r = -(pretrigger+1).*ones(max_nb_s2s,nb_evts_in_file);
ft50r = -(pretrigger+1).*ones(max_nb_s2s,nb_evts_in_file);
ft1 = -(pretrigger+1).*ones(max_nb_s2s,nb_evts_in_file);
% t2l = pts * ones(itt_max,nb_evts_in_file);
% t2h = pts * ones(itt_max,nb_evts_in_file);
pikheight = -1 * ones(max_nb_pulses,nb_evts_in_file); % pikheight = -1
peak_height = -1 * ones(max_nb_pulses,nb_evts_in_file);
%       pikheight = -1 * ones(chs,itt_max,nb_evts_in_file); %new field 20090724, for looking at saturation from peak height - CHF
whichpik = zeros(max_nb_pulses,nb_evts_in_file);
whichpeak = zeros(max_nb_pulses,nb_evts_in_file);
whichpulse = zeros(max_nb_pulses,nb_evts_in_file);
s2candidate_areadiffs = zeros(max_nb_pulses,nb_evts_in_file);
%a0 = -1 * ones(1,nb_evts_in_file); % sum of whole trace
%a00 = -1 * ones(chs,nb_evts_in_file); % sum of each channel's trace
evt_sum_per_ch = -1 * ones(chs,nb_evts_in_file); % sum of each channel's trace
%p_p = -1 * ones(2,nb_evts_in_file); %
evt_mean = -1 * ones(1,nb_evts_in_file);
evt_std = -1 * ones(1,nb_evts_in_file);
%a1 = -1 * ones(max_nb_pulses,nb_evts_in_file); % a1 = -1
%a2 = -1 * ones(max_nb_pulses,nb_evts_in_file);
t_mean = -(pretrigger+1).*ones(max_nb_pulses,nb_evts_in_file);
t_std = -(pretrigger+1).*ones(max_nb_pulses,nb_evts_in_file);
pikarea_per_ch = zeros(max_nb_pulses , chs , nb_evts_in_file); % sum peak regions for each PMT
peak_area_per_ch = zeros(max_nb_pulses , chs , nb_evts_in_file); % sum peak regions for each PMT
pikheight_per_ch = -1.*ones(max_nb_pulses, chs, nb_evts_in_file); % height at t1 for each PMT
peak_height_per_ch = -1.*ones(max_nb_pulses, chs, nb_evts_in_file); % height at t1 for each PMT
pikheight_mV_per_ch = -1.*ones(max_nb_pulses, chs, nb_evts_in_file);
peak_height_mV_per_ch = -1.*ones(max_nb_pulses, chs, nb_evts_in_file);
prompt_fraction = -1 * ones(max_nb_pulses,nb_evts_in_file); %
preS1 = -1 * ones(max_nb_pulses,nb_evts_in_file); %
postS1 = -1 * ones(max_nb_pulses,nb_evts_in_file);
% qqq = zeros(itt_max, chs, nb_evts_in_file); % height at t0 for each PMT
% xxx = zeros(itt_max, chs, nb_evts_in_file); % naughtiest bit of each peak region for each PMT
% xx = zeros(chs, evts);           % a nice beer that I wish I had instead of each PMT
%bs = zeros(chs,nb_evts_in_file);
baseline_mV = -9999.*ones(max_nb_pulses,chs,nb_evts_in_file);
%digi = -1 * ones(max_nb_pulses,nb_evts_in_file);
%flattened = false(chs,nb_evts_in_file);
%coin = false(max_nb_pulses,chs,nb_evts_in_file);
%outlier_bl = zeros(3,chs,nb_evts_in_file); % new field for data quality monitoring 070125 pfs
%summed_bl = zeros(4,nb_evts_in_file); % new field for data quality monitoring 070125 pfs
timestamp = -1.*ones(1,nb_evts_in_file);
npts_above_thresh = -1.*ones(1,nb_evts_in_file);
head = -9999.*ones(max_nb_pulses,chs,nb_evts_in_file); % modified to be n_preceeding and n_following from POD mode per channel
tail = -9999.*ones(max_nb_pulses,chs,nb_evts_in_file);
sat = zeros(chs,nb_evts_in_file);
event_areas.clock_time_sec = 0;

