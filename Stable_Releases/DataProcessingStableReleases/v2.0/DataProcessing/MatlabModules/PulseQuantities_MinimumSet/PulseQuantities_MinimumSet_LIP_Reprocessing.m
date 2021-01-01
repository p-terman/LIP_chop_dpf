function status = PulseQuantities_MinimumSet_LIP_Reprocessing(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)
%
% This is a PulseQuantities module that calculates the MinimumSet in Matlab
%
% status = PulseQuantities_MinimumSet(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)
%
%
%
%
%
%
% Required RQs:
%
%
%
%
%
% Versioning:
%
%   v1.0 20121115 CHF - Created
%        20121213 CHF - Fixed issue with indexing when same channel fired more
%                  than once for the same pulse
%        20121218 JRV - Changed signature to use paths.
%                  Now does not load any DP or IQ settings from admin!
%   v2.0 20130212 CHF - Big changes.
%                * Changed pulse_start_timestamp_samples (and end) to
%                  pulse_start_samples, which now uses trigger-based timing,
%                  instead of global timestamp.
%                * Added fields pulse_start_pod, pulse_end_pod
%                * Changed all references from cvt_struct(evt).chsum.pod()
%                  to cvt_struct(evt).sumpod(), and fixed new subfield
%                * Re-did logic for per-channel RQs; using new pulse_start_pod
%                  and pulse_start_pod simplified it.
%                * Treating per ch RQs with multiple pod contributions to
%                  sumpod as a larger pod by clustering them together.
%                * Updated name sumpod_data_phe to sumpod_data_phe_per_sample.
%        20130216 CHF - Minor fixes.
%        20130219 CHF - Added event_timestamp_samples rq.
%        20130304 JRV - Added if statements enforcing length around pre/post
%                baseline calculations
%        20130305 CHF - Changed pulse_detect_... -> peak_detect_..
%   v1.5 20130319 CHF - Fixed prompt fraction - there was a problem with time
%                definitions that resulted in most values being 0.
%                LUXPromptFraction_framework expects a sample number relative to
%                start of pulse for t10l, while it was being provided the
%                time relative to trigger.
%        20130322 CHF - Minor improvements to per ch computations. No actual RQ
%                changes.
%        20130325 pfs - Minor changes. modified comment on line 178 to read "pulse" 
%                instead of "pod." Also modified variable name in prompt fraction
%                code, for consistency; this is name-only, not algorithmic.
%              - Major change: added -1 to tleft variable for prompt
%                fraction calculation. it looks more correct, based on the data.
%       20130327 pfs - wrapped a try-catch around LUXExponentialFit_framework. this is pretty
%                ghetto, but I don't understand the error yet and don't have the time now.
%       20130328 pfs - fixed major issue with per-channel RQs. see line ~258 and below.
%       20130413 CHF - dp.peak_detect... were being initialized incorrectly
%                as dp.pulse_detect...
%       20130315 JRV - Now skips files with 0 events
%       20130315 pfs - requires new xml parameter thr_phe (should have done this ~weeks ago)
%                      default value 0.05 based on nothing other than past precedent.
%       20130503 JRV - Added AC's rq.rms_width_samples RQ
%                      Commented out PFS's thr_phe RQs until XML setings
%                      are updated
%       20130523 JJC - uncommented PFS's thr_phe RQs because XML settings
%       appear to be updated, and this RQ is darned nice to have. Also add
%       section in Bookkeeping that will give this thr_phe value a default
%       of 0.05 if it is not found in the xml settings file (should have
%       done this to begin with).
%       20130529 JRV - Fits now accept tolerance parameters from XML
%       20130626 pfs - added skinny_peak_height_mV
%       20140204 SS - Added calculations from TRC: pre/post pulse PMT sums
%                amis1_fraction, n_samples_in_evt, full_evt_area_phe
%                Also redefined pulse_length_samples as did not include
%                empty sumpods, added new RQ pulse_length_sumpods to take
%                original definition of pulse_length_samples
%       20140405 AC  - Modified the loop over channels to process only
%       channels that are non-empty in event_struct (loads 1:122,hard-coded) 
%       *and* cvt_struct (might not), since code in the loop needs both.
%       20140502 AC - Fixed rms_width_samples to be the rms (not ms) 
%       deviation of pulse area in time; added skewness, kurtosis in time
%       20140605 SS - Added calculation of skinny pulse start and end times
%                     for timing module compatibility
%       20140805 SS - Added protection against pulses with index_kept_sumpods
%                     = 0 for module compatibility
%       20140528 AL - The type of each RQ is now clearly defined, instead of
%                     using the default (double).
%                     - Changed initialisation values for some of the RQs, mainly
%                     to move away from NANs.
%                     - event_timestamp_samples and pulse_height_timing_samples
%                     are no longer created (they were just commented out in
%                     the code for now).
%                     - Avoided situation in which skinny_peak_height_mV could
%                     be a NAN - when there is no pod data for a channel inside
%                     the skinny_peak interval). It now gets 0.
%       20140616 AC - If an acquisition was taken with delay_buffer=0,
%       subtract the per-channel spurious area, as reported by IQ
%       per_pod_spurious_area, from peak_area_phe and pulse_area_phe 

%20180423 PAT made this version that will only process things that need to
%be reprocessed for LIP analysis. 
% RQ versions:
%
%
%
%% Load .rq file

status = [];

dp = LUXLoadRQ1s_framework(filename_rq, data_path_rq);

settings.evt_settings = dp.admin.evt_settings;
settings.daq_settings = dp.admin.daq_settings;
settings.filename_prefix = dp.admin.filename_prefix;

event_number = dp.event_number;
livetime = dp.admin.livetime;

delay_buffer = dp.admin.daq_settings.sis3301.global.delay_buffer;

%% Bookkeeping

myname = 'PulseQuantities_MinimumSet';
fprintf('\n\n *** Starting module %s\n',myname);

if isempty(event_number)
    fprintf('\n\n *** Skipping module (no events in file) %s\n',myname);
    return
end

dp_settings_xml = XMLReader_framework(data_processing_xml_path);
lug_iqs_xml = XMLReader_framework(iq_xml_path);

module_names = {dp_settings_xml.data_processing_settings.module.module_name};
index_temp = strfind(module_names,myname);
index_module = find(not(cellfun('isempty', index_temp)));

mymodule_settings = dp_settings_xml.data_processing_settings.module(index_module).parameters;

max_num_pulses = dp_settings_xml.data_processing_settings.global.max_num_pulses;

%% Load .cvt file or calculate it

filename_cvt = strrep(filename_evt,'evt','cvt');

if ~exist([data_path_evt filesep filename_cvt],'file')
    fprintf('Did not find .cvt file. Running Summer Module\n');
    status = PulseCalibration_BaselineZen(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path);
    status = PODSummer_LUXSumPOD(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path);
end

fprintf('Loading sum from %s\n',filename_cvt);
[cvt_struct settings] = LUXCVTLoader_framework(data_path_evt,strrep(filename_evt,'evt','cvt'));

% We still need the mV data...
event_struct = LUXEventLoader_framework(data_path_evt, filename_evt);


%% Get IQs to calculate spurious_area in phe

% Figure out which iq has the pmt_gains - assuming only ONE iq was returned
% for each type
pmt_gains_mVns_per_phe = [];

for ii = 1:length(lug_iqs_xml.iq)
    if isfield(lug_iqs_xml,'iq') && isfield(lug_iqs_xml.iq(ii),'global') && isfield(lug_iqs_xml.iq(ii).global,'iq_type')
        if strcmp(lug_iqs_xml.iq(ii).global.iq_type,'pmt_gains') == 1
            pmt_gains_mVns_per_phe = [lug_iqs_xml.iq(ii).fit.channel.mVns_per_phe];
            break
        end
    end
end

if isempty(pmt_gains_mVns_per_phe);
    error('*** ERROR: No PMT Gain Calibrations were provided!');
end

% Figure out which iq has the per_pod_spurious_area - assuming 1 iq was returned
% for each type
per_pod_spurious_area_phe = [];

for ii = 1:length(lug_iqs_xml.iq)
    if isfield(lug_iqs_xml,'iq') && isfield(lug_iqs_xml.iq(ii),'global') && isfield(lug_iqs_xml.iq(ii).global,'iq_type')
        if strcmp(lug_iqs_xml.iq(ii).global.iq_type,'per_pod_spurious_area') == 1
            per_pod_spurious_area_phe = [lug_iqs_xml.iq(ii).fit.channel.per_pod_spurious_area_mVns]./pmt_gains_mVns_per_phe;
            break
        end
    end
end

if isempty(per_pod_spurious_area_phe);
    error('*** ERROR: No spurious-area IQ provided!');
end

%% Initialize variables

pmt_chs = 1:122;

N = length(dp.event_number); %changed PAT 

pulse_event_size = [max_num_pulses N];
per_channel_event_size = [max_num_pulses length(pmt_chs) N];

% should make these sparse!


% Per sum
dp.pulse_area_phe                       = zeros(pulse_event_size,'single');
dp.skinny_pulse_area_phe                = zeros(pulse_event_size,'single');
dp.pulse_area_positive_phe              = zeros(pulse_event_size,'single');
dp.pulse_area_negative_phe              = zeros(pulse_event_size,'single');
dp.pulse_height_phe_per_sample          = zeros(pulse_event_size,'single');
dp.prompt_fraction                      = zeros(pulse_event_size,'single');
dp.prompt_fraction_tlx                  = zeros(pulse_event_size,'single');
dp.top_bottom_ratio                     = zeros(pulse_event_size,'single');
dp.top_bottom_asymmetry                 = zeros(pulse_event_size,'single');

dp.exp_fit_amplitude_phe_per_sample     = zeros(pulse_event_size,'single');
dp.exp_fit_tau_fall_samples             = zeros(pulse_event_size,'single');
dp.exp_fit_time_offset_samples          = zeros(pulse_event_size,'single');
dp.exp_fit_tau_rise_samples             = zeros(pulse_event_size,'single');
dp.exp_fit_chisq                        = zeros(pulse_event_size,'single');
dp.exp_fit_dof                          = zeros(pulse_event_size,'int32');

dp.gaus_fit_amplitude_phe_per_sample    = zeros(pulse_event_size,'single');
dp.gaus_fit_mu_samples                  = zeros(pulse_event_size,'single');
dp.gaus_fit_sigma_samples               = zeros(pulse_event_size,'single');
dp.gaus_fit_chisq                       = zeros(pulse_event_size,'single');
dp.gaus_fit_dof                         = zeros(pulse_event_size,'int32');

dp.s2filter_max_area_diff               = zeros(pulse_event_size,'single');
dp.s2filter_max_s2_area                 = zeros(pulse_event_size,'single');
dp.s2filter_max_s1_area                 = zeros(pulse_event_size,'single');

dp.pulse_mean_samples                   = zeros(pulse_event_size,'single');
dp.rms_width_samples                    = zeros(pulse_event_size,'single');
dp.pulse_skewness                       = zeros(pulse_event_size,'single');
dp.pulse_kurtosis                       = zeros(pulse_event_size,'single');

dp.pulse_length_samples                 = zeros(pulse_event_size,'uint32');

dp.skinny_pulse_start_samples           = zeros(pulse_event_size,'int32');
dp.skinny_pulse_end_samples             = zeros(pulse_event_size,'int32');
% initialise values to -999999
dp.skinny_pulse_start_samples(:)       = -999999;
dp.skinny_pulse_end_samples(:)         = -999999;

preBins     = mymodule_settings.prompt_fraction.preBins_samples;
windowBins  = mymodule_settings.prompt_fraction.windowBins_samples;

s1window_samples = mymodule_settings.s2filter.s1window_samples;
s2window_samples = mymodule_settings.s2filter.s2window_samples;

% Per channel
dp.peak_area_phe                        = zeros(per_channel_event_size,'single');
dp.peak_area_negative_phe               = zeros(per_channel_event_size,'single');
dp.peak_height_phe_per_sample           = zeros(per_channel_event_size,'single');
dp.mean_first_last_pts_phe_per_sample   = zeros(per_channel_event_size,'single');
dp.peak_height_mV                       = zeros(per_channel_event_size,'single');
dp.daq_saturation_flag                  = zeros(per_channel_event_size,'uint8');
dp.pmt_2pct_saturation_flag             = zeros(per_channel_event_size,'uint8');
dp.peak_detect_pretrigger_mean_phe      = zeros(per_channel_event_size,'single');
dp.peak_detect_posttrigger_mean_phe     = zeros(per_channel_event_size,'single');
dp.peak_detect_pretrigger_std_phe       = zeros(per_channel_event_size,'single');
dp.peak_detect_posttrigger_std_phe      = zeros(per_channel_event_size,'single');
dp.baseline_daq_mV                      = zeros(per_channel_event_size,'single');

% new RQs --pfs
dp.peak_area_thr_phe                    = zeros(per_channel_event_size,'single');
dp.skinny_peak_height_mV                = zeros(per_channel_event_size,'single');
dp.skinny_peak_area_phe                 = zeros(per_channel_event_size,'single');
  
% RQs originally in TRC  - SS
dp.pre_pulse_area_positive_phe          = zeros(pulse_event_size,'single');
dp.pre_pulse_area_negative_phe          = zeros(pulse_event_size,'single');
dp.post_pulse_area_positive_phe         = zeros(pulse_event_size,'single');
dp.post_pulse_area_negative_phe         = zeros(pulse_event_size,'single');
dp.amis1_fraction                       = zeros(pulse_event_size,'single');
dp.n_samples_in_evt                     = zeros(1,N,'uint32');
dp.full_evt_area_phe                    = zeros(1,N,'single');

dp.pulse_length_sumpods                 = zeros(pulse_event_size,'uint32');


%% Compute RQs

fprintf('\n');

% Removing this RQ as it is already calculated in the InitialiseRQFile module.
% Will just leave it commented out for now, to be deleted in future revisions.
% AL - 140528
%dp.event_timestamp_samples = [event_struct(:).timestamp];

for i = 1:length(dp.event_number)% changed from 1:N - PAT, now only does ones that need reprocessing
    evt =  dp.event_number(i);
    if mod(evt,100) == 1
        fprintf('\n')
    end
    
    fprintf('.');
    
    % Check if fields exist and has data
    if isfield(cvt_struct(evt),'sumpod_data_phe_per_sample') && cvt_struct(evt).empty == 0
        loop_max = min([sum(isfinite(dp.index_kept_sumpods(:,i))) max_num_pulses]); %PAT changed i from evt

        % for every pulse (NOT pod)
        for pp = 1:loop_max
            %dis('pp=%d,evt=%d ... index_kept=%d\n',pp,evt,dp.index_kept_sumpods(pp,evt))
            % -1 and +1 signs are needed to get the right time cut
            %pulse_cut = inrange(cvt_struct(evt).sumpod_time_samples,dp.pulse_start_samples(pp,evt)-1,dp.pulse_end_samples(pp,evt)+1);

            if dp.index_kept_sumpods(pp,i) == 0 %PAT changed to i from evt
                continue
            end

            % use explicit range cuts, not "inrange" function
            pulse_cut = (cvt_struct(evt).sumpod_time_samples >= dp.pulse_start_samples(pp,i)) ...
                      & (cvt_struct(evt).sumpod_time_samples <= dp.pulse_end_samples(pp,i)); %PAT changed evt to i 
                  
            if sum(pulse_cut)==0 % if-then statment added 160301 by PAT this excludes pulses that have no pulse data in a time period. 
                dp.pulse_start_samples(pp,i)=999999;
                dp.pulse_end_samples(pp,i)=999999; %PAT changed evt to i
                continue
                
            end
            
            pulse_data_phe = cvt_struct(evt).sumpod_data_phe_per_sample(pulse_cut);
            pulse_data_phe_full = cvt_struct(evt).sumpod_data_phe_per_sample(cvt_struct(evt).sumpod_time_samples<1e5);
            pulse_time_samples = cvt_struct(evt).sumpod_time_samples(pulse_cut);
            pulse_time_samples_full = cvt_struct(evt).sumpod_time_samples(cvt_struct(evt).sumpod_time_samples<1e5);
            pulse_length = dp.pulse_end_samples(pp,i) - dp.pulse_start_samples(pp,i)+1; %changed to i fro evt PAT
            
            dp.n_samples_in_evt (i) = length(pulse_time_samples_full); %PAT changed
            dp.full_evt_area_phe (i) = sum(pulse_data_phe_full);
          
           
                % ***Calculate PULSE quantities per pulse here***
            
            % Pulse area -- see below, it is calculated from peak_area_phe 
%            dp.pulse_area_phe(pp,evt) = sum(pulse_data_phe);
%            dp.pulse_area_thr_phe(pp,evt) = sum(pulse_data_phe(pulse_data_phe>mymodule_settings.thr_phe));
%            dp.pulse_area_positive_phe(pp,evt) = sum(pulse_data_phe(pulse_data_phe >= 0));
%            dp.pulse_area_negative_phe(pp,evt) = sum(pulse_data_phe(pulse_data_phe < 0));
            
                       
          %Skinny Boxcar
          fullBoxArea = 0;
          for ii = 1:pulse_length
            if ii <= numel(pulse_data_phe)
              fullBoxArea = fullBoxArea + pulse_data_phe(ii);
            else 
              break
            end
          end
          % Initiliaze values
          maxSkinnyBoxArea = 0;
          tt_maxSkinnyBoxArea = 0;
          
          % Find skinny box with maximum area
          for tt = 1:(pulse_length-mymodule_settings.skinnyBoxSamples+1);
            skinnyBoxArea = 0;
            for ii=tt:tt+mymodule_settings.skinnyBoxSamples-1
              if ii <= numel(pulse_data_phe)
                skinnyBoxArea = skinnyBoxArea + pulse_data_phe(ii);
              else 
                break
              end
            end
              if skinnyBoxArea > maxSkinnyBoxArea
                maxSkinnyBoxArea = skinnyBoxArea;
                tt_maxSkinnyBoxArea = tt;
              end 
          end
            
          % When pulse is shorter than the skinny box:
          if mymodule_settings.skinnyBoxSamples > pulse_length
             tt_maxSkinnyBoxArea = 1;
             maxSkinnyBoxArea = fullBoxArea;     
          end
           
          % Record values
          dp.skinny_pulse_start_samples(pp,i) = dp.pulse_start_samples(pp,i) + tt_maxSkinnyBoxArea-1;
          dp.skinny_pulse_end_samples(pp,i)   = dp.pulse_start_samples(pp,i) + tt_maxSkinnyBoxArea + mymodule_settings.skinnyBoxSamples-1;
    
          % amis1 - measure of S1-ness 0=>S2, 1=>S1 (ish)
          dp.amis1_fraction(pp,i) = maxSkinnyBoxArea/fullBoxArea;
          if isnan(dp.amis1_fraction(pp,i))
            dp.amis1_fraction(pp,i) = 0; % reasonable asymptote
          end
          if isinf(dp.amis1_fraction(pp,i))
            dp.amis1_fraction(pp,i) = -1; % reasonable asymptote 
          end
          if dp.amis1_fraction(pp,i) > 10 || dp.amis1_fraction(pp,i) < -10
            dp.amis1_fraction(pp,i) = -2; % error code
          end
          
            
            % Pulse height
            [dp.pulse_height_phe_per_sample(pp,i) loc] = max(pulse_data_phe);

            % Removing this RQ. Will just leave it commented out for now,
            % to be deleted in future revisions. AL - 140528
            %dp.pulse_height_timing_samples(pp,evt) = loc;
            
            % Prompt fraction
            tleft = find(pulse_time_samples == (dp.aft_t0_samples(pp,i))); % I added the -1, based on the event profiles. we may need to change this as algorithms evolve ..pfs
            tlx = find(pulse_time_samples == (dp.aft_tlx_samples(pp,i)));
            dp.prompt_fraction(pp,i) = LUXPromptFraction_framework(pulse_data_phe,preBins,windowBins,tleft);
            dp.prompt_fraction_tlx(pp,i) = LUXPromptFraction_framework(pulse_data_phe,preBins,windowBins,tlx);

            try
            % Exponential fits
            if isfield(mymodule_settings,'exponential_fit')
                efit_settings = mymodule_settings.exponential_fit;
            else
                efit_settings = [];
            end
            [efit_params, dp.exp_fit_chisq(pp,i) dp.exp_fit_dof(pp,i)] = LUXExponentialFit_framework(1:length(pulse_data_phe),pulse_data_phe,efit_settings);
            catch;efit_params=[0 0 0 0];end; % kludge since I don't understand the error (yet)

            dp.exp_fit_amplitude_phe_per_sample(pp,i) = efit_params(1);
            dp.exp_fit_tau_fall_samples(pp,i) = efit_params(2);
            dp.exp_fit_time_offset_samples(pp,i) = efit_params(3); % this is relative to the start of the PULSE!
            dp.exp_fit_tau_rise_samples(pp,i) = efit_params(4);
            clear efit_params;
            
            % Gaussian fits
            if isfield(mymodule_settings,'gaussian_fit')
                gfit_settings = mymodule_settings.gaussian_fit;
            else
                gfit_settings = [];
            end
            [gfit_params, dp.gaus_fit_chisq(pp,i) dp.gaus_fit_dof(pp,i)] = LUXGaussianFit_framework(1:length(pulse_data_phe),pulse_data_phe,gfit_settings);
            
            dp.gaus_fit_amplitude_phe_per_sample(pp,i) = gfit_params(1);
            dp.gaus_fit_mu_samples(pp,i) = gfit_params(2);
            dp.gaus_fit_sigma_samples(pp,i) = gfit_params(3);
            clear gfit_params;
            
            % S2Filter
            ext_pulse = zeros(1,length(pulse_data_phe) + 2*s2window_samples); % extend pulse so that filters are not longer than data
            ext_pulse((s2window_samples+1) : (s2window_samples+length(pulse_data_phe))) = pulse_data_phe;
            
            [max_area_diff max_s2_area max_s1_area] = LUXS2FilterMatlab_framework(ext_pulse,s2window_samples,s1window_samples);
            
            dp.s2filter_max_area_diff(pp,i) = max_area_diff;
            dp.s2filter_max_s2_area(pp,i) = max_s2_area;
            dp.s2filter_max_s1_area(pp,i) = max_s1_area;
            
            % sample central moments of the distribution in time
            % calculate only for samples from aft_t0 to
            % aft_t2
            rms_width_range = inrange(pulse_time_samples,[dp.aft_t0_samples(pp,i), dp.aft_t2_samples(pp,i)]');
            dp.pulse_mean_samples(pp,i) = sum(pulse_data_phe(rms_width_range).*pulse_time_samples(rms_width_range))./ sum(pulse_data_phe(rms_width_range));
            dp.rms_width_samples(pp,i) = sum(pulse_data_phe(rms_width_range).*(pulse_time_samples(rms_width_range) - dp.pulse_mean_samples(pp,i)).^2) ./ sum(pulse_data_phe(rms_width_range));
            if(dp.rms_width_samples(pp,i)>0)
                dp.rms_width_samples(pp,i)=sqrt(dp.rms_width_samples(pp,i));
            else
                dp.rms_width_samples(pp,i)=-1.0*sqrt(-1.0*dp.rms_width_samples(pp,i));
            end
            dp.pulse_skewness(pp,i) = sum(pulse_data_phe(rms_width_range).*(pulse_time_samples(rms_width_range) - dp.pulse_mean_samples(pp,i)).^3) ./ sum(pulse_data_phe(rms_width_range));
            dp.pulse_skewness(pp,i)=dp.pulse_skewness(pp,i)/(dp.rms_width_samples(pp,i)^3);
            dp.pulse_kurtosis(pp,i) = sum(pulse_data_phe(rms_width_range).*(pulse_time_samples(rms_width_range) - dp.pulse_mean_samples(pp,i)).^4) ./ sum(pulse_data_phe(rms_width_range));
            dp.pulse_kurtosis(pp,i)=dp.pulse_kurtosis(pp,i)/(dp.rms_width_samples(pp,i)^4)-3.0;
            

            
            %Summing of pre and post pulse PMT signals
         
            tmp_pre_t0 = max(pulse_time_samples_full(1), dp.pulse_start_samples(pp,i)-mymodule_settings.preBoxSamples);  
            tmp_pre_range = find((pulse_time_samples_full >= tmp_pre_t0) & (pulse_time_samples_full <= (dp.pulse_start_samples(pp,i)-1)));
          
            tmp_post_t2 = min(pulse_time_samples_full(end),dp.pulse_end_samples(pp,i) + mymodule_settings.preBoxSamples);
            tmp_post_range = find((pulse_time_samples_full >= (dp.pulse_end_samples(pp,i)+1)) & (pulse_time_samples_full<=tmp_post_t2));    
  
            tmp_pre_pulse_data_phe = pulse_data_phe_full(tmp_pre_range);
            tmp_post_pulse_data_phe = pulse_data_phe_full(tmp_post_range);
 
            dp.pre_pulse_area_positive_phe(pp,i) = sum(tmp_pre_pulse_data_phe(tmp_pre_pulse_data_phe>=0));
            dp.pre_pulse_area_negative_phe(pp,i) = sum(tmp_pre_pulse_data_phe(tmp_pre_pulse_data_phe<0));
            dp.post_pulse_area_positive_phe(pp,i) = sum(tmp_post_pulse_data_phe(tmp_post_pulse_data_phe>=0));
            dp.post_pulse_area_negative_phe(pp,i) = sum(tmp_post_pulse_data_phe(tmp_post_pulse_data_phe<0));
            
            % Pulse Length
            dp.pulse_length_samples(pp,i) = pulse_length;
            
            % If the following rq is shorter than pulse_length_samples there are missing sumpods in pulse_data_phe
            dp.pulse_length_sumpods(pp,i) = length(pulse_data_phe);  
            
            total_sum_peak_area_phe = 0;
            
            
            % PER CH QUANTITIES
            ch_map = find(~[cvt_struct(evt).ch(:).empty]);
            ch_map = intersect(ch_map,find(~[event_struct(evt).ch(:).empty])); %need the phe and the mV scales
            
            if ~isempty(ch_map)
                
                for ch = ch_map
                    
                        % Make slice of data
                        peak_cut = inrange(cvt_struct(evt).ch(ch).pod_time_samples,dp.pulse_start_samples(pp,i)-1,dp.pulse_end_samples(pp,i)+1);

                    if ~isempty(peak_cut) && sum(peak_cut) > 0
                        
                        peak_data_phe = cvt_struct(evt).ch(ch).pod_data_phe_per_sample(peak_cut);
                        peak_time_samples = cvt_struct(evt).ch(ch).pod_time_samples(peak_cut);
                        peak_data_mV = event_struct(evt).ch(ch).pod_data_mV(peak_cut);
                        
                        % This is to find which baseline_mV to use. Since
                        % we are cutting pods in time, it's not
                        % straightforward to know which pod baseline 'twas                            
                        % But this statement takes care of that
                        peak_inds = unique(cumsum(ismember(cvt_struct(evt).ch(ch).pod_time_samples,cvt_struct(evt).ch(ch).pod_start_samples)));
                        
                        dp.baseline_daq_mV(pp,ch,i) = mean(event_struct(evt).ch(ch).pod_baseline_mV(peak_inds));
                        
                        % Peak area, subtracting spurious_area where it
                        % exists if delay_buffer==0
                        dp.peak_area_phe(pp,ch,i) = sum(peak_data_phe);
                        if(ch<=length(per_pod_spurious_area_phe)&&delay_buffer==0);
                            dp.peak_area_phe(pp,ch,i) = dp.peak_area_phe(pp,ch,i)-per_pod_spurious_area_phe(ch);
                        end;
                        dp.peak_area_negative_phe(pp,ch,i) = sum(peak_data_phe(peak_data_phe<0));
                        dp.peak_area_thr_phe(pp,ch,i) = sum(peak_data_phe(peak_data_phe>mymodule_settings.thr_phe));

						% new RQ: skinnyBox peak area
						skinny_peak_cut = inrange(cvt_struct(evt).ch(ch).pod_time_samples,dp.skinny_pulse_start_samples(pp,i)-1,dp.skinny_pulse_end_samples(pp,i)+1);				
						dp.skinny_peak_area_phe(pp,ch,i) = sum(cvt_struct(evt).ch(ch).pod_data_phe_per_sample(skinny_peak_cut));
						tmp = max(event_struct(evt).ch(ch).pod_data_mV(skinny_peak_cut));
						if isempty(tmp); tmp = 0; end;
						dp.skinny_peak_height_mV(pp,ch,i) = tmp;

                        % For debugging
                        total_sum_peak_area_phe = total_sum_peak_area_phe + sum(peak_data_phe);
                        
                        % Peak height
                        dp.peak_height_phe_per_sample(pp,ch,i) = max(peak_data_phe);
                        dp.peak_height_mV(pp,ch,i) = max(peak_data_mV);
                        
                        dp.mean_first_last_pts_phe_per_sample(pp,ch,i) = mean(peak_data_phe([1 end]));
                        
			% Pre and post trigger baseline rqs. Corrected to use POD array instead of peak array - AL - 2013/08/02
                        if length(cvt_struct(evt).ch(ch).pod_data_phe_per_sample) >= 24 % added to fix occational crash - JRV - 2013-03-04
                            dp.peak_detect_pretrigger_mean_phe(pp,ch,i) = mean(cvt_struct(evt).ch(ch).pod_data_phe_per_sample(1:24));
                            dp.peak_detect_pretrigger_std_phe(pp,ch,i)  = std(cvt_struct(evt).ch(ch).pod_data_phe_per_sample(1:24));
                        end
                        
                        if length(cvt_struct(evt).ch(ch).pod_data_phe_per_sample) >= 31 % added to fix occational crash - JRV - 2013-03-04
                            dp.peak_detect_posttrigger_mean_phe(pp,ch,i) = mean(cvt_struct(evt).ch(ch).pod_data_phe_per_sample((end-31+1):end));
                            dp.peak_detect_posttrigger_std_phe(pp,ch,i) = std(cvt_struct(evt).ch(ch).pod_data_phe_per_sample((end-31+1):end));
                        end
                        
                        % Saturation flags
                        dp.daq_saturation_flag(pp,ch,i) = dp.peak_height_mV(pp,ch,i) > mymodule_settings.saturation_flags.daq_saturation_threshold_mV;
                        dp.pmt_2pct_saturation_flag(pp,ch,i) = dp.peak_height_mV(pp,ch,i) > mymodule_settings.saturation_flags.pmt_saturation_2pct_threshold_mV;
                    end % fi


                end % for ch

                dp.pulse_area_phe(pp,i) = squeeze(sum(dp.peak_area_phe(pp,:,i),2));
                dp.skinny_pulse_area_phe(pp,i) = squeeze(sum(dp.skinny_peak_area_phe(pp,:,i),2));
                
                % For debugging sum(peak_area_phe) ~= pulse_area_phe
                if 0
                    if ~inrange(total_sum_peak_area_phe./dp.pulse_area_phe(pp,i),0.99,1.001)
                        
                        figure(43423); clf;
                        h1 = subplot(2,1,1); hold on
                        plot(pulse_time_samples,pulse_data_phe,'k.-')
                        plot(cvt_struct(evt).sumpod_time_samples,cvt_struct(evt).sumpod_data_phe_per_sample,'ro-')
                        title(sum(pulse_data_phe));
                        
                        h2 = subplot(2,1,2); hold on
                        
                        for ch = 1:122
                            peak_cut = inrange(cvt_struct(evt).ch(ch).pod_time_samples,dp.pulse_start_samples(pp,i)-1,dp.pulse_end_samples(pp,i)+1);
                            if sum(peak_cut); fprintf('%d ',ch); end
                            peak_data_phe = cvt_struct(evt).ch(ch).pod_data_phe_per_sample(peak_cut);
                            peak_time_samples = cvt_struct(evt).ch(ch).pod_time_samples(peak_cut);
                            plot(peak_time_samples,peak_data_phe,'k.-')
                            plot(cvt_struct(evt).ch(ch).pod_time_samples,cvt_struct(evt).ch(ch).pod_data_phe_per_sample,'ro-')
                        end
                        title(total_sum_peak_area_phe);
                        linkaxes([h1 h2],'xy');
                        xlim([dp.pulse_start_samples(pp,i)-10 dp.pulse_end_samples(pp,i)+10])
                        
                        
                    end
                end
                
                
            end % fi            
            
            top = sum(squeeze(dp.peak_area_phe(pp,[1:60 121],i)),2);
            bot = sum(squeeze(dp.peak_area_phe(pp,[61:120 122],i)),2);
            
            dp.top_bottom_ratio(pp,i) =  (top)./(bot);
            dp.top_bottom_asymmetry(pp,i) =  (top-bot)./(top+bot);
            
            
           
            
            
 
       
            
            
        end % for pulse pp

    end % fi sumpod
    
end % for event i / evt


 
fprintf('Done! \n');


%% Write Output File

status = LUXBinaryRQWriter_framework(settings, dp, filename_rq, data_path_rq, event_number, livetime); % Should add livetime input at the end, remove event_number as settings field
