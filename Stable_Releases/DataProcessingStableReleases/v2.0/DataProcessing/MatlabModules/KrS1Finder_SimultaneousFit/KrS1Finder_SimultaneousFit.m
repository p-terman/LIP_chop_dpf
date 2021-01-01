function status = KrS1Finder_SimultaneousFit(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)
%
% This is a PulseQuantities module that fits to the two S1s in Kr83 events
% It does this by making a single summed trace from the start of the
% largest S1 to 1us later, then fitting a 2-S1 double double-expo to this.
%
% status = PulseQuantities_MinimumSet(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)
%
% Required RQs:
%   pulse_classification
%   xyz_corrected_pulse_area_bot_phe
%   pulse_start_samples
%
%
% Versioning:
%
%   v0.9 20140826 CHF&AC derived from Pulse_Quantities_Minimum_Set
%   v1.0 20140827 CHF&AC added minimum S1 size and fixed various bugs
%   v1.1 20140828 AC added 'fraction of area in max PMT for S1 a and b'RQs
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

myname = 'KrS1Finder_SimultaneousFit';
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
[cvt_struct settings2] = LUXCVTLoader_framework(data_path_evt,strrep(filename_evt,'evt','cvt'));

%% Initialize variables

N = length(cvt_struct);
min_s2bc_phe      = mymodule_settings.s2_selection.min_xyz_corrected_pulse_area_bot_phe;
max_s2bc_phe      = mymodule_settings.s2_selection.max_xyz_corrected_pulse_area_bot_phe;
max_roi_samples   = mymodule_settings.s1_selection.max_roi_samples;
excluded_channels = mymodule_settings.excluded_channels;
risetime_samples  = mymodule_settings.S1_pulse_shape.risetime_samples;
falltime_samples  = mymodule_settings.S1_pulse_shape.falltime_samples;
minFitAmp         = mymodule_settings.minFitAmp_phesamp;

% Fit options
fit_options = optimset('lsqcurvefit');
fit_options.Algorithm = 'trust-region-reflective';
fit_options.Display   = 'off';
fit_options.TolFun      = mymodule_settings.fit_opts.TolFun;
fit_options.TolX        = mymodule_settings.fit_opts.TolX;
fit_options.MaxFunEvals = mymodule_settings.fit_opts.MaxFunEvals;
fit_options.MaxIter     = mymodule_settings.fit_opts.MaxIter;


%Per-event RQs, to be filled when we find an event with Kr-like S2
dp.Kr83fit_s1a_area_phe          = zeros(1,N,'single');
dp.Kr83fit_s1b_area_phe          = zeros(1,N,'single');
dp.Kr83fit_dt_samples            = zeros(1,N,'single');
dp.Kr83fit_chisq                 = zeros(1,N,'single');
dp.Kr83fit_dof                   = zeros(1,N,'single');
dp.Kr83fit_exitflag              = zeros(1,N,'single');
dp.Kr83fit_s1a_max_pmt_fraction  = zeros(1,N,'single');
dp.Kr83fit_s1b_max_pmt_fraction  = zeros(1,N,'single');


%% Compute RQs

fprintf('\n');


for evt = 1:N
    
    if mod(evt,100) == 1
        fprintf('\n')
    end
    
    
    % ID active region Kr83s with 1 S2 in area range and >0 S1s before
    is_candidate_s2 = (dp.pulse_classification(:,evt)==2) ... 
        & inrange(dp.xyz_corrected_pulse_area_bot_phe(:,evt),[min_s2bc_phe,max_s2bc_phe]);
    if sum(is_candidate_s2)~=1
        continue
    end
    is_candidate_s1 = dp.pulse_classification(1:find(is_candidate_s2),evt)==1 ...
        & dp.pulse_area_phe(1:find(is_candidate_s2),evt)>mymodule_settings.s1_selection.min_area_phe;
    if sum(is_candidate_s1)==0
        continue
    end
    fprintf('.');
    
    % make a channel map that excludes APing PMTs
    ch_map = 1:122;
    ch_map = ch_map(~ismember(ch_map,excluded_channels));
    cvt_struct = LUXSumPOD_framework_new(cvt_struct,ch_map);
    
    % define ROI as [ start of first S1, min(max_roi_samples later, start of S2) )
    roi_start_samples=min(dp.pulse_start_samples(is_candidate_s1,evt));
    roi_end_samples=min(roi_start_samples+max_roi_samples,dp.pulse_start_samples(is_candidate_s2,evt));
    roi_length_samples = roi_end_samples-roi_start_samples;
    roi_cut = cvt_struct(evt).sumpod_time_samples >= roi_start_samples ...
                      & (cvt_struct(evt).sumpod_time_samples < roi_end_samples);
    roi_data_phe_per_sample = cvt_struct(evt).sumpod_data_phe_per_sample(roi_cut);
    roi_time_samples = cvt_struct(evt).sumpod_time_samples(roi_cut)';
    
    
    %-------- call the smoothing, preliminary pulsefinding, and fitting function ---------%
    
    fit_result = Kr83_Double_S1_Fitter(roi_time_samples, roi_data_phe_per_sample, minFitAmp, risetime_samples, falltime_samples, evt, fit_options);
    
    % fill RQs with fit parameters
    dp.Kr83fit_s1a_area_phe(evt)      = fit_result.Kr83fit_s1a_area_phe;
    dp.Kr83fit_s1b_area_phe(evt)      = fit_result.Kr83fit_s1b_area_phe;
    dp.Kr83fit_dt_samples(evt)        = fit_result.Kr83fit_dt_samples;
    dp.Kr83fit_chisq(evt)             = fit_result.Kr83fit_chisq;
    dp.Kr83fit_dof(evt)               = fit_result.Kr83fit_dof;
    dp.Kr83fit_exitflag(evt)          = fit_result.Kr83fit_exitflag;

    s1a_t_samples = fit_result.s1a_t_samples;
    s1b_t_samples = fit_result.s1b_t_samples;
    
    %------- calculate the `fraction in single PMT' rq -------------------%
    
    % calculate the maximum fractional area in any one PMT for S1a,b
    % define individual ROI_a and ROI_b +/- 2 rise/fall times
    window=max(2*risetime_samples,2);
    roi_a_start_samples = s1a_t_samples-window;
    window=max(2*falltime_samples,3);
    roi_a_end_samples = s1a_t_samples+window;
    
    window=max(2*risetime_samples,2);
    roi_b_start_samples = s1b_t_samples-window;
    window=max(2*falltime_samples,3);
    roi_b_end_samples = s1b_t_samples+window;
    
    % fill arrays with each channel's contribution to the two ROIs
    s1a_peak_area_phe=zeros(1,122,'single');
    s1b_peak_area_phe=zeros(1,122,'single');
    if ~isempty(ch_map)
      for ch = ch_map
        peak_cut = inrange(cvt_struct(evt).ch(ch).pod_time_samples,roi_a_start_samples,roi_a_end_samples+1);
        s1a_peak_area_phe(ch)=sum(cvt_struct(evt).ch(ch).pod_data_phe_per_sample(peak_cut));
        peak_cut = inrange(cvt_struct(evt).ch(ch).pod_time_samples,roi_b_start_samples,roi_b_end_samples+1);
        s1b_peak_area_phe(ch)=sum(cvt_struct(evt).ch(ch).pod_data_phe_per_sample(peak_cut));
      end
    end
    dp.Kr83fit_s1a_max_pmt_fraction(evt)=max(s1a_peak_area_phe)/sum(s1a_peak_area_phe);
    dp.Kr83fit_s1b_max_pmt_fraction(evt)=max(s1b_peak_area_phe)/sum(s1b_peak_area_phe);
    
    
end % for event evt


 
fprintf('Done! \n');


%% Write Output File

status = LUXBinaryRQWriter_framework(settings, dp, filename_rq, data_path_rq, event_number, livetime); % Should add livetime input at the end, remove event_number as settings field
