function status = PulseTiming_PerusePeeksMatlab(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)
%
% This is a PulseTiming module using the PerusePeeksMatlab algorithm.
%
% function status = PulseTiming_PerusePeeksMatlab(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)
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
%   v1.0 20121025 CHF - Created
%        20121218 JRV - Changed signature to use paths for xmls.
%                       Now does not load any DP or IQ settings from admin!
%   v2.0 20130212 CHF - Big changes.
%                * Changed pulse_start_timestamp_samples (and end) to 
%                  pulse_start_samples, which now uses trigger-based timing, 
%                  instead of global timestamp.
%                * Added fields pulse_start_pod, pulse_end_pod
%                * Changed all references from cvt_struct(evt).chsum.pod() 
%                  to cvt_struct(evt).sumpod(), and fixed new subfield names
%                * Made sure all calls are to _framework code versions
%                * Updated name sumpod_data_phe to sumpod_data_phe_per_sample.       
%        20130216 CHF - Minor fixes.
%   v2.5 20130401 CHF - Changed variable names to conform to standards
%                (aft_t0_samples, aft_t1_samples, etc).
%        20130315 JRV - Now skips files with 0 events
%
% RQ versions:
%
%
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

%% Bookkeeping

myname = 'PulseTiming_PerusePeeksMatlab';
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
    fprintf('*** WARNING: Did not find .cvt file. Running Summer Module\n');
    status = PulseCalibration_BaselineZen(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path);
    status = PODSummer_LUXSumPOD(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path);    
end

fprintf('Loading sum from %s\n',filename_cvt);
[cvt_struct settings] = LUXCVTLoader_framework(data_path_evt,strrep(filename_evt,'evt','cvt'));

%% Initialize variables

N = length(cvt_struct);

pulse_event_size = [max_num_pulses N];

dp.aft_box_open_samples = zeros(pulse_event_size);
dp.aft_t0_samples = zeros(pulse_event_size);
dp.aft_t1_samples = zeros(pulse_event_size);
dp.aft_t2_samples = zeros(pulse_event_size);
dp.aft_box_close_samples = zeros(pulse_event_size);
dp.preBoxArea = zeros(pulse_event_size);
dp.postBoxArea = zeros(pulse_event_size);

% Shortcuts

BoxWidthSamples = mymodule_settings.fullBoxSamples;
edgeFraction = mymodule_settings.edgeFraction;

%% Compute timing RQs

% Assuming dp already has following RQs:
%
% dp.pulse_start_samples(pulse,evt) - start of pulse range (inclusive), in samples
% dp.pulse_end_samples(pulse,evt)   - end of pulse range (inclusive), in samples

for evt = 1:N
    
    if mod(evt,100) == 0
        fprintf('\n')
    end
    
    fprintf('.'); drawnow();
    
    % Check if fields exist and has data
    if isfield(cvt_struct(evt),'sumpod_data_phe_per_sample') && cvt_struct(evt).empty == 0
        
        loop_max = min([sum(isfinite(dp.index_kept_sumpods(:,evt))) max_num_pulses]);
        % for every POD dd
        for pp = 1:loop_max
            
            % -1 and +1 signs are needed to get the right time cut
            pulse_cut = inrange(cvt_struct(evt).sumpod_time_samples,dp.pulse_start_samples(pp,evt)-1,dp.pulse_end_samples(pp,evt)+1);
            pulse_data_phe = cvt_struct(evt).sumpod_data_phe_per_sample(pulse_cut);
            pulse_time_samples = cvt_struct(evt).sumpod_time_samples(pulse_cut);
            
            [aft_box_open aft_left_edge aft_center aft_right_edge aft_box_close...
               dp.preBoxArea(pp,evt) dp.postBoxArea(pp,evt)] = LUXPerusePeeksMatlab(pulse_data_phe,BoxWidthSamples,edgeFraction);
            
            dp.aft_box_open_samples(pp,evt) = pulse_time_samples(aft_box_open);
            dp.aft_t0_samples(pp,evt) = pulse_time_samples(aft_left_edge);
            dp.aft_t1_samples(pp,evt) = pulse_time_samples(aft_center);
            dp.aft_t2_samples(pp,evt) = pulse_time_samples(aft_right_edge);
            dp.aft_box_close_samples(pp,evt) = pulse_time_samples(aft_box_close);
            
            % Plot to check result
            if 0
               figure(32193); clf
               plot(cvt_struct(evt).sumpod_time_samples,cvt_struct(evt).sumpod_data_phe_per_sample,'k.-'); hold on
               plot(dp.aft_box_open_samples(pp,evt),0,'ro','markers',10);
               plot(dp.aft_t0_samples(pp,evt),0,'rs','markers',10);
               plot(dp.aft_t1_samples(pp,evt),0,'rx','markers',10);
               plot(dp.aft_t2_samples(pp,evt),0,'rs','markers',10);
               plot(dp.aft_box_close_samples(pp,evt),0,'ro','markers',10);
               xlim([dp.pulse_start_samples(pp,evt)-30  dp.pulse_end_samples(pp,evt)+30])
               keyboard
            end
            
        end % for pulse pp
    end % fi
end % for event evt

fprintf('Done!\n')

%% Write Output File

status = LUXBinaryRQWriter_framework(settings, dp, filename_rq, data_path_rq, event_number, livetime); % Should add livetime input at the end, remove event_number as settings field


