function status = PulseFinder_Naive(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)
%
% This pulse finder operates in "naive" mode, simpy giving the POD edges as
% pulses, with a 1:1 conversion from PODs to pulses. This ensures that the
% required RQs are present for other modules to use. Replace this function
% with an actual pulse finder when one is ready for prime time.
%
% status = PulseFinder_Naive(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)
%
%
%
% Required RQs:
%
%
% Versioning
%   20121024 - CHF - Created
%   20121116 - CHF - Changed inputs. It doesn't need settings as inputs
%                    anymore. It grabs them from dp.admin in RQ.
%   20121218 - JRV - Changed signature to use paths for xmls.
%                    Now does not load any DP or IQ settings from admin!
%   20130212 CHF - Big changes.
%                  * Changed pulse_start_timestamp_samples (and end) to 
%                    pulse_start_samples, which now uses trigger-based timing, 
%                    instead of global timestamp.
%                  * Added fields pulse_start_pod, pulse_end_pod
%                  * Changed all references from cvt_struct(evt).chsum.pod() 
%                    to cvt_struct(evt).sumpod(), and fixed new subfield names
%                  * Made sure all calls are to _framework code versions
%   20130216 CHF - Minor fixes.
%   20130304 JRV - Fixed bug #15; added line ~100 check for empty event
%   20130315 JRV - Now skips files with 0 events
%
% RQ versions:
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

myname = 'PulseFinder_Naive';
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

%% Prepare

N = length(cvt_struct);

pulse_event_size = [max_num_pulses N];

dp.pulse_start_samples = nan(pulse_event_size);
dp.pulse_end_samples   = nan(pulse_event_size);
dp.index_kept_sumpods  = nan(pulse_event_size);

settings.evt_settings  = dp.admin.evt_settings;
settings.daq_settings  = dp.admin.daq_settings;


%% Compute timing RQs

for evt = 1:N    
    
    if mod(evt,100) == 0
        fprintf('\n')
    end
    
    fprintf('.'); drawnow();
    
    % Check if fields exist and the event isn't empty
    if isfield(cvt_struct(evt),'sumpod_data_phe_per_sample') && ~cvt_struct(evt).empty
        
        edges = find(diff(cvt_struct(evt).sumpod_time_samples) > 1)';
        ind_start_samples_all = [1 edges+1];
        ind_end_samples_all = [edges length(cvt_struct(evt).sumpod_time_samples)];
        
        pulse_start_samples_all = cvt_struct(evt).sumpod_time_samples(ind_start_samples_all);
        pulse_end_samples_all = cvt_struct(evt).sumpod_time_samples(ind_end_samples_all);
        
        % If we need to trim to keep only max_num_pulses
        M = length(pulse_start_samples_all);
        
        if M > max_num_pulses
            
            ss = cumsum(cvt_struct(evt).sumpod_data_phe_per_sample);
            areas_phe_all = diff([0 ss(ind_end_samples_all)']);

            [areas_phe_sorted index_kept_sumpods_raw] = sort(areas_phe_all,'descend');
            dp.index_kept_sumpods(:,evt) = sort(index_kept_sumpods_raw(1:max_num_pulses));
            
            dp.pulse_start_samples(:,evt) = pulse_start_samples_all(dp.index_kept_sumpods(:,evt));
            dp.pulse_end_samples(:,evt) = pulse_end_samples_all(dp.index_kept_sumpods(:,evt));            
            
        else
            
            dp.index_kept_sumpods(1:M,evt) = 1:M;
            dp.pulse_start_samples(1:M,evt) = pulse_start_samples_all;
            dp.pulse_end_samples(1:M,evt) = pulse_end_samples_all;
            
        end
        
        
    end % fi
    

    
end % for event evt



fprintf('Done!\n')

%% Write Output File

status = LUXBinaryRQWriter_framework(settings, dp, filename_rq, data_path_rq, event_number, livetime); % Should add livetime input at the end, remove event_number as settings field

