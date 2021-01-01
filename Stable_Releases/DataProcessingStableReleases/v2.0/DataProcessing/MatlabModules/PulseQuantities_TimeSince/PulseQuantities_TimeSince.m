function status = PulseQuantities_TimeSince(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)
% 
% This function determines how long it has been since an event with
% full_event_area_phe greater than a threshold value (4 different
% thresholds are used, so 4 RQs are created).
%
% Required RQs: full_evt_area_phe, event_timestamp_samples
%
% Versioning:
%   v1.0 20130716 EP - Created 
%   v1.1 20130912 EP - Eliminated 1e4, 2e5, 4e5, 1e6 phe thresholds; added
%                      8e5 phe threshold (now 1e5, 3e5, 5e5, and 8e5)
%
% RQ versions:
%   v1.0 time_since_1e4_phe
%   v1.0 time_since_1e5_phe
%   v1.0 time_since_2e5_phe
%   v1.0 time_since_3e5_phe
%   v1.0 time_since_4e5_phe
%   v1.0 time_since_5e5_phe
%   v1.0 time_since_1e6_phe
%
%   v1.1 time_since_1e5_phe
%   v1.1 time_since_3e5_phe
%   v1.1 time_since_5e5_phe
%   v1.1 time_since_8e5_phe
%
%% THIS STUFF YOU ALWAYS NEED
% Load .rq file

dp = LUXLoadRQ1s_framework(filename_rq, data_path_rq);

settings.evt_settings = dp.admin.evt_settings;
settings.daq_settings = dp.admin.daq_settings;
settings.filename_prefix = dp.admin.filename_prefix;

event_number = dp.event_number;
livetime = dp.admin.livetime;

%% Bookkeeping - access parameters

myname = 'TimeSince';
fprintf('\n\n *** Starting module %s\n',myname);

dp_settings_xml = XMLReader_framework(data_processing_xml_path);
lug_iqs_xml = XMLReader_framework(iq_xml_path);

module_names = {dp_settings_xml.data_processing_settings.module.module_name};
index_temp = strfind(module_names,myname);
index_module = find(not(cellfun('isempty', index_temp)));

mymodule_settings = dp_settings_xml.data_processing_settings.module(index_module).parameters;

max_num_pulses = dp_settings_xml.data_processing_settings.global.max_num_pulses;

%% Initialize variables

pmt_chs = 1:122;

N = length(dp.event_timestamp_samples);
pulse_event_size = [max_num_pulses N];

dp.time_since_1e5_phe = zeros(1,N,'single');
dp.time_since_3e5_phe = zeros(1,N,'single');
dp.time_since_5e5_phe = zeros(1,N,'single');
dp.time_since_8e5_phe = zeros(1,N,'single');
last_holdoff_time_1e5_phe = dp.event_timestamp_samples(1);
last_holdoff_time_3e5_phe = dp.event_timestamp_samples(1);
last_holdoff_time_5e5_phe = dp.event_timestamp_samples(1);
last_holdoff_time_8e5_phe = dp.event_timestamp_samples(1);
counts = zeros(1,4);

%% Loop per event and assign time_since value for each threshold

for evt = 1:N;
    if counts(1)==0
        dp.time_since_1e5_phe(evt)=1;
    else
        dp.time_since_1e5_phe(evt) = dp.event_timestamp_samples(evt)-last_holdoff_time_1e5_phe;
    end
    if counts(2)==0
        dp.time_since_3e5_phe(evt)=1;
    else
        dp.time_since_3e5_phe(evt) = dp.event_timestamp_samples(evt)-last_holdoff_time_3e5_phe;
    end
    if counts(3)==0
        dp.time_since_5e5_phe(evt)=1;
    else
        dp.time_since_5e5_phe(evt) = dp.event_timestamp_samples(evt)-last_holdoff_time_5e5_phe;
    end
    if counts(4)==0
        dp.time_since_8e5_phe(evt)=1;
    else
        dp.time_since_8e5_phe(evt) = dp.event_timestamp_samples(evt)-last_holdoff_time_8e5_phe;
    end

    if dp.full_evt_area_phe(evt)>100000
        counts(1) = 1;
        last_holdoff_time_1e5_phe = dp.event_timestamp_samples(evt);
         
        if dp.full_evt_area_phe(evt)>300000
            counts(2) = 1;
            last_holdoff_time_3e5_phe = dp.event_timestamp_samples(evt);
                    
            if dp.full_evt_area_phe(evt)>500000
                counts(3) = 1;
                last_holdoff_time_5e5_phe = dp.event_timestamp_samples(evt);
                
                if dp.full_evt_area_phe(evt)>800000
                    counts(4) = 1;
                    last_holdoff_time_8e5_phe = dp.event_timestamp_samples(evt);
                end
            end
        end
    end

end


%% Write Output File

status = LUXBinaryRQWriter_framework(settings, dp, filename_rq, data_path_rq, event_number, livetime); % Should add livetime input at the end, remove event_number as settings field
