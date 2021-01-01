function status = PulseQualityCheck_Basic(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)
%
% status = PulseQualityCheck_Basic(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)
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
%   20121115 CHF - Created
%   20130212 CHF - Made sure all calls are to _framework code versions
%   20130315 JRV - Now skips files with 0 events
%
%
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

%% Bookkeeping

myname = 'PulseQualityCheck_Basic';
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


%% Initialize variables

pmt_chs = 1:122; % must take from daq settings

N = length(dp.event_number);

pulse_event_size = [max_num_pulses N];

dp.negative_to_positive_area_ratio = zeros(pulse_event_size);
dp.pulse_quality_flag = zeros(pulse_event_size);

livetime = LUXGetLivetime_framework(filename_evt,data_path_evt);

%% For now

dp.negative_to_positive_area_ratio = dp.pulse_area_negative_phe ./ dp.pulse_area_positive_phe;

% must use double() since logicals don't play well with binary writer
dp.pulse_quality_flag = double(dp.negative_to_positive_area_ratio <= mymodule_settings.negative_to_positive_area_ratio);

%% Write Output File

status = LUXBinaryRQWriter_framework(settings, dp, filename_rq, data_path_rq, dp.event_number, livetime); % Should add livetime input at the end, remove event_number as settings field


