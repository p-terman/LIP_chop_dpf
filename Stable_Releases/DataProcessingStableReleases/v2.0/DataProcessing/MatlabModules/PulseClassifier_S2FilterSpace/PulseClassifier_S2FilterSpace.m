function status = PulseClassifier_S2FilterSpace(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)
%
% status = PulseClassifier_S2FilterSpace(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)
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

myname = 'PulseClassifier_S2FilterSpace';
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

pmt_chs = 1:122;

N = length(dp.event_number);
pulse_event_size = [max_num_pulses N];

% Initialize as none of the above (Category 5)
dp.pulse_classification = nan(pulse_event_size);

livetime = LUXGetLivetime_framework(filename_evt,data_path_evt);

%% Categories

% 0 means no pulse (empty entry in matrix)
% 1 means S1
% 2 means S2
% 3 means single phe
% 4 means single electron S2
% 5 means none of the above (nota)

cut_no_pulse = dp.pulse_area_phe == 0;

cut_s1 = dp.s2filter_max_area_diff./dp.pulse_area_phe <= mymodule_settings.fractional_s2filter_division;
cut_s2 = dp.s2filter_max_area_diff./dp.pulse_area_phe > mymodule_settings.fractional_s2filter_division;

% These are placeholders
cut_sphe = zeros(pulse_event_size);
cut_single_electron_s2 = zeros(pulse_event_size);
cut_nota = ~cut_s1 & ~cut_s2 & ~cut_no_pulse; % (none of the above)

% Put it all together. Remaining zeros means there is no pulse in that location.
dp.pulse_classification = double(cut_s1) + 2*double(cut_s2) + ...
                        3*double(cut_sphe) + 4*double(cut_single_electron_s2) + ...
                        5*double(cut_nota);

%% Write Output File

status = LUXBinaryRQWriter_framework(settings, dp, filename_rq, data_path_rq, dp.event_number, livetime); % Should add livetime input at the end, remove event_number as settings field

