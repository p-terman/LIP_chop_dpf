function status = PulseClassifier_TransparentCornea(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)
%
% status = PulseClassifier_TransparentCornea(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)
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
%   20130618 pfs - Created by copying template from PulseClassifier_S2FilterSpace
%                  and dumping in a collage of selection criteria
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

myname = 'PulseClassifier_TransparentCornea';
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

%mymodule_settings = dp_settings_xml.data_processing_settings.module(index_module).parameters;
livetime = LUXGetLivetime_framework(filename_evt,data_path_evt);


%% Initialize variables

%max_num_pulses = dp_settings_xml.data_processing_settings.global.max_num_pulses;
max_num_pulses = size(dp.peak_area_phe,1);
N = size(dp.peak_area_phe,3);
pulse_event_size = [max_num_pulses N];

% Initialize as 0
dp.pulse_classification_tc = zeros(pulse_event_size);


%% Categories

% 0 means no pulse (empty entry in matrix)
% 10 means primary S1
% 11 means secondary S1
% 20 means primary S2
% 21 means secondary S2
% 30 means single phe
% 40 means single electron S2
% 5 means none of the above (nota)



%% decisions

% NOTE: all params will need to be dumped into xml eventually, for now can hard-code

dp.coin = squeeze(sum(dp.skinny_peak_area_phe>0.5,2)); % coincidence for S1

s2mask = ... 
		dp.amis1_fraction<=0.5 ... maximum box area ratio
		& (dp.aft_t2_samples-dp.aft_t0_samples)>50 ... minimum width
		& squeeze(sum(dp.peak_area_phe,2))>5 ... minimum area
		;
dp.pulse_classification_tc = s2mask.*21; % label all S2 secondary
[s2p,c2] = max(squeeze(sum(dp.peak_area_phe,2)) .* s2mask);
for ii=1:size(dp.pulse_classification_tc,2);
	dp.pulse_classification_tc(c2(ii),ii) = 20; % overwrite label for the primary S2 
end

s1mask = ...
		dp.amis1_fraction>0.5  ... minimum box area ratio
		...& (dp.aft_t2_samples-dp.aft_t0_samples)<=100 ... above line effectively checks width, no need to repeat?
		& (dp.aft_t1_samples-dp.aft_t0_samples)<=10 ... maximum aft risetime
		& dp.coin>=2 ... minimum coincidence in skinny box
		;
[s1p,c1] = max(squeeze(sum(dp.peak_area_phe,2)) .* s1mask);
for ii=1:size(dp.pulse_classification_tc,2);
	dp.pulse_classification_tc(c1(ii),ii) = 10; % label the primary S1 
end



%% Write Output File

status = LUXBinaryRQWriter_framework(settings, dp, filename_rq, data_path_rq, dp.event_number, livetime); % Should add livetime input at the end, remove event_number as settings field

