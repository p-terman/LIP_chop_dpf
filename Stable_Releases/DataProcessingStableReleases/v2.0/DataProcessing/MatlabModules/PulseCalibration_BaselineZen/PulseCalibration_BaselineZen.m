function status = PulseCalibration_BaselineZen(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)
%
%
%
% Required RQs:
%
%
% Versioning
%   20130212 - Created
%   20130305 - Added error trapping for IQs that don't have appropriate
%              fields: iq.global.iq_type
%   20130328 pfs - recording adc rqs from BaselineZen
%   20130315 JRV - Now skips files with 0 events
%
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
amp_gain = dp.admin.daq_settings.global.preamp .* dp.admin.daq_settings.global.postamp;

%% Bookkeeping

myname = 'PulseCalibration_BaselineZen';
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

if ~isempty(index_module)
    dp_settings_xml.data_processing_settings.module(index_module);
    mymodule_settings = dp_settings_xml.data_processing_settings.module(index_module).parameters;
else
%     error(sprintf('*** ERROR: Module was not found in settings file:\n%s\n',data_processing_xml_path));
end

max_num_pulses = dp_settings_xml.data_processing_settings.global.max_num_pulses;

%% Load .evt file, event list and livetime

event_struct = LUXEventLoader_framework(data_path_evt,filename_evt); % CAN WE STANDARDIZE THE ORDER OF INPUTS PLEASE??? CHF

%% Calibrate Pulses

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

% THIS NEEDS UPDATING ??
fprintf('\nGiving Zen to Baselines... ')
event_struct = LUXBaselineZen_framework(event_struct);
% here we record the rqs from LUXBaselineZen
dp.adc_ppe = zeros(1,length(dp.event_number),'uint8');
dp.adc_sds = zeros(1,length(dp.event_number),'uint8');
dp.zen_applied = zeros(1,length(dp.event_number),'uint8');
for evt=1:length(dp.event_number)
	dp.adc_ppe(:,evt) = any(event_struct(evt).adc_ppe);
	dp.adc_sds(:,evt) = any(event_struct(evt).adc_sds);
	dp.zen_applied(:,evt) = any(event_struct(evt).zen_applied);
end

%% Write output rq file since we have defined 3 rqs at this point

status = LUXBinaryRQWriter_framework(settings, dp, filename_rq, data_path_rq, dp.event_number, livetime); % Should add livetime input at the end, remove event_number as settings field

fprintf('Calibrating Pulses... ')
event_struct(1).thr = mymodule_settings.flatten_thr; % simple way to port this value into the function
event_struct = LUXCalibratePulses_framework(event_struct,pmt_gains_mVns_per_phe,amp_gain);

%% Write Output File

filename_cvt_fullpath = [data_path_evt filesep strrep(filename_evt,'evt','cvt')];

fprintf('Saving CVT file... ');
status = LUXCVTWriter_framework( event_struct, settings, livetime, filename_cvt_fullpath );
fprintf('Done\n');


