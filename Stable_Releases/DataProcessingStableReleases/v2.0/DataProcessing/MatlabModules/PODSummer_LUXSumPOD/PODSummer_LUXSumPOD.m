function status = PODSummer_LUXSumPOD(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)
%
%
%
% Required RQs:
%
%
% Versioning
%
%
% RQ versions:
%
% 20130315 JRV - Now skips files with 0 events
%
%% Load .rq file

status = [];

myname = 'PODSummer_LUXSumPOD';
fprintf('\n\n *** Starting module %s\n',myname);

dp = LUXLoadRQ1s_framework(filename_rq, data_path_rq);

event_number = dp.event_number;
livetime = dp.admin.livetime;

if isempty(event_number)
    fprintf('\n\n *** Skipping module (no events in file) %s\n',myname);
    return
end

%% Load .cvt file

[cvt_struct settings] = LUXCVTLoader_framework(data_path_evt,strrep(filename_evt,'evt','cvt'));

%% grab module settings
dp_settings_xml = XMLReader_framework(data_processing_xml_path);
module_names = {dp_settings_xml.data_processing_settings.module.module_name};
index_temp = strfind(module_names,myname);
index_module = find(not(cellfun('isempty', index_temp)));
mymodule_settings = dp_settings_xml.data_processing_settings.module(index_module).parameters;
%cvt_struct(1).thr = mymodule_settings.flatten_thr; % just stick it in cvt_struct(1)

%% Calculate Ch Sum
tic
fprintf('Summing Pulses... ')
cvt_struct = LUXSumPOD_framework(cvt_struct);
toc

%% trust, but verify
if 0
   figure(4);%clf;
   plot(cvt_struct(1).sumpod_time_samples,cvt_struct(1).sumpod_data_phe_per_sample,'r-'); 
   hold on;
   %plot(cvt_struct(1).sumpod_time_samples,cvt_struct(1).sumpod_data_thr_phe_per_sample,'r-'); 
   keyboard;
end

%% Write Output File

filename_cvt_fullpath = [data_path_evt filesep strrep(filename_evt,'evt','cvt')];

fprintf('Saving CVT file... ');
status = LUXCVTWriter_framework( cvt_struct, settings, livetime, filename_cvt_fullpath );
fprintf('Done\n');


