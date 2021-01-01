function status = InitializeRQFile_Default(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)
% This is a InitializeRQFile modules. It produces an rq file with evt_num.
% status = InitializeRQFile_Default(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)
% 
% status = InitializeRQFile_Default(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)
%
% Required RQs:
%   None
%
% Versioning
%   20121024 CHF - Created
%   20121218 JRV - Changed signature to use paths for xmls
%                    DP and IQ settings are no longer written to admin!
%   20130212 CHF - Made sure all calls are to _framework code versions
%
%
% RQ versions:
%
%%

myname = 'InitializeRQFile_Default';
fprintf('\n\n *** Starting module %s\n',myname);

filename_prefix = filename_evt(1:19);

[settings, ~] = LUXSuperLoader_framework(filename_evt,data_path_evt);
    
settings.filename_prefix = filename_prefix;

%%

dp.event_number = LUXNumberEvents_framework(filename_evt,data_path_evt);
dp.file_number = ones(size(dp.event_number))'*str2double(filename_evt(22:27));
livetime = LUXGetLivetime_framework(filename_evt,data_path_evt);

%% Write Output File

status = LUXBinaryRQWriter_framework(settings, dp, filename_rq, data_path_rq, dp.event_number, livetime); % Should add livetime input at the end, remove event_number as settings field


