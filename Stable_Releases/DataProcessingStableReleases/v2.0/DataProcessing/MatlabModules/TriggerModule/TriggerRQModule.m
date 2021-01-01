function status = TriggerRQModule(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)
% 
% This module propagates all trigger xlm information stored in evt file to
% rqs which also convert them to more user-friendly form. This module also
% looks at the tte waveform to determine the type of the trigger associated
% to the events.
%
%
% Required RQs: none, but need to load evt file once again.
%
% Versioning:
%   v1.0 20131213 MM - Initial released
%
% RQ versions:
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

myname = 'TriggerRQModule';
fprintf('\n\n *** Starting module %s\n',myname);

dp_settings_xml = XMLReader_framework(data_processing_xml_path);
lug_iqs_xml = XMLReader_framework(iq_xml_path);

module_names = {dp_settings_xml.data_processing_settings.module.module_name};
index_temp = strfind(module_names,myname);
index_module = find(not(cellfun('isempty', index_temp)));

mymodule_settings = dp_settings_xml.data_processing_settings.module(index_module).parameters;

max_num_pulses = dp_settings_xml.data_processing_settings.global.max_num_pulses;

%% Load events
Debug.Mode = dp_settings_xml.data_processing_settings.global.debug_mode;
Debug.tok = 0;
Debug.PauseTime = 1;
options.load_tte_ch = 1;
options.load_xlm_info = 1;
options.load_xenon_chs = 1;
Debug.tok = fprintf(1, '%s\n%s',char(8*ones(1,Debug.tok)),'Loading evt file.                                       ') - Debug.tok;
event_struct = LUXEventLoader_framework(data_path_evt, filename_evt, options);
Debug.tok = fprintf(1, '\n%s\n%s',char(8*ones(1,Debug.tok)),'                                                      ') - Debug.tok;
Debug.tok = fprintf(1, '%s\n%s',char(8*ones(1,Debug.tok)),'Finish loading evt file.                                ') - Debug.tok;
TotEvent = length(event_struct);

%% Load dp and initialize the RQs
if TotEvent > 0
    dp.xlm_timestamp_samples                    = zeros(1,TotEvent,'int64');
    dp.xlm_trigger_number                       = zeros(1,TotEvent,'int32');
    dp.xlm_max_filter_output_value_mVns         = zeros(1,TotEvent,'single');
    dp.xlm_max_filter_output_channel            = zeros(1,TotEvent,'int32');
    dp.xlm_s1_hit_vector                        = zeros(16,TotEvent,'int32');
    dp.xlm_s2_hit_vector                        = zeros(16,TotEvent,'int32');

    dp.tte_sat_value_mV                         = zeros(1,TotEvent,'single');
    dp.tte_sat_flatness_mV                      = zeros(1,TotEvent,'single');
    dp.tte_type                                 = zeros(1,TotEvent,'int8');
    dp.tte_comment                              = zeros(1,TotEvent,'int32');
    if Debug.Mode
        Debug.tok = fprintf(1, '%s\n%s',char(8*ones(1,Debug.tok)),'Memory preallocated.                                    ') - Debug.tok;
        pause(Debug.PauseTime)
    end
else
    if Debug.Mode
        Debug.tok = fprintf(1, '%s\n%s',char(8*ones(1,Debug.tok)),'Total event is zero.                                    ') - Debug.tok;
        pause(Debug.PauseTime)
    end
end
%% Load other reference
if Debug.Mode
    Debug.tok = fprintf(1, '%s\n%s',char(8*ones(1,Debug.tok)),'Loading reference files.                                ') - Debug.tok;
    pause(Debug.PauseTime)
end
TrgRef = load(strrep(which(myname),[myname '.m'],'TTEOutputTemplate.mat'));
if Debug.Mode
    Debug.tok = fprintf(1, '%s\n%s',char(8*ones(1,Debug.tok)),'Reference files loaded.                                 ') - Debug.tok;
    pause(Debug.PauseTime)
end
%% Calculate xlm and tte

for evt = 1:TotEvent
    if Debug.Mode
        Debug.tok = fprintf(1, '%s\n%s',char(8*ones(1,Debug.tok)),'Calling XLM_TTE_Decoder.                                ') - Debug.tok;
        pause(Debug.PauseTime)
    end
    [xlm tte] = XLM_TTE_Decoder(event_struct(evt),TrgRef.TTETemplate,vertcat(event_struct.timestamp),options,Debug);
    if xlm.success == 1
        dp.xlm_timestamp_samples(1,evt)                = xlm.xlm_timestamp_samples;
        dp.xlm_trigger_number(1,evt)                   = xlm.trigger_number;
        dp.xlm_max_filter_output_value_mVns(1,evt)     = xlm.max_filter_output_value_mVns;
        dp.xlm_max_filter_output_channel(1,evt)        = xlm.max_filter_output_channel;
        dp.xlm_s1_hit_vector(:,evt)                    = xlm.s1_hit_vector;
        dp.xlm_s2_hit_vector(:,evt)                    = xlm.s2_hit_vector;
    else
        dp.xlm_timestamp_samples(1,evt)                = -999999;
        dp.xlm_trigger_number(1,evt)                   = -999999;
        dp.xlm_max_filter_output_value_mVns(1,evt)     = -999999;
        dp.xlm_max_filter_output_channel(1,evt)        = -1;
        dp.xlm_s1_hit_vector(:,evt)                    = -1;
        dp.xlm_s2_hit_vector(:,evt)                    = -1;
    end
    if tte.success == 1
        dp.tte_sat_value_mV(1,evt)                     = tte.tte_sat_value_mV;
        dp.tte_sat_flatness_mV(1,evt)                  = tte.tte_sat_flatness_mV;
        dp.tte_type(1,evt)                             = tte.tte_type;
        dp.tte_comment(1,evt)                          = tte.tte_comment;
    else
        dp.tte_sat_value_mV(1,evt)                     = -999999;
        dp.tte_sat_flatness_mV(1,evt)                  = -999999;
        dp.tte_type(1,evt)                             = -6;
        dp.tte_comment(1,evt)                          = -6;
    end
    if Debug.Mode
        Debug.tok = fprintf(1, '%s\n%s',char(8*ones(1,Debug.tok)),['RQ values of event ' num2str(evt,'%d') ' assigned. ']) - Debug.tok;
        pause(Debug.PauseTime)
    else
        Debug.tok = fprintf(1, '%s\nEvents completed : [%d/%d]',char(8*ones(1, Debug.tok)), evt, TotEvent) - Debug.tok;
    end
end

%% Write Output File
status = LUXBinaryRQWriter_framework(settings, dp, filename_rq, data_path_rq, event_number, livetime); % Should add livetime input at the end, remove event_number as settings field









