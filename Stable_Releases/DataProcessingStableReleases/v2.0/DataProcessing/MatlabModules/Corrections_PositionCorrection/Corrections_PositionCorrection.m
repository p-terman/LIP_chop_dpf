function status = Corrections_PositionCorrection(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)
% Corrections_PositionCorrection module. It corrects the positions obtained
% using the position reconstruction of the S2 signal.
%
% status = Corrections_PositionCorrection(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)
%
% Required RQs: s1s2_pairing, z_drift_samples,
% pulse_classification, (x_cm, y_cm) or (x_cm_tmplt, y_cm_tmplt)
%
% RQ versions:
% 1.0 Initial Release
% 1.1 20130729 - CFPS A translation along the plane XY was added to ensure that each quadrant of the chamber has the same number of events for the Krypton callibration
% 1.10.1 20130731 - CFPS A patch related with the name of the module
% 1.2 20130903 - CFPS Corrects both TAXY and Mercury reconstruction
% 1.3.1 % 20140723 CFPS - Variable type changed from double to float

module_version = 1.2;

%% Load .rq file

status = [];

dp = LUXLoadRQ1s_framework(filename_rq, data_path_rq);

%% Definition of the settings
if isfield(dp.admin, 'evt_settings')
    settings.evt_settings = dp.admin.evt_settings;%settings.evt_settings = dp.admin.settings.event_builder_settings;
end
if isfield(dp.admin, 'daq_settings')  
    settings.daq_settings = dp.admin.daq_settings;
end
settings.filename_prefix = dp.admin.filename_prefix;
event_number = dp.event_number;

%% Read the LUG settings and the table with the corrections

myname = 'Corrections_PositionCorrection';
disp(sprintf('\n *** Starting module %s - version %.2f***\n',myname, module_version));

if isempty(event_number)
    fprintf('\n\n *** Skipping module (no events in file) %s\n',myname);
    return
end

dp_settings_xml = XMLReader_framework(data_processing_xml_path);
lug_iqs_xml = XMLReader_framework(iq_xml_path);

for ii = 1:length(lug_iqs_xml.iq)
    if isfield(lug_iqs_xml.iq(ii).global, 'iq_type') 
        if strcmp(lug_iqs_xml.iq(ii).global.iq_type,'xy_rec_cor')
            lrf_iq = lug_iqs_xml.iq(ii);
        end
    end
end

if ~isstruct(lrf_iq)
    disp(fprintf('\n\n %s: The IQ xy_rec_cor was not found\n\n************* Fatal Error*************\n\n',myname));
    return
end

position_correction_path = which('Corrections_PositionCorrection');
 
if ~isfield(lrf_iq, 'file_table_Mercury') & isfield(lrf_iq, 'file_table')
    lrf_iq.file_table_Mercury = lrf_iq.file_table;
end

if isfield(lrf_iq, 'file_table_Mercury')
    position_correction_table_path = [position_correction_path(1:(end-numel(myname)-2)) lrf_iq.file_table_Mercury];
    if ~any(position_correction_table_path)
        disp(fprintf('\n\n %s: The file %s was not found in %s\n\n*************\n\n',myname, lrf_iq.file_table, position_correction_path(1:(end-numel(myname)-2))));
    else
        table_corrections = load(position_correction_table_path);
        %% Run the function. 
        dp = Corrections_PositionCorrection_Function(dp, table_corrections);
        
    end
    dp.x_corrected = single(dp.x_corrected);
    dp.y_corrected = single(dp.y_corrected);

end

%{
 if isfield(lrf_iq, 'file_table_TAXY')
    position_correction_table_path = [position_correction_path(1:(end-numel(myname)-2)) lrf_iq.file_table_TAXY];
    if ~any(position_correction_table_path)
        disp(fprintf('\n\n %s: The file %s was not found in %s\n\n*************\n\n',myname, lrf_iq.file_table, position_correction_path(1:(end-numel(myname)-2))));
    else
        table_corrections = load(position_correction_table_path);
        %% Run the function. 
        dp = Corrections_PositionCorrection_Function(dp, table_corrections);        
    end
    dp.x_tmplt_corrected = single(dp.x_tmplt_corrected);
    dp.y_tmplt_corrected = single(dp.y_tmplt_corrected);
    dp.xy_sigma_corrected = single(dp.xy_sigma_corrected);
end
%}

livetime = LUXGetLivetime_framework(filename_evt,data_path_evt);
status = LUXBinaryRQWriter_framework(settings, dp, filename_rq, data_path_rq, dp.event_number, livetime);


