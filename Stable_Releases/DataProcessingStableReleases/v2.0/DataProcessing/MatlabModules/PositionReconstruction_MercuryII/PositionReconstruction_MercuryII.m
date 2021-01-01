function status = PositionReconstruction_MercuryII(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)
% 
% This is a PositionReconstruction module using the MercuryII algorithm.
%
% status = PositionReconstruction_MercuryII(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)
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
% 20121025 CHF - Created
% 20121218 JRV - Changed signature to use paths for xmls.
%                  Now does not load any DP or IQ settings from admin!
% 20130120 CHF - Minor changes. Should now work with new LRF format in
%                  LUG (2-dim params, instead of 7 deg polynomial).
% 20130212 CHF - Made sure all calls are to _framework code versions
% 20130315 JRV - Now skips files with 0 events
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

myname = 'PositionReconstruction_MercuryII';
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

dp.x_cm = nan(pulse_event_size);
dp.y_cm = nan(pulse_event_size);

top = [1:60 121];


livetime = LUXGetLivetime_framework(filename_evt,data_path_evt);

%% Prepare LRFs

% Figure out which iq has the pmt_gains - assuming only ONE iq was returned
% for each type
for ii = 1:length(lug_iqs_xml.iq)
    if strcmp(lug_iqs_xml.iq(ii).global.iq_type,'lrfs') == 1

        LRF_matrix.xx = lug_iqs_xml.iq(ii).LRF.channel(1).data_radius_vector_cm;

        M = length(LRF_matrix.xx);
        LRF_matrix.yy = zeros(122,M);
        LRF_matrix.params = zeros(122,2);
        
        qq = 1;
        for pp = top
           
            if pp == lug_iqs_xml.iq(ii).LRF.channel(qq).PMT_number
                
                LRF_matrix.yy(pp,:) = lug_iqs_xml.iq(ii).LRF.channel(qq).data_LRF_vector;
                LRF_matrix.params(pp,:) = lug_iqs_xml.iq(ii).LRF.channel(qq).params;
            else
                fprintf('LRF indexing mismatch!\n');
            end
            
            qq = qq + 1; % this index needs to go from 1:61 (instead of 121)
        end
    
        
        break
    end
end

%% Perform position reconstruction

for ii = 1:max_num_pulses
    
    signal = squeeze(dp.peak_area_phe(ii,:,:));
    [dp.x_cm(ii,:) dp.y_cm(ii,:)] = LUXMercuryPositionReconstruction_framework(signal,LRF_matrix);
    
end

%% Write Output File

status = LUXBinaryRQWriter_framework(settings, dp, filename_rq, data_path_rq, dp.event_number, livetime); % Should add livetime input at the end, remove event_number as settings field



