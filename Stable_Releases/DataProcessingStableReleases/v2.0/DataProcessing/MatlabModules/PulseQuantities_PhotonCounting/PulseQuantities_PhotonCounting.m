function status = PulseQuantities_PhotonCounting(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)
%
% This is a PulseQuantities module that calculates the Photon
% Counting quantitiesin Matlab
%
% status = PulseQuantities_PhotonCounting(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)
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
%
%   v1.0 20140415 CFPS - Photon counting created from the PulseQuantities_MinimumSet
%   v1.1 20140421 CFPS - The module can use peak_height_mV or peak_height_phe
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

%% Bookkeeping

myname = 'PulseQuantities_PhotonCounting';
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

%% Load .cvt file or calculate it

filename_cvt = strrep(filename_evt,'evt','cvt');

if ~exist([data_path_evt filesep filename_cvt],'file')
    fprintf('Did not find .cvt file. Running Summer Module\n');
    status = PulseCalibration_BaselineZen(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path);
    status = PODSummer_LUXSumPOD(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path);
end

fprintf('Loading sum from %s\n',filename_cvt);
[cvt_struct settings] = LUXCVTLoader_framework(data_path_evt,strrep(filename_evt,'evt','cvt'));

tic
disp(tic)
% We still need the mV data...
event_struct = LUXEventLoader_framework(data_path_evt, filename_evt);
disp(toc)
disp(tic)

%% Initialize variables

pmt_chs = 1:122;

N = length(cvt_struct);

pulse_event_size = [max_num_pulses N];
per_channel_event_size = [max_num_pulses length(pmt_chs) N];

dp.peak_width_samples                   = zeros(per_channel_event_size,'uint32');
dp.spike_count                          = zeros(per_channel_event_size,'uint32');
% new RQs --pfs
%% Compute RQs

fprintf('\n');

dp.event_timestamp_samples = [event_struct(:).timestamp];

for evt = 1:N
    
    if mod(evt,100) == 1
        fprintf('\n')
    end
    
    fprintf('.');
    
    % Check if fields exist and has data
    if isfield(cvt_struct(evt),'sumpod_data_phe_per_sample') && cvt_struct(evt).empty == 0
        loop_max = min([sum(isfinite(dp.index_kept_sumpods(:,evt))) max_num_pulses]);

        % for every pulse (NOT pod)
        for pp = 1:loop_max
            
            
            % PER CH QUANTITIES
            ch_map = find(~[cvt_struct(evt).ch(:).empty]);
            ch_map = intersect(ch_map,find(~[event_struct(evt).ch(:).empty])); %need the phe and the mV scales
            if ~isempty(ch_map)
                
                for ch = ch_map
                    
                    % Make slice of data
                    peak_cut = inrange(cvt_struct(evt).ch(ch).pod_time_samples,dp.pulse_start_samples(pp,evt)-1,dp.pulse_end_samples(pp,evt)+1);
                    
                    if ~isempty(peak_cut) && sum(peak_cut) > 0
                        
                        if strcmp(mymodule_settings.height_units, 'mV') == 1                          
                            peak_data_height = event_struct(evt).ch(ch).pod_data_mV(peak_cut);
                        elseif strcmp(mymodule_settings.height_units, 'phe') == 1
                            peak_data_height = cvt_struct(evt).ch(ch).pod_data_phe_per_sample(peak_cut);
                        end
                        peak_time_samples = cvt_struct(evt).ch(ch).pod_time_samples(peak_cut);

                        % This is to find which baseline_mV to use. Since
                        % we are cutting pods in time, it's not
                        % straightforward to know which pod baseline 'twas                            
                        % But this statement takes care of that
                        peak_inds = unique(cumsum(ismember(cvt_struct(evt).ch(ch).pod_time_samples,cvt_struct(evt).ch(ch).pod_start_samples)));
                        
                        % new RQs related with the photon couting module
                        %test 0
                        
                        peak_above_threshold = peak_time_samples(peak_data_height>mymodule_settings.threshold);
                        
                        if ~isempty(peak_above_threshold)
                            if numel(peak_data_height) > 2
                                dp.peak_width_samples(pp,ch,evt) = max(peak_above_threshold)-min(peak_above_threshold)+1;
                                peak_height_phe = peak_data_height>mymodule_settings.threshold;
                                Variacao = (peak_data_height(2:end)-peak_data_height(1:end-1)).*peak_height_phe(1:end-1);
                                dp.spike_count(pp,ch,evt) = sum(Variacao(2:end)<0 & Variacao(1:end-1)>=0);
                            end
                        end
                    
                    end % fi
                end % for ch
            end % fi            
        end % for pulse pp
    end % fi sumpod
end % for event evt

fprintf('Done!\n')

%% Write Output File

status = LUXBinaryRQWriter_framework(settings, dp, filename_rq, data_path_rq, event_number, livetime); % Should add livetime input at the end, remove event_number as settings field
disp(toc)

