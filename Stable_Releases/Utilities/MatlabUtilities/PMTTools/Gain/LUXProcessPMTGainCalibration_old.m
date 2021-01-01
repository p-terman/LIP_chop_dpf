% function LUXProcessPMTGainCalibration(filename,path,local_acquisition_flag)
%
% This function will process the raw built data for PMT gain calibrations and output
% a .pmtcal file with the areas, peak heights and peak times.
%
% Inputs:
%       filename - dataset name (e.g. lux01_20090101T1234)
%       path     - [optional] folder where the dataset lies. Do not include filename in it!
%       local_acquisition_flag - if 1, use local acquisition settings
%       
% 20091002 - CHF, created. Using LUX binary writer code backbone

function LUXProcessPMTGainCalibration(filename,path,local_acquisition_flag)

%% Get settings

if nargin < 3
    local_acquisition_flag = 0;
end

% XML settings
xml_settings = LUXLoadSettings(filename,path);
fprintf('XML Settings Loaded\n')

% Number of events
nb_evts_in_file = xml_settings.nb_evts;
event_list = 1:nb_evts_in_file;

if local_acquisition_flag
    [PMT_serial PMT_bias_V amplification] = LUXGetLocalAcquisitionSettings(filename);
else
    amplification = LUX01LoadAmp(filename,1); %must change how channel is handled - default, 1
end


%%
% Load data
fprintf('Loading %d events\n',numel(event_list))
event_traces = LUXEvent2Trace_old(xml_settings,event_list); %can we change this to custom LUXEventLoader?

%% Compute areas

dt = 10/1.05;
window = 5:40;
pulses = squeeze(-1*event_traces.yy_mV(window,:,:))/amplification; %size [pts ch events]

calibration.areas_mVns = squeeze(sum(pulses,1)).*dt; %size [ch events]
calibration.peak_height_mV = squeeze(max(pulses,[],1));
[temp time_index] = max(pulses,[],1);
calibration.peak_time_ns = squeeze(time_index)*10/1.05;

%% Write Gain Calibration RQs in Binary File

fprintf('Writing Binary RQ File for PMT Gain Calibration\n');

gain_filename = sprintf('%s/%s/%s.pmtcal',path,filename,filename);
bin_fid = fopen(gain_filename, 'wb');

if isempty(bin_fid) || bin_fid<0
    fprintf('Could not open binary file to write\n');
    return
end

%%% Block 1 %%%

% Write Block 1 Header %
first_evt_in_file = event_list(1);

endianness = hex2dec('01020304');
fwrite(bin_fid, endianness, 'uint32');

b1_header_string = 'dataset_name;char;19;first_evt_in_file;uint32;1;nb_evts_in_file;uint32;1;';
b1_header_string_size = length(b1_header_string);
fwrite(bin_fid, b1_header_string_size,'uint16');
fwrite(bin_fid, b1_header_string, 'char');
b1_nb_data_lines = 1;
fwrite(bin_fid, b1_nb_data_lines, 'int32');

% Write Block 1 Data Lines %
fwrite(bin_fid, filename, 'char');
fwrite(bin_fid, first_evt_in_file, 'uint32');
fwrite(bin_fid, nb_evts_in_file, 'uint32');

%%% Block 2: Calibration RQs %%

% Write Block 2 Header %
var_name = fieldnames(calibration); % get field names
nb_RQs = length(var_name);
per_evt = zeros(nb_RQs,1);

for rq=1:nb_RQs
    
    thisvarname = char(var_name(rq));
    var_types{rq} = class(calibration.(thisvarname));
    if size(calibration.(thisvarname),ndims(calibration.(thisvarname))) == nb_evts_in_file
        dimz_end = (ndims(calibration.(thisvarname))-1);
        per_evt(rq) = 1;
    else
        dimz_end = ndims(calibration.(thisvarname));
        per_evt(rq) = 0;
    end
    for dimz=1:dimz_end
        var_size{rq}(dimz) = size(calibration.(thisvarname),dimz);
    end
    var_size_str{rq} = '';
    for dimz=1:dimz_end
        var_size_str{rq} = strcat(var_size_str{rq}, num2str(var_size{rq}(dimz)),',');
    end
    var_size_str{rq} = var_size_str{rq}(1:(end-1));
    if per_evt(rq) == 1
        if strcmp(var_types(rq),'logical')
            var_types{rq} = 'uint8';
        end
        var_header{rq} = sprintf('%s;%s;%s;', var_name{rq}, var_types{rq}, var_size_str{rq});
    end
end

b2_header_string = [var_header{:}];
b2_header_string_size = length(b2_header_string);
fwrite(bin_fid, b2_header_string_size, 'uint16');
fwrite(bin_fid, b2_header_string, 'char');
fwrite(bin_fid, nb_evts_in_file, 'int32');

%%% Write Block 2 Data %%%

for evt=1:nb_evts_in_file

    % Loop over every RQ %
    for rq=1:nb_RQs
        
        thisvarname = char(var_name(rq));
        RQ_data = calibration.(thisvarname);
        if per_evt(rq) == 1

            if length(var_size{rq})==1
                for dimz=1:var_size{rq}(1)
                    fwrite(bin_fid, RQ_data(dimz,evt), var_types{rq});
                end
            else % if multi-dimensional
                for dimz1=1:var_size{rq}(end) % last dimension

                    for dimz2=1:var_size{rq}(end-1) % next to last dimension

                        if length(var_size{rq})==2
                            fwrite(bin_fid, RQ_data(dimz2,dimz1,evt), var_types{rq});
                        else
                            for dimz3=1:var_size{rq}(end-2) % next to next to last dimension - only good for three dimensions + events right now
                                fwrite(bin_fid, RQ_data(dimz3,dimz2,dimz1,evt), var_types{rq});
                            end
                        end

                    end

                end

            end
        end



    end

end

% Set file permissions to 755 for lux user
% [status, result] = unix(['chmod 755 ' gain_filename]);
% if status ~= 0
%     disp(['Warning: unable to set file permissions for ' gain_filename]);
%     fprintf('Message: %s', result);
% end

clear var_size;
clear per_evt;
clear var_header;
clear var_size_str;
clear dimz*;
clear var_name;
clear var_types;
clear RQ_data;
fclose(bin_fid);
