% function amp = LUX01LoadAmp(filename_prefix,amp_ch_out)
% 
% This function loads the total amplification from the LUG settings for the dataset filename_prefix
% 
% Inputs:
%         filename_prefix - dataset name (e.g. lux01_20090328T2032)
%         amp_ch_out      - amplification ch output. NOT THE PMT CHANNEL
% 
% Outputs:
%         amp             - total amplification (preamp * out_height * out_area)
% 

function amp = LUX01LoadAmp(filename_prefix,amp_ch_out)

%% Load Settings

if nargin < 2
    amp_ch_out = 1; %using the 1st amp out by default
end

xmlsettings = XMLReader('defaultLUGSettings.xml');
filename = [char(39) filename_prefix char(39)];
settings.filename = filename_prefix;

daq_acquisition_fields      = {'filename','daq_acquisition_id','source','run_number','daq_version'...
                               'amplification_state_id','PMT_hv_state_id','grid_hv_state_id'};
                           
amplification_fields        = {'preamp','postamp_out1_height','postamp_out1_area',...
                               'postamp_out2_height','postamp_out2_area','postamp_out3_height','postamp_out3_area'};

                           
%% DAQ Acquisition Settings

clear query_string
query_string = {'use lug'};

for ii = 1:numel(daq_acquisition_fields)
        query_string{end+1} = sprintf('select %s from daq_acquisition where filename = %s',daq_acquisition_fields{ii},filename);
end

query_string = String_Replacement_Table(query_string);

[result_cells query_time_LUG url_alive] = LUGQuery(query_string,xmlsettings);

amplification_state_id = str2double(result_cells{7});

%% Amplification

clear query_string result_cells
query_string = {'use lug'};

for ii = 1:numel(amplification_fields)
        query_string{end+1} = sprintf('select %s from amplification_state where amp_state_id = %d',amplification_fields{ii},amplification_state_id);
end

query_string = String_Replacement_Table(query_string);

[result_cells query_time_LUG url_alive] = LUGQuery(query_string,xmlsettings);

name_height = sprintf('out%d_height',amp_ch_out);
name_area = sprintf('out%d_area',amp_ch_out);

settings.amplification.preamp = str2double(result_cells{2});
settings.amplification.postamp.(name_height) = str2double(result_cells{1 + 2*amp_ch_out});
settings.amplification.postamp.(name_area) = str2double(result_cells{2 + 2*amp_ch_out});

amp = settings.amplification.preamp * ...
            settings.amplification.postamp.(name_height) * ...
            settings.amplification.postamp.(name_area);
