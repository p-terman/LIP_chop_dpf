% function channels = LUXGetChannels(filename_prefix)
% 
% This function outputs the channels used for the acquisition filename_prefix
% 
% Inputs:
%         filename_prefix - dataset name (e.g. lux01_20090328T2032)
% 
% Outputs:
%         channels.serial - goes for 1 to number of channels used in acquisition. Each
%                           entry containes the PMT serial number of that channel. If the
%                           channel was off, the serial will be empty
% 
%
%

function channels = LUXGetChannels(filename_prefix)

%% Load Settings

xmlsettings = XMLReader('defaultLUGSettings.xml');
filename = [char(39) filename_prefix char(39)];
settings.filename = filename_prefix;

daq_acquisition_fields      = {'filename','daq_acquisition_id','source','run_number','daq_version'...
                               'amplification_state_id','PMT_hv_state_id','grid_hv_state_id'};
                           
PMT_ch_assignment_fields    = {'ch1','ch2','ch3','ch4'};


%% DAQ Acquisition Settings

clear query_string
query_string = {'use lug'};

for ii = 1:numel(daq_acquisition_fields)
        query_string{end+1} = sprintf('select %s from daq_acquisition where filename = %s',daq_acquisition_fields{ii},filename);
end

query_string = String_Replacement_Table(query_string);

[result_cells query_time_LUG url_alive] = LUGQuery(query_string,xmlsettings);

settings.run_number = str2double(result_cells{5});

%% PMT Channel Assignment ID

clear query_string result_cells
query_string = {'use lug'};

query_string{end+1} = sprintf('select PMT_ch_assignment_id from runs where run_number = %d',settings.run_number);
query_string = String_Replacement_Table(query_string);
[result_cells query_time_LUG url_alive] = LUGQuery(query_string,xmlsettings);

PMT_ch_assignment_id = str2double(result_cells{2});

%% PMT Channel Assignments

clear query_string result_cells
query_string = {'use lug'};

for ii = 1:numel(PMT_ch_assignment_fields)
        query_string{end+1} = sprintf('select %s from lux01_PMT_ch_assignments where ch_assignment_id = %d',PMT_ch_assignment_fields{ii},PMT_ch_assignment_id);
end

query_string = String_Replacement_Table(query_string);

[result_cells_temp query_time_LUG url_alive] = LUGQuery(query_string,xmlsettings);

%length(result_cells_temp)

for ii = 2:length(result_cells_temp)
    result_cells{ii-1} = result_cells_temp{ii};
end;

m = length([result_cells]);
ch_ind = zeros(1,m);

for ii = 1:m
    if ~isempty(result_cells{ii})
        ch_ind(ii) = true;
    end
end

number_channels = sum(ch_ind);


for num = 1:m
    if m
        serial = result_cells{num};
        channels(num).serial = serial(1:end-2);
    else
        channels(num).serial = [];
    end
end
