%{ 

function amp = LUX01LoadPMTHV(filename_prefix,ch)

This function loads the PMT HV value from the LUG settings for the channel ch in dataset filename_prefix

Inputs:
        filename_prefix - dataset name (e.g. lux01_20090328T2032)
        ch              - PMT channel (starts at 1) for finding hv bias value

Outputs:
        hv             - PMT HV value for given channel ch (remember: biases are negative!)


%}

function hv = LUX01LoadPMTHV(filename_prefix,ch)

%% Load Settings

xmlsettings = XMLReader('defaultLUGSettings.xml');
filename = [char(39) filename_prefix char(39)];
settings.filename = filename_prefix;

daq_acquisition_fields      = {'filename','daq_acquisition_id','source','run_number','daq_version'...
                               'amplification_state_id','PMT_hv_state_id','grid_hv_state_id'};
                           
PMT_hv_fields               = {'bias_ch1','bias_ch2','bias_ch3','bias_ch4'};

                           
%% DAQ Acquisition Settings

clear query_string
query_string = {'use lug'};

for ii = 1:numel(daq_acquisition_fields)
        query_string{end+1} = sprintf('select %s from daq_acquisition where filename = %s',daq_acquisition_fields{ii},filename);
end

query_string = String_Replacement_Table(query_string);

[result_cells query_time_LUG url_alive] = LUGQuery(query_string,xmlsettings);

PMT_hv_state_id = str2double(result_cells{8});

%% PMT HV

clear query_string result_cells
query_string = {'use lug'};

for ii = 1:numel(PMT_hv_fields)
        query_string{end+1} = sprintf('select %s from lux01_PMT_hv_state where hv_state_id = %d',PMT_hv_fields{ii},PMT_hv_state_id);
end

query_string = String_Replacement_Table(query_string);

[result_cells query_time_LUG url_alive] = LUGQuery(query_string,xmlsettings);

ind = 1;
while isempty(result_cells{ind})
    ind = ind + 1;
end

hv = str2double(result_cells{ch+ind-1});

