% [PMT_serial PMT_bias_V amplification] = LUXGetLocalAcquisitionSettings(filename_prefix)
%
% This function will read off the PMT and acquisition settings for a local acquisition (single-channel PMT testing)
% 
% Inputs:
%           filename_prefix - dataset name
% Outputs:
%           PMT_serial      - serial of PMT being tested
%           PMT_bias_V      - bias voltage (V). Remember, biases are negative!
%           amplification   - total amplification value of dataset
%
% 20091029 CHF - Created
% 20091030 CHF - Added amplification output
%
function [PMT_serial PMT_bias_V amplification] = LUXGetLocalAcquisitionSettings(filename_prefix)

xmlsettings = XMLReader('defaultLUGSettings.xml');
filename = [char(39) filename_prefix char(39)];
settings.filename = filename_prefix;

PMT_local_acquisition_fields = {'pmt_serial','pmt_bias_V','amplification'};

%% DAQ Acquisition Settings

clear query_string
query_string = {'use lug'};

for ii = 1:numel(PMT_local_acquisition_fields)
        query_string{end+1} = sprintf('select %s from PMT_local_acquisition where filename = %s',PMT_local_acquisition_fields{ii},filename);
end

query_string = String_Replacement_Table(query_string);

[result_cells query_time_LUG url_alive] = LUGQuery(query_string,xmlsettings);

serial_temp = char(result_cells{2});

PMT_serial = serial_temp(1:6);
PMT_bias_V = str2num(result_cells{3});
amplification = str2num(result_cells{4});

                           