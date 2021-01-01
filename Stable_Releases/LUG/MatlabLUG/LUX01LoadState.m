% function settings = LUX01LoadState(filename_prefix)
%
%   This function loads the LUG settings for the dataset filename_prefix
%
% Inputs:
%       filename_prefix - dataset name (e.g. lux01_20090328T2032)
%
% Outputs:
%       settings
%               .filename                     - same as filename_prefix
%               .source                       - data source (e.g. LED, Co-57)
%               .run_number                   - run number
%               .daq_version                  - DAQ software version (e.g. 6.1)
%               .ch                             [1x4 struct]
%                  .bias_V                    - Channel's bias voltage
%                  .serial                    - PMT serial number in channel
%                  .sphe_area_mVns            - Area of single photoelectron
%                  .gain                      - PMT Gain (over 25 Ohms)
%                  .sphe_resolution_percent   - Single pe resolution (mu/sigma)*100
%                  .PMT_calibration_filename  - Which dataset was used to
%                                               calibrate the PMT gain
%               .amplification                  [1x1 struct]
%                  .preamp                    - Preamp amplification value
%                  .postamp                     [1x1 struct]
%                     .out1_height            - Postamp OUT1 height amplification
%                     .out1_area              - Postamp OUT1 area amplification
%                     .out2_height            - Postamp OUT2 height amplification
%                     .out2_area              - Postamp OUT2 area amplification
%                     .out3_height            - Postamp OUT3 height amplification
%                     .out3_area              - Postamp OUT3 area amplification
%               .grid_HV                        [1x1 struct]
%                  .T_kV                      - Top grid voltage in kV
%                  .A_kV                      - Anode grid voltage in kV
%                  .G_kV                      - Gate grid voltage in kV
%                  .C_kV                      - Cathode grid voltage in kV
%                  .B_kV                      - Bottom grid voltage in kV
%
% 20090623 CHF
%
function settings = LUX01LoadState(filename_prefix)
%% Get LUG Settings

xmlsettings = XMLReader('defaultLUGSettings.xml');
filename = [char(39) filename_prefix char(39)];
settings.filename = filename_prefix;

%% Prepare Strings

daq_acquisition_fields      = {'filename','daq_acquisition_id','source','run_number','daq_version'...
                               'amplification_state_id','PMT_hv_state_id','grid_hv_state_id'};
PMT_hv_fields               = {'bias_ch1','bias_ch2','bias_ch3','bias_ch4'};
PMT_ch_assignment_fields    = {'ch1','ch2','ch3','ch4'};
amplification_fields        = {'preamp','postamp_out1_height','postamp_out1_area',...
                               'postamp_out2_height','postamp_out2_area','postamp_out3_height','postamp_out3_area'};
grid_hv_fields              = {'T_kV','A_kV','G_kV','C_kV','B_kV'};
PMT_gain_calibration_fields = {'filename','sphe_ch1_mVns','sphe_ch2_mVns','sphe_ch3_mVns','sphe_ch4_mVns'...
                               'sphe_resolution_ch1','sphe_resolution_ch2','sphe_resolution_ch3','sphe_resolution_ch4'};

%% DAQ Acquisition Settings

clear query_string
query_string = {'use lug'};

for ii = 1:numel(daq_acquisition_fields)
        query_string{end+1} = sprintf('select %s from daq_acquisition where filename = %s',daq_acquisition_fields{ii},filename);
end

query_string = String_Replacement_Table(query_string);

[result_cells query_time_LUG url_alive] = LUGQuery(query_string,xmlsettings);

amplification_state_id = str2double(result_cells{7});
PMT_hv_state_id = str2double(result_cells{8});
grid_hv_state_id = str2double(result_cells{9});

source = result_cells{4};
settings.source = source(1:end-2);
settings.run_number = str2double(result_cells{5});
settings.daq_version = str2double(result_cells{6});


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

for num = 1:4
    settings.ch(num).bias_V = str2double(result_cells{num+ind-1});
end


%% Get PMT Channel Assignment ID for this Run

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

[result_cells query_time_LUG url_alive] = LUGQuery(query_string,xmlsettings);

ind = 1;
while isempty(result_cells{ind})
    ind = ind + 1;
end

for num = 1:4
    serial = result_cells{num+ind-1};
    settings.ch(num).serial = serial(1:end-2);
end


%% Amplification

clear query_string result_cells
query_string = {'use lug'};

for ii = 1:numel(amplification_fields)
        query_string{end+1} = sprintf('select %s from amplification_state where amp_state_id = %d',amplification_fields{ii},amplification_state_id);
end

query_string = String_Replacement_Table(query_string);

[result_cells query_time_LUG url_alive] = LUGQuery(query_string,xmlsettings);

settings.amplification.preamp = str2double(result_cells{2});
settings.amplification.postamp.out1_height = str2double(result_cells{3});
settings.amplification.postamp.out1_area = str2double(result_cells{4});
settings.amplification.postamp.out2_height = str2double(result_cells{5});
settings.amplification.postamp.out2_area = str2double(result_cells{6});
settings.amplification.postamp.out3_height = str2double(result_cells{7});
settings.amplification.postamp.out4_area = str2double(result_cells{7});


%% Grid HV

clear query_string result_cells
query_string = {'use lug'};

for ii = 1:numel(grid_hv_fields)
        query_string{end+1} = sprintf('select %s from grid_hv_state where grid_hv_state_id = %d',grid_hv_fields{ii},grid_hv_state_id);
end

query_string = String_Replacement_Table(query_string);

[result_cells query_time_LUG url_alive] = LUGQuery(query_string,xmlsettings);

settings.grid_HV.T_kV = str2double(result_cells{2});
settings.grid_HV.A_kV = str2double(result_cells{3});
settings.grid_HV.G_kV = str2double(result_cells{4});
settings.grid_HV.C_kV = str2double(result_cells{5});
settings.grid_HV.B_kV = str2double(result_cells{6});

%% PMT Gain Calibration

for num = 1:4
    gain_calibration_id = LUX01FindBestPMTGainCalibrationDataset(num,settings.ch(num).bias_V,filename,settings.run_number);
    % gain_calibration_id = 1; %temporary

    clear query_string result_cells
    query_string = {'use lug'};

    for ii = 1:numel(PMT_gain_calibration_fields)
        query_string{end+1} = sprintf('select %s from lux01_PMT_gain_calibration where gain_calibration_id = %d',PMT_gain_calibration_fields{ii},gain_calibration_id);
    end

    query_string = String_Replacement_Table(query_string);

    [result_cells query_time_LUG url_alive] = LUGQuery(query_string,xmlsettings);

    PMT_gain_dataset = result_cells{2};

    settings.ch(num).sphe_area_mVns = str2double(result_cells{2+num});
    settings.ch(num).gain = mVns_to_gain(settings.ch(num).sphe_area_mVns);
    settings.ch(num).sphe_resolution_percent = str2double(result_cells{6+num});
    settings.ch(num).PMT_calibration_filename = PMT_gain_dataset(1:end-2);
end
