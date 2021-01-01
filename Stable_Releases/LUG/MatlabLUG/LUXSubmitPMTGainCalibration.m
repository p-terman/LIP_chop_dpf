% function status = LUXSubmitPMTGainCalibration(username,password,filename,ch,sphe_mVns,sphe_resolution,comments)
% 
% This function 
%
% Inputs:
%
%
% Outputs:
%
%

function status = LUXSubmitPMTGainCalibration(username,password,data)

%% Get values from data
filename_prefix = data.filename_prefix;
ch = data.ch;
sphe_area_mVns = data.sphe_area_mVns;
sphe_resolution = data.sphe_resolution;
comments = data.comments;

%% Load Settings and Setup Fields

xmlsettings = XMLReader('insertLUGSettings.xml');
xmlsettings.password = password;

filename = sprintf('''%s''',filename_prefix);
settings.filename = filename_prefix;

users_fields                     = {'username'};

daq_acquisition_fields           = {'daq_acquisition_id','PMT_hv_state_id'};
                           
lux_PMT_hv_state_fields          = {'ch_assignment_id',sprintf('bias_ch%d',ch)};

lux_PMT_ch_assignments_fields    = {sprintf('ch%d',ch)};
                           
lux_PMT_gain_calibration_fields  = {'entry_group_id'};

%% Check that user account exist, otherwise exit with error

clear query_string
query_string = {'use lug'};

query_string{end+1} = sprintf('select %s from users',users_fields{1});

query_string = String_Replacement_Table(query_string);

[result_cells query_time_LUG url_alive] = LUGQuery(query_string,xmlsettings);

user_list = result_cells{2};

check = strfind(user_list,username);

if isempty(check)
    error('Username not found')
end


%% DAQ Acquisition Settings

clear query_string
query_string = {'use lug'};

for ii = 1:numel(daq_acquisition_fields)
        query_string{end+1} = sprintf('select %s from daq_acquisition where filename = %s',daq_acquisition_fields{ii},filename);
end

query_string = String_Replacement_Table(query_string);

[result_cells query_time_LUG url_alive] = LUGQuery(query_string,xmlsettings);

amplification_state_id = str2double(result_cells{2});
PMT_hv_state_id = str2double(result_cells{3});


%% Get Ch Assignment and Bias from PMT HV State Data

clear query_string
query_string = {'use lug'};

for ii = 1:numel(lux_PMT_hv_state_fields)
    query_string{end+1} = sprintf('select %s from lux01_PMT_hv_state where hv_state_id = %d',lux_PMT_hv_state_fields{ii},PMT_hv_state_id);
end

query_string = String_Replacement_Table(query_string);

[result_cells query_time_LUG url_alive] = LUGQuery(query_string,xmlsettings);

lux_PMT_ch_assignments_id = str2double(result_cells{2});
bias_V = str2double(result_cells{3});


%% Get serial number

clear query_string
query_string = {'use lug'};

query_string{end+1} = sprintf('select %s from lux01_PMT_ch_assignments where ch_assignment_id = %d',lux_PMT_ch_assignments_fields{1},lux_PMT_ch_assignments_id);

query_string = String_Replacement_Table(query_string);

[result_cells query_time_LUG url_alive] = LUGQuery(query_string,xmlsettings);

temp = result_cells{2};
PMT_serial = char(temp(1:6));


%% Get max entry_group_id from lux_PMT_gain_calibration

clear query_string
query_string = {'use lug'};

query_string{end+1} = sprintf('select %s from lux_PMT_gain_calibration',lux_PMT_gain_calibration_fields{1});

query_string = String_Replacement_Table(query_string);

[result_cells query_time_LUG url_alive] = LUGQuery(query_string,xmlsettings);

old_entry_group_id = max(str2num(result_cells{2}));

if isempty(old_entry_group_id)
    entry_group_id = 1;
else
    entry_group_id = old_entry_group_id + 1;
end
%increase current entry_group_id by 1. This will be the entry_group_id for this submission. If it's the first time, set it to 1

%% Send information to LUG

clear query_string
query_string = {'use lug'};

action_date = 0;

query_string{end+1} = sprintf('insert into lux_PMT_gain_calibration ...(action_date,action_user,entry_user,entry_group_id,daq_acquisition_id,filename,ch,bias_V,PMT_serial,sphe_area_mVns,sphe_resolution,comments) VALUES (''%s'',''%s'',''%s'',%d,%d,''%s'',%d,%d,''%s'',%3.2f,%.2.1f,''%s'');',action_date,username,username,entry_group_id,daq_acquisition_id,filename_prefix,ch,bias_V,PMT_serial,sphe_area_mVns,sphe_resolution,comments);

query_string = String_Replacement_Table(query_string);

[result_cells query_time_LUG url_alive] = LUGQuery(query_string,xmlsettings);




