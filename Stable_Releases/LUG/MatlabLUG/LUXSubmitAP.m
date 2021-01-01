%
% function [] = LUXsubmitAP(info_to_LUG, channel,lug_table,local_flag)
%
% Submits PMT afterpulsing data to the LUG.
% 
%
% Inputs: 
%
%                   info_to_LUG - structure containing the info. to be sent to the lug from PMTHealthMonitor.m
%                              .num_pC_main - cell array of main pulse sizes in picoCoulombs
%                              .APR - cell array of afterpulsing ratios, area afterpulses / area main
%                              .APratio - cel array containing vectors of afterpusing perecentage for each element
%                              .tus - cell array containing time data for each pmt (us)
%                              .ymean - cell array containing y-data for each element
%                              .range - cell array containing structs of the AP ranges for each element per pmt
%                              .filename - dataset name
%                              .PMTList - cell array of pmt serials
%                              .comment - comment on submission
%                              .PMTbias - cell array of bias for each pmt (V)
%                              .gid - group id of submission
%                              .username - username for lug submission
%                              .password - password for lug submission
%
%
%                     lug_table - submit to this table in the lug (AP or ap_software_testing)   
%                     channel - submit data from this PMT channel to LUG
%                     local_flag - 1 if acquisition was taken single-channel and stored in PMTdb
%
% JRV - 20091027

function [] = LUXSubmitAP(info_to_LUG,channel,lug_table,local_flag)

date_string = sprintf('%04d-%02d-%02d %02d:%02d:%02d',round(clock));

% format entry_date string
xmlsettings = XMLReader('insertLUGSettings.xml');
xmlsettings.password = info_to_LUG.password;

if ~isfield(info_to_LUG,'comment')
    info_to_LUG.comment = '';
end

if local_flag
    
    result_dai = '0';
    % zero means it's a local acquisition
    
    % get action_date from lug
    [result_ad query_time_ad url_alive_ad] = LUGQuery(sprintf('select action_date from lug.PMT_local_acquisition where filename = ''%s''',info_to_LUG.filename),xmlsettings);
    result_ad = regexprep(result_ad, '\n', '');

    % get action_user from lug
    [result_au query_time_au url_alive_au] = LUGQuery(sprintf('select action_user from lug.PMT_local_acquisition where filename = ''%s''',info_to_LUG.filename),xmlsettings);
    result_au = regexprep(result_au, '\n', '');
    
else
    % get daq_acquisition_id from LUG
    [result_dai query_time_dai url_alive_dai] = LUGQuery(sprintf('select daq_acquisition_id from lug.daq_acquisition where filename = ''%s''',info_to_LUG.filename),xmlsettings);
    result_dai = regexprep(result_dai, '\n', '');

    % get action_date from lug
    [result_ad query_time_ad url_alive_ad] = LUGQuery(sprintf('select action_date from lug.daq_acquisition where filename = ''%s''',info_to_LUG.filename),xmlsettings);
    result_ad = regexprep(result_ad, '\n', '');

    % get action_user from lug
    [result_au query_time_au url_alive_au] = LUGQuery(sprintf('select action_user from lug.daq_acquisition where filename = ''%s''',info_to_LUG.filename),xmlsettings);
    result_au = regexprep(result_au, '\n', '');
end


clear query;

AP_fields = '(action_date, action_user, entry_date, entry_user, entry_group_id, daq_acquisition_id, filename, comments, pmt_bias_V, pmt_serial, MP_pC, AP_H, AP_He, AP_NO, AP_Ar, AP_Xe, AP_total)';

AP_insert = {result_ad,result_au,date_string,info_to_LUG.username,num2str(info_to_LUG.gid),result_dai,info_to_LUG.filename,info_to_LUG.comment,...
    num2str(info_to_LUG.PMTbias{channel}), info_to_LUG.PMTList{channel}, num2str(info_to_LUG.num_pC_main{channel}),num2str(info_to_LUG.APratio{channel}.individual(1)),...
    num2str(info_to_LUG.APratio{channel}.individual(2)), num2str(info_to_LUG.APratio{channel}.individual(3)),...
    num2str(info_to_LUG.APratio{channel}.individual(4)), num2str(info_to_LUG.APratio{channel}.individual(5)),...
    num2str(info_to_LUG.APratio{channel}.main_total)};

AP_insert_string = '(';    % format AP_insert for mySQL input
for ii = 1:numel(AP_insert)
    if ii == 1
        AP_insert_string = strcat(AP_insert_string,'''',  AP_insert{ii},'''' );
    else
        AP_insert_string = strcat(AP_insert_string, ', ','''', AP_insert{ii},'''' );
    end
end
AP_insert_string = strcat(AP_insert_string, ')');

query{1} = sprintf('INSERT INTO lug.%s %s VALUES %s',lug_table,AP_fields,AP_insert_string{1});

% submit query to lug
[result_cells query_time_LUG url_alive] = LUGQuery(query,xmlsettings);


