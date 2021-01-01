% gain_calibration_id = LUX01FindBestPMTGainCalibrationDataset(Ch,bias_V,filename)
%
%
%
%
%
%
function gain_calibration_id = LUX01FindBestPMTGainCalibrationDataset(Ch,bias_V,filename,run_num)
%
%
%% Find Gain Calibrations Entry List
clear query_string result_cells
query_string = {'use lug'};

query_string{end+1} = sprintf('select hv_state_id from lux01_PMT_hv_state where bias_ch%d = %d',Ch,bias_V);

query_string = String_Replacement_Table(query_string);

xmlsettings = XMLReader('defaultLUGSettings.xml');
[result_cells query_time_LUG url_alive] = LUGQuery(query_string,xmlsettings);

hv_entry_list = str2num(result_cells{2:end});

%%
ii = 1;
for hv_set = hv_entry_list'
    % Filter out previous runs as well!
    clear query_string result_cells
    query_string = {'use lug'};
    query_string{end+1} = sprintf('select gain_calibration_id from lux01_PMT_gain_calibration where PMT_hv_state_id = %d and run_number = %d and ignore_calibration = 0',hv_set,run_num);
    query_string{end+1} = sprintf('select filename from lux01_PMT_gain_calibration where PMT_hv_state_id = %d and run_number = %d and ignore_calibration = 0',hv_set,run_num);
    query_string{end+1} = sprintf('select sphe_ch%d_mVns from lux01_PMT_gain_calibration where PMT_hv_state_id = %d and run_number = %d and ignore_calibration = 0',Ch,hv_set,run_num);
    query_string{end+1} = sprintf('select sphe_resolution_ch%d from lux01_PMT_gain_calibration where PMT_hv_state_id = %d and run_number = %d and ignore_calibration = 0',Ch,hv_set,run_num);
%     query_string{end+1} = sprintf('select ignore_calibration from lux01_PMT_gain_calibration where PMT_hv_state_id = %d',hv_set);
    
    query_string = String_Replacement_Table(query_string);
    [result_cells query_time_LUG url_alive] = LUGQuery(query_string,xmlsettings);
    
    gain_entry = str2num(result_cells{2});
%     ignore_calibration = str2num(result_cells{6});
    
    if ~isempty(gain_entry)
        
        N = numel(gain_entry);
        filename_candidates = result_cells{3};
        
        a = 1:18:((N+1)*18);
        for kk = 1:N
            datenumber(kk) = filename2datenum(filename_candidates(a(kk):(a(kk+1))));
            a = a + 3;
        end

        
        
        sphe_mVns = str2num(result_cells{4});
        sphe_res = str2num(result_cells{5});

        gain_entry_list(ii:(ii+N-1)) = gain_entry;
%         datenum_list(ii:(ii+N-1)) = datenumber;
%         sphe_mVns_list(ii:(ii+N-1)) = sphe_mVns;
%         sphe_res_list(ii:(ii+N-1)) = sphe_res;

        ii = ii + N;
    end
end

if ~exist('gain_entry_list','var')
    error(sprintf('\n\n********************************\nERROR: No PMT gain calibrations have been performed for Ch.%d HV = %d configuration!\nPlease run a PMT gain calibration and try again\n********************************',Ch,bias_V))
else
    datenum_dataset = filename2datenum(filename);
    [a b] = min(abs(datenumber - datenum_dataset));
    gain_calibration_id = gain_entry_list(b);
end