function dp = PositionReconstruction_Function(dp,data_processing_xml_path,iq_xml_path)
% PositionReconstruction module using the MercuryI algorithm.
% status = PositionReconstruction_Function(dp,data_processing_xml_path,iq_xml_path)
%
%
%
% Inputs:
% dp - The data structure with the following RQs (peak_area, peak_height_mV, event_number)
% The following inputs are optional:
% data_processing_xml_path - Settings for the reconstruction. 
% iq_xml_path - The LRF
% The program will load default files if these files are not provided.
%
%Outputs:
% dp - Output data structure
%
%EXAMPLE OF USAGE:    
%
% dp = PositionCorrection_Function(dp)    
% Required RQs:
%(peak_area, peak_height_mV, event_number)
%
% Versioning:
% 20130731 - CFPS Created as an independent function for the
% position reconstruction.
% 
% RQ versions:
% 1.2 
%

module_version = 1.2;

%% Load .rq file

status = [];

myname = 'PositionReconstruction_MercuryI';
myname_fun = 'PositionReconstruction_Function';

if nargin < 3
    position_reconstruction_path = which(myname_fun);
    IQs = dir([position_reconstruction_path(1:(end-numel(myname_fun)-2)) 'LUGSettings*.xml']);
    if numel(IQs)<1
        disp(sprintf(['***\nPositionReconstruction_Mercury: MAJOR ERROR\n ' ...
                      '- No IQs for this function.\n Please check ' ...
                      'if the settings (a xml file with the PMT gains ' ...
                      'and the lrfs and starting with LUGSettings) is in this folder\n %s\n***'], position_reconstruction_path(1:(end-numel(myname_fun)-2))));
        return
    else
        iq_xml_path = [position_reconstruction_path(1:(end-numel(myname_fun)-2)) IQs(end).name];
        disp(sprintf('***\nPositionReconstruction_Function: You have not indicated the IQs for this module.\n Loading the default IQ with the name %s and within the following path\n %s\n***', IQs(end).name, position_reconstruction_path(1:(end-numel(myname_fun)-2))));
    end
end

if nargin < 2
    myname_fun = 'PositionReconstruction_Function';
    position_reconstruction_path = which(myname_fun);
    sets = dir([position_reconstruction_path(1:(end-numel(myname_fun)-2)) 'MercurySettings_Default*.xml']);
    if numel(sets)<1
        disp(sprintf(['***\nPositionReconstruction_Mercury: MAJOR ERROR\n ' ...
                      '- No settings files for this function.\n Please check ' ...
                      'if the settings (a xml file with the settings needed for this function and starting with MercurySettings_Default) is in this folder\n %s\n***'], position_reconstruction_path(1:(end-numel(myname_fun)-2))));
        return
    else
        data_processing_xml_path = [position_reconstruction_path(1:(end-numel(myname_fun)-2)) sets(end).name];
        disp(sprintf('***\nPositionReconstruction_Function: You have not indicated the settings for this module.\n Loading the default settings with the name %s and within the following path\n %s\n***', sets(end).name, position_reconstruction_path(1:(end-numel(myname_fun)-2))));
    end
end

%%
%% Starting the file PositionReconstruction_MercuryI from the Line 66
%%

%% Bookkeeping
%%

dp_settings_xml = XMLReader_framework(data_processing_xml_path);
lug_iqs_xml = XMLReader_framework(iq_xml_path);

module_names = {dp_settings_xml.data_processing_settings.module.module_name};
index_temp = strfind(module_names,myname);
index_module = max(find(not(cellfun('isempty', index_temp))));

if isstruct(dp_settings_xml.data_processing_settings.module(index_module).parameters)
    rec_set = dp_settings_xml.data_processing_settings.module(index_module).parameters;
else
    rec_set.start = 1;
end

%%
%% Input checking for the structure rec_set.
%%

if ~isfield(rec_set, 'simulated'),               rec_set.simulated = 0; end
if ~isfield(rec_set, 'getapproximatedposition'), rec_set.getapproximatedposition = 'GuessTri'; end    
if ~isfield(rec_set, 'minimization_method'),     rec_set.minimization_method = 'fminsearch'; end
if ~isfield(rec_set, 'reconstruction_method'),   rec_set.reconstruction_method = 'lookuptable'; end
if ~isfield(rec_set, 'algorithm'),               rec_set.algorithm = 'mercuryi_functional'; end
if ~isfield(rec_set, 'compute_sd'),              rec_set.compute_sd = 0; end
if ~isfield(rec_set, 'cutzerophotonsdata'),      rec_set.cutzerophotonsdata = 0; end    
if ~isfield(rec_set, 'verbose'),                 rec_set.verbose = 1; end
if ~isfield(rec_set, 'maximum_distance_PMTevent'), rec_set.maximum_distance_PMTevent = 10; end
if ~isfield(rec_set, 'pulse_classifier'),        rec_set.pulse_classifier = 1;  end
if ~isfield(rec_set, 'min_top_pmts_with_signal'), rec_set.min_top_pmts_with_signal = 7;  end
if ~isfield(rec_set, 'border_minimization'),     rec_set.border_minimization = 'Multiple'; end
if ~isfield(rec_set, 'multiple_minimization_start_point'), rec_set.multiple_minimization_start_point = 25; end      
if ~isfield(rec_set, 'plot_positions'),          rec_set.plot_positions = 'No'; end   
if ~isfield(rec_set, 'cut_hits'),                rec_set.cut_hits = 1; end    
if ~isfield(rec_set, 'cut_first_pulse'),         rec_set.cut_first_pulse = 1; end    
if ~isfield(rec_set, 'chi2_MultipleMinimization_Point'),   rec_set.chi2_MultipleMinimization_Point = -1; end    
if ~isfield(rec_set, 'KeepAmplitude'),           rec_set.KeepAmplitude = 1; end    
if ~isfield(rec_set, 'select_all'),              rec_set.select_all = 1; end
if ~isfield(rec_set, 'PMTS_To_Use'),             rec_set.PMTS_To_Use= ...
        ones(1, 122); rec_set.PMTS_To_Use(32) = 0;  rec_set.PMTS_To_Use(5) = 0; rec_set.PMTS_To_Use(93) = 0; rec_set.PMTS_To_Use(41) = 0; end
if ~isfield(rec_set, 'energy_minimization'),     rec_set.energy_minimization = 1; end
if ~isfield(rec_set, 'minimization_start_radius'),  rec_set.minimization_start_radius = -1; end
if ~isfield(rec_set, 'only_radial_component'),   rec_set.only_radial_component = 0; end
if ~isfield(rec_set, 'get_estimated_peak_areas'), rec_set.get_estimated_peak_areas = 1; end
if ~isfield(rec_set, 'energy_minimum_for_minimization'), rec_set.energy_minimum_for_minimization = 80000; end
if ~isfield(rec_set, 'tests'),                  rec_set.tests = 0; end
if ~isfield(rec_set, 'remove_existing_fields'), rec_set.remove_existing_fields = 1; end
if ~isfield(rec_set, 'saturated_pmt_phe_limit'), rec_set.saturated_pmt_phe_limit = 10000; end % Above this
if ~isfield(rec_set, 'min_num_of_phe'),          rec_set.min_num_of_phe = 10; end    
if ~isfield(rec_set, 'PMT_MinNum'),             rec_set.PMT_MinNum = 21; end
if ~isfield(rec_set, 'MLM_maxphe'),              rec_set.MLM_maxphe = 1000; end
if ~isfield(rec_set, 'chib'),              rec_set.chib = 0; end

if ~isfield(rec_set, 'file_lrfmat'),
    %positions = findstr(data_path_evt, '/');    %% Comented from the original file
    rec_set.file_lrfmat = 'PositionReconstion_SaveFile_rec_fun.mat'; %% rec_set.file_lrfmat = [char(data_path_evt(1:positions(end))) '.' filename_evt(1:end-20) '_rec_fun.mat']; Line Modified
end

%% Remove the fields
if rec_set.remove_existing_fields
    if isfield(dp, 'x_cm'), dp = rmfield(dp, 'x_cm'); end
    if isfield(dp, 'y_cm'), dp = rmfield(dp, 'y_cm'); end
    if isfield(dp, 'chi2'), dp = rmfield(dp, 'chi2'); end
    if isfield(dp, 's2_rec'), dp = rmfield(dp, 's2_rec'); end
    if isfield(dp, 'rec_dof'), dp = rmfield(dp, 'rec_dof'); end
end

%% Prepare LRFs

id_prog = 'PositionReconstruction_MercuryI: '; 

for ii = 1:length(lug_iqs_xml.iq)
    if isfield(lug_iqs_xml.iq(ii).global, 'iq_type') 
        if strcmp(lug_iqs_xml.iq(ii).global.iq_type,'lrfs') & isfield(lug_iqs_xml.iq(ii).global, 'algorithm_name')
            if strcmp(lug_iqs_xml.iq(ii).global.algorithm_name,'functional') == 1 | strcmp(lug_iqs_xml.iq(ii).global.algorithm_name,rec_set.algorithm) 
                if lug_iqs_xml.iq(ii).global.algorithm_version <= module_version
                    lrf_iq = lug_iqs_xml.iq(ii);
                else
                    fprintf([' PositionReconstruction_MercuryI: I found the version %.1f but I need the version %.1f\n'], lug_iqs_xml.iq(ii).global.algorithm_version, module_version);
                end
            end
        end
    end
end

if ~exist('lrf_iq')
    fprintf('The LRF IQ was not found or it was found the wrong version. We need the LRF with the algorithm %s and version %.1f\n', rec_set.algorithm, module_version);
    fprintf(['Bye Bye. Please contact Claudio Silva if you cannot find the correct version\n']);
    return
end
    

%% Get PMT Gains
for ii = 1:length(lug_iqs_xml.iq)
    if isfield(lug_iqs_xml.iq(ii).global, 'iq_type') 
        if strcmp(lug_iqs_xml.iq(ii).global.iq_type,'pmt_gains')
            lrf_iq.sphe = [lug_iqs_xml.iq(ii).fit.channel(1:122).sigma_mVns_per_phe]./[lug_iqs_xml.iq(ii).fit.channel(1:122).mVns_per_phe];
            lrf_iq.sphe(~inrange(lrf_iq.sphe, 0.1, 2.)) = 0.6;
            if sum(~inrange(lrf_iq.sphe, 0.1, pi))>0
                fprintf('\n\n *** Some sphe curves width have strange values -- please check that %s\n',myname);            
            end
        end
        %% This removes non-sensical values
    end
end
if ~isfield(lrf_iq, 'sphe')
    fprintf(['PositionReconstruction_MercuryI. We did not found the PMT gains in the IQ settings file.\n']);
    fprintf(['The chi2 method will be used instead.\n']);  
    rec_set.MLM_maxphe = -100000;
end


%% Some information printing

if rec_set.verbose > 1
    fprintf('\n %sLRF Type used in the reconstruction: %s\n', id_prog, lrf_iq.global.iq_type);
    fprintf(' PositionReconstruction_MercuryI: Algorithm version: %s\n', lrf_iq.global.algorithm_version);
    fprintf(' PositionReconstruction_MercuryI: Filenamename used to get the LRF: %s\n', lrf_iq.global.filename_prefix);
    fprintf(' PositionReconstruction_MercuryI: Computed by: %s on %d\n', lrf_iq.global.computed_by, lrf_iq.global.computed_date);
end

%% Selection of the events for the reconstruction. 

[dp rec_set] = MercurySelectEventsForReconstruction(dp, rec_set);

%%% 
%%% Array equalization and the preliminary estimative of the position of the event the position of the event.
%%%

if isfield(dp, 'pulse_area_phe')
    [dp rec_set] = MercuryPrepareMinimization(dp, rec_set, lrf_iq);
    dp.reconstructed = double(dp.reconstructed);
    dp = PositionReconstruction_MercuryGetXY(dp, rec_set, lrf_iq);
    dp = rmfield(dp, {'x_cm_old', 'y_cm_old', 'minimization_flag', 'rec_energy_flag','pulse_area_eq'});
else
    disp(fprintf('\n **** %s 0 Events in the rq file %s. **** \n',myname, filename_rq));
    dp.x_cm = 100*ones(size(dp.pulse_area_phe)); % a very out of range value
    dp.y_cm = 100*ones(size(dp.pulse_area_phe)); % a very out of range value
    dp.minimization_flag = -1*ones(size(dp.pulse_area_phe)); % a failure value
    dp.chi2 = -1*ones(size(dp.pulse_area_phe)); % a failure value
end
    

if rec_set.get_estimated_peak_areas==0
    dp = rmfield(dp, 'peak_area_rec');
end
