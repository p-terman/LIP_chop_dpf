function status = PositionReconstruction_MercuryI(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)
% PositionReconstruction module using the MercuryI algorithm.
% status = PositionReconstruction_MercuryI(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)
%
%
%s
% Required RQs:
%
% Versioning:
% 20130129 - CFPS Created
% 20130321 - CFPS New Versioning 
% 20130325 - CHF  Minor changes to conform to Module format.
% 20130403 - CFPS Pulse area included in the minimization; PMTs
% that are saturated are not used; improved performance
% 20130403 - pfs Multiple patches to obtain smooth operation in the new DP:
%            lines 47-58 (added)
%            line 65-69 (modified logic)
%            line 102-103 (error-trap for complex values in x_cm_old and y_cm_old)
% 20130404 - CFPS More corrections and a better energy callibration
% 20130503 - pfs added default value for rec_set.simulated=0
%            also cast dp.rec_energy to int8 for compatibility with LUXBinaryRQWriter_framework
% 20130511 - CFPS some modifications for the initial release
% 20130414 - JRV Now skips files with 0 events
%                Now saves the pre-calcuated table as hidden in the evt dir
% 20130703 - CFPS Version 1.2 with several modifications
% 20130717 - CFPS Version 1.2.1 Maximum likelihood method  implemented. Cut of the PMT pulses with a peak height smaller  than 1.5 mV.
% 20130808 - CFPS Version 1.2.2 Bug - In some situations the variable ENERGY_MINIMIZED was not defined
% 20130910 - CFPS Version 1.3 Small errors corrected and a new method based on photon counting introduced
% 20140507 - CFPS Version 2.0 Several improvement in the LRF curves
% 20140625 - CFPS Version 2.0.1 Small modifications to ensure the
% stability of the module
% and new photon counting reconstructions
% 20140723 - CFPS Version 2.1.3 Modification related with the
% variable types

module_version = 2.0;


%%
%% Load .rq file
%%

status = [];

dp = LUXLoadRQ1s_framework(filename_rq, data_path_rq);

%%
%% Definition of the settings
%%
 
if isfield(dp.admin, 'evt_settings')
    settings.evt_settings = dp.admin.evt_settings;%settings.evt_settings = dp.admin.settings.event_builder_settings;
end
if isfield(dp.admin, 'daq_settings')  
    settings.daq_settings = dp.admin.daq_settings;
end
settings.filename_prefix = dp.admin.filename_prefix;
event_number = dp.event_number;

%%
%% Bookkeeping
%%

myname = 'PositionReconstruction_MercuryI';
fprintf('\n\n *** Starting module %s - version %.1f\n',myname, module_version);

if isempty(event_number)
    fprintf('\n\n *** Skipping module (no events in file) %s\n',myname);
    return
end

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

if ~isfield(rec_set, 'file_lrfmat'),
    positions = findstr(data_path_evt, '/');    
    rec_set.file_lrfmat = [char(data_path_evt(1:positions(end))) '.' filename_evt(1:end-20) '_rec_fun.mat'];
end
rec_set.PMTS_Working = rec_set.PMTS_To_Use==1; rec_set.search_for_double_scatterer = 0; rec_set.compute_map = 0; rec_set.tests = 0; rec_set.only_radial_component = 0; rec_set.get_estimated_peak_areas = 1;

%rec_set.saturated_pmt_phe_limit
rec_set.saturated_pmt_phe_limit (1:61)= 100000; %[10000  12000  10000   7000   8000   7000  12000  10000   2500   9000   7000   8000  12000   7000   6000   6000   6000   6000   6000   8000  10000   8000 12000   6000   8000  12000   7000   8000   7000   6000  12000   8000   6000  12000  12000   7000   8000   7000   5000   5500  12000  12000   7000   8000   7000  11000  12000  10000  12000   6000   6000   4000   5000  12000   7000  12000  12000  12000  12000   4000  12000];
 
%% Remove the fields
if isfield(dp,'x_cm'), dp = rmfield(dp,'x_cm'); end
if isfield(dp,'y_cm'), dp = rmfield(dp,'y_cm'); end
if isfield(dp,'chi2'), dp = rmfield(dp,'chi2'); end
if isfield(dp,'s2_rec'), dp = rmfield(dp,'s2_rec'); end
if isfield(dp,'rec_dof'), dp = rmfield(dp,'rec_dof'); end
if isfield(dp,'reconstructed'), dp = rmfield(dp,'reconstructed'); end 


%% Prepare LRFs

id_prog = 'PositionReconstruction_MercuryI: '; 

for ii = 1:length(lug_iqs_xml.iq)
    if isfield(lug_iqs_xml.iq(ii).global, 'iq_type') 
        if strcmp(lug_iqs_xml.iq(ii).global.iq_type,'lrfs') & isfield(lug_iqs_xml.iq(ii).global, 'algorithm_name')
            if strcmp(lug_iqs_xml.iq(ii).global.algorithm_name,'functional') == 1 | strcmp(lug_iqs_xml.iq(ii).global.algorithm_name,rec_set.algorithm) 
                lrf_iq = lug_iqs_xml.iq(ii);
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
if rec_set.MLM_maxphe > 0
    for ii = 1:length(lug_iqs_xml.iq)
        if isfield(lug_iqs_xml.iq(ii).global, 'iq_type') 
            if strcmp(lug_iqs_xml.iq(ii).global.iq_type,'pmt_gains')
                lrf_iq.sphe = [lug_iqs_xml.iq(ii).fit.channel(1:122).sigma_mVns_per_phe]./[lug_iqs_xml.iq(ii).fit.channel(1:122).mVns_per_phe];
                lrf_iq.sphe(~inrange(lrf_iq.sphe, 0.1, 2.)) = 0.6;
            end
        end
    end
    if ~isfield(lrf_iq, 'sphe')
        fprintf(['PositionReconstruction_MercuryI. We did not found the PMT gains in the IQ settings file.\n']);
        fprintf(['The chi2 method will be used instead.\n']);  
        rec_set.MLM_maxphe = -100000;
    end
end

%% Some information printing

if rec_set.verbose > 1
    fprintf('\n %sLRF Type used in the reconstruction: %s\n', id_prog, lrf_iq.global.iq_type);
    fprintf(' PositionReconstruction_MercuryI: Algorithm version: %s\n', lrf_iq.global.algorithm_version);
    fprintf(' PositionReconstruction_MercuryI: Filenamename used to get the LRF: %s\n', lrf_iq.global.filename_prefix);
    fprintf(' PositionReconstruction_MercuryI: Computed by: %s on %d\n', lrf_iq.global.computed_by, lrf_iq.global.computed_date);
end

%% Selection of the events for the reconstruction. Previously in an
%% independent file
    
NEVTS = size(dp.peak_area_phe, 3); NPULS = size(dp.peak_area_phe, 1);

dp.reconstructed  = ones([NPULS, NEVTS],'uint8');

NumberOfTopPMTsWithSignal = squeeze(sum(dp.peak_area_phe(:, [1:60 121],:)>0, 2));
dp.reconstructed(squeeze(sum(dp.peak_area_phe,2))<rec_set.min_num_of_phe) = 0;
dp.rec_energy_flag = zeros(NPULS, NEVTS);

if rec_set.energy_minimization
    dp.rec_energy_flag = squeeze(sum(dp.peak_area_phe,2))>rec_set.energy_minimum_for_minimization;
else
    dp.rec_energy_flag = dp.reconstructed.*0;
end



%%% 
%%% Array equalization and the preliminary estimative of the position of the event the position of the event.
%%%

if ~isfield(lrf_iq, 'sigma_zero') lrf_iq.sigma_zero = 0.4; end

if isfield(dp, 'pulse_area_phe')
    [dp rec_set] = MercuryPrepareMinimization(dp, rec_set, lrf_iq);
    dp.reconstructed = double(dp.reconstructed);
    dp = PositionReconstruction_MercuryGetXY(dp, rec_set, lrf_iq);
    dp = rmfield(dp, {'x_cm_old', 'y_cm_old', 'minimization_flag', 'rec_energy_flag', 'peak_area_eq', 'pulse_area_eq'});
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

%%
%% Variable conversion 
%%

dp.chi2 = single(dp.chi2);
dp.x_cm = single(dp.x_cm);
dp.y_cm = single(dp.y_cm);
dp.s2_rec = single(dp.s2_rec);
dp.peak_area_rec = single(dp.peak_area_rec);
dp.rec_dof = uint8(dp.rec_dof);
dp.reconstructed = uint8(dp.reconstructed);

%livetime = LUXGetLivetime_framework(filename_evt,data_path_evt);
livetime=dp.admin.livetime;

status = LUXBinaryRQWriter_framework(settings, dp, filename_rq, data_path_rq, dp.event_number, livetime);
