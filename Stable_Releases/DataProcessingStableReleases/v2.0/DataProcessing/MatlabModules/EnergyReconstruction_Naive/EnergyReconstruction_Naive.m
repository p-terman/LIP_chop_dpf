function status = EnergyReconstruction_Naive(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)
% 
% This is a EnergyReconstruction_Naive module.
%
% status = EnergyReconstruction_Naive(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)
%
% DESCRIPTION:
%    This module will do a very naive energy reconstruction
%
%  1) Calibrate the n_gamma in the S1 via the g1 / PDE found off the 164 keV Xenon
%  activation line (If the fields change this needs a recalibration)
%  2) Using an estimate of the 1eS2 size and the extraction efficiency
%  estimate the S2 n_e
%
% The data is stored in the following format:
%
% * energy_found_1eS2_flag (binary) - Did this dataset find a 1eS2 size out
% of the LUG (1) or did it use the default numbers (0).
% * energy_bottom_poor_fit (binary) - Was the bottom array 1eS2 considered
% a poor fit (1) or not (0) based on was it within a window of 10.5 phe
% the window is defined based on the difference between 1eS2 area (all) and
% average 1eS2 area (24.5 phe)
% * num_electrons (num_pulses,num_evts). Number of liquid electrons in all S2
% pulses (as defined by pulse_classification rq).  Outputs 0 for all non S2
% pulses.  Done over all PMTs.
% * num_electrons_bot (num_pulses,num_evts). Number of liquid electrons in all S2
% pulses (as defined by pulse_classification rq).  Outputs 0 for all non S2
% pulses.  Done over bot PMTs.
% * num_photons (num_pulses,num_evts).  Number of produced photons in all S1
% pulses (as defined by pulse_classificaiton rq.) Outputs 0 for all non S1
% pulses. Done over all PMTs.
%
% * energy_keVee_all (num_evts).  Energy reconstructed for single
% scatter events from all PMTs (keVee).
% * energy_keVee_bot (num_evts).  Energy reconstructed for single scatter
% events from just bottom PMTs (keVee).
% * energy_keVnr_all (num_evts). Energy reconstructed for single scatter
% events from all PMTs (keVnr).
% * energy_keVnr_bot (num_evts). Energy reconstructed for single scatter
% events from just the bottom PMTs (keVnr).
% * refined_energy_all (num_evts). Refined energy reconstructed for
%     single scatter events from all PMTS (keVee)
% * refined_energy_bot (num_evts). Refined energy reconstructed for
%     single scatter events from bot PMTS (keVee)
%
%   This set of outputs should allow a user to also treat multiple-scattes
%   if they so choose.
%
%
% Required RQs:
% * pulse_classification - Need to know which pulsea are S1 and S2s
% * xyz_corrected_pulse_area_all_phe
% * xyz_corrected_pulse_area_bot_phe
%
% Versioning:
%   v3.0 20140718 TPB
%      - Fixed a bug with the keVnr energy where energy_keVnr_all
%      RQ was being overwritten by energy_keVnr_bot and the energy_keVnr_bot
%      RQ was then not set.
%      - Also, new single electron calibrations don't include "adjRsqr" tag so an
%      an if statement was added to maintain compatibility with the old and new 
%      IQs
%   v2.0 20130906 PHP - Add refined energy and work towards 1eS2 size = IQ
%   v1.0 20130702 PHP - Created
%
% RQ versions:
%   v1.0 num_electrons
%   v1.0 num_photons
%   v1.0 energy_keVee_all
%   v1.0 energy_keVee_bot
%   v1.0 energy_keVnr_all
%   v1.0 energy_keVnr_bot
%
%% THIS STUFF YOU ALWAYS NEED
% Load .rq file

dp = LUXLoadRQ1s_framework(filename_rq, data_path_rq);

settings.evt_settings = dp.admin.evt_settings;
settings.daq_settings = dp.admin.daq_settings;
settings.filename_prefix = dp.admin.filename_prefix;

event_number = dp.event_number;
livetime = dp.admin.livetime;

%% Bookkeeping - access parameters

myname = 'EnergyReconstruction_Naive';
fprintf('\n\n *** Starting module %s\n',myname);

dp_settings_xml = XMLReader_framework(data_processing_xml_path);
lug_iqs_xml = XMLReader_framework(iq_xml_path);

module_names = {dp_settings_xml.data_processing_settings.module.module_name};
index_temp = strfind(module_names,myname);
index_module = find(not(cellfun('isempty', index_temp)));

mymodule_settings = dp_settings_xml.data_processing_settings.module(index_module).parameters;

max_num_pulses = dp_settings_xml.data_processing_settings.global.max_num_pulses;

%% Initialize variables

pmt_chs = 1:122;

% These should all move to iqs in version 2+
        pde = mymodule_settings.pde;
        extractioneff = mymodule_settings.extractioneff;
        NR_conversion_A = mymodule_settings.NR_conversion_A;
        NR_conversion_B = mymodule_settings.NR_conversion_B;
        W = mymodule_settings.W;

[aa b] = size(dp.event_number);
dp.energy_keVee_all = zeros(aa,b);
dp.energy_keVee_bot = zeros(aa,b);
dp.energy_keVnr_all = zeros(aa,b);
dp.energy_keVnr_bot = zeros(aa,b);

[c d] = size(dp.pulse_classification);
dp.num_electrons = zeros(c,d);
dp.num_electrons_bot = zeros(c,d);
dp.num_photons = zeros(c,d);

[dim1 dim2] = size(dp.pulse_classification);

single_scatter = repmat((sum(dp.pulse_classification==2).*sum(dp.pulse_classification==1))==1,[dim1,1]);

%% IQ fetching (here we go!)

single_e_area = [];
single_e_area_bot = [];
found_it = 0;

[a num_iqs] = size(lug_iqs_xml.iq);

 for i=1:num_iqs
  if isfield(lug_iqs_xml.iq(i),'fit')==1   

    if (isfield(lug_iqs_xml.iq(i).fit,'allpmt') && isfield(lug_iqs_xml.iq(i).fit,'botpmt'))==1        
           single_e_area = lug_iqs_xml.iq(i).fit.allpmt.mean;
           single_e_area_bot = lug_iqs_xml.iq(i).fit.botpmt.mean;
           if isfield(lug_iqs_xml.iq(i).fit.botpmt, 'adjrsquare') == 1
                   adjRsqr = lug_iqs_xml.iq(i).fit.botpmt.adjrsquare;
           else
                   adjRsqr = 1.0;
           end
           filename_prefixs = lug_iqs_xml.iq(i).global.filename_prefix;
           found_it = 1;
           dp.energy_found_1eS2_flag = true(aa,b);
           dp.energy_bottom_poor_fit = false(aa,b);
    end
  end  
 end
 
 %% Check things are kosher
 if (found_it == 0)
     % We didn't find our 1eS2 size
     dp.energy_found_1eS2_flag = false(aa,b);
     dp.energy_bottom_poor_fit = false(aa,b);
     single_e_area = mymodule_settings.single_e_area;
     single_e_area_bot = mymodule_settings.single_e_area_bot;
%  elseif ~strcmpi(filename_prefixs,settings.filename_prefix)
%      % We didn't find our 1eS2 size
%      dp.energy_found_1eS2_flag = 0;
%      dp.energy_bottom_poor_fit = 0;
%      single_e_area = mymodule_settings.single_e_area;
%      single_e_area_bot = mymodule_settings.single_e_area_bot;
%  elseif isempty(single_e_area)
%      % We didn't find our 1eS2 size
%      dp.energy_found_1eS2_flag = 0;
%      dp.energy_bottom_poor_fit = 0;
%      single_e_area = mymodule_settings.single_e_area;
%      single_e_area_bot = mymodule_settings.single_e_area_bot;
 elseif (adjRsqr < 0.97)
     single_e_area_bot = mymodule_settings.single_e_area_bot;
     dp.energy_bottom_poor_fit = true(aa,b);
 end
%% Calculate Quanta

    dp.num_photons = (dp.pulse_classification==1).*dp.xyz_corrected_pulse_area_all_phe./pde;
    dp.num_electrons = (dp.pulse_classification==2).*dp.xyz_corrected_pulse_area_all_phe./(single_e_area.*extractioneff);
    dp.num_electrons_bot = (dp.pulse_classification==2).*dp.xyz_corrected_pulse_area_bot_phe./(single_e_area_bot.*extractioneff);

%% Calculate Energy

    dp.energy_keVee_all = sum(dp.num_electrons.*single_scatter + dp.num_photons.*single_scatter,1)*W;
    dp.energy_keVee_bot = sum(dp.num_electrons_bot.*single_scatter + dp.num_photons.*single_scatter,1)*W;
    dp.energy_keVnr_all = (dp.energy_keVee_all/NR_conversion_A).^(1/NR_conversion_B);
    dp.energy_keVnr_bot = (dp.energy_keVee_bot/NR_conversion_A).^(1/NR_conversion_B);
    
%% Refine Energy (keVee)

    dp.refined_energy_all = sum(repmat(Beta(dp.energy_keVee_all),dim1,1).*dp.num_photons.*single_scatter + dp.num_electrons.*single_scatter,1).*W.*Alpha(dp.energy_keVee_all);
    dp.refined_energy_bot = sum(repmat(Beta(dp.energy_keVee_bot),dim1,1).*dp.num_photons.*single_scatter + dp.num_electrons_bot.*single_scatter,1).*W.*Alpha(dp.energy_keVee_bot);
    
%% Write Output File

status = LUXBinaryRQWriter_framework(settings, dp, filename_rq, data_path_rq, event_number, livetime); % Should add livetime input at the end, remove event_number as settings field

end

function betas = Beta(raw_energy_all)
    betas = 0.83827*ones(size(raw_energy_all)) - 0.83336 * exp(-0.46154 * raw_energy_all);
end

function alphas = Alpha(raw_energy_all)
    alphas =1.0504*ones(size(raw_energy_all)) + 1.0054*Beta(raw_energy_all) ...
        - 1.8906*(Beta(raw_energy_all)).^2 + 1.277*(Beta(raw_energy_all)).^3 ...
        - 0.41817*(Beta(raw_energy_all)).^4;
end



