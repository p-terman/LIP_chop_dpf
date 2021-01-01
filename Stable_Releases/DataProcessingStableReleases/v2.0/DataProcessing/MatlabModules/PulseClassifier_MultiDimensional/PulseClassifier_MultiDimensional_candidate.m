function status = PulseClassifier_MultiDimensional_candidate(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)
%
% status = PulseClassifier_MultiDimensional_candidate(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)
%
%
%
%
%
%
% Required RQs:
%
% peak_area_phe
% pulse_area_phe
% s2filter_max_s2_area
% s2filter_max_area_diff
% prompt_fraction
% prompt_fraction_tlx
% top_bottom_asymmetry
% skinny_pulse_area_phe
%
%
%
% Versioning:
%   20121115 CHF - Created
%   20130212 CHF - Made sure all calls are to _framework code versions
%   20130409 SAH - overhauled Carlos' PulseClassifier_s2filterspace to
%                  start including more variables and more finely-tuned
%                  cuts.  Very sloppy now, though.
%   20130415 JRV - Renamed to PulseClassifier_MultiDimensional
%                  Now skips files with 0 events
%   20130621 KOS - Updated with Lea's cuts
%   20130702 LR  - Added addtional s1 coincidence cut using skinny_area_phe
%   20130702 KOS - Adjusted S1 cuts to handle 83mKr events
%   20130821 LR  - Changed S2 threshold to 100 phe (~4 e-)
%   20130906 AL  - Moved Stable Release v1.2 to Trunk (replacing r4272)
%   20130906 AL  - Replaced prompt_fraction cut for S1 selection by
%                  prompt_fraction_tlx,copied from r4272
%   20130910 KOS - Added final 83mKr S1 allowances
%   20130911 LR  - Changed significantly - updated to a candidate version - added skinny box car cut and S1 rise time cut, removed s2filter boxcar and change some of the pre-existing cuts
%
% RQ versions:
%
%
%
%
%
%% Load .rq file

status = [];

dp = LUXLoadRQ1s_framework(filename_rq, data_path_rq);

settings.evt_settings = dp.admin.evt_settings;
settings.daq_settings = dp.admin.daq_settings;
settings.filename_prefix = dp.admin.filename_prefix;

event_number = dp.event_number;

%% Bookkeeping

myname = 'PulseClassifier_MultiDimensional';
fprintf('\n\n *** Starting module %s\n',myname);

if isempty(event_number)
    fprintf('\n\n *** Skipping module (no events in file) %s\n',myname);
    return
end

dp_settings_xml = XMLReader_framework(data_processing_xml_path);
lug_iqs_xml = XMLReader_framework(iq_xml_path);

module_names = {dp_settings_xml.data_processing_settings.module.module_name};
index_temp = strfind(module_names,myname);
index_module = find(not(cellfun('isempty', index_temp)));

mymodule_settings = dp_settings_xml.data_processing_settings.module(index_module).parameters;

max_num_pulses = dp_settings_xml.data_processing_settings.global.max_num_pulses;

%% Initialize variables

pmt_chs = 1:122;

[A B] = size(dp.pulse_area_phe);

N = length(dp.event_number);
pulse_event_size = [max_num_pulses N];

livetime = dp.admin.livetime;


%% Define Categories %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  This is where you define the cuts for S1 (cut cS1), S2 (cS2),
%  sphe (cut_sphe), SE (cut_SE) and other (cut_nota). 
%
%


% no pulse %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cut_no_pulse = dp.pulse_area_phe == 0;



% S1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% band 1: prompt fraction at 10% (tlx)
%doesnt need Kr allowance

pf_min_S1_1 = 0.7+(-1*exp(-0.27*(dp.pulse_area_phe-1)));
pf_min_S1_2 = 0.75*(1-0.8*sigmf(log10(dp.pulse_area_phe), [3 2]));
pf_energy = 18.3;

cS1pf_1 = ((dp.prompt_fraction_tlx >  pf_min_S1_1) & (dp.pulse_area_phe<pf_energy));
cS1pf_2 = ((dp.prompt_fraction_tlx >  pf_min_S1_2) & (dp.pulse_area_phe>=pf_energy));


cS1pf = cS1pf_1 | cS1pf_2;


%band 2: ratios of width

w_max_S1 = exp(-0.04*dp.pulse_area_phe)+0.3-8e-1*(1./(exp(-1*(log10(dp.pulse_area_phe)-6))));
%Kr allowance from a previous version but no change in aft_x_samples variables
w_max_S1(inrange(dp.pulse_area_phe, [150, 500])) = .34;

cS1w = (dp.aft_t1_samples-dp.aft_t0_samples)./(dp.aft_t2_samples-dp.aft_t0_samples)<w_max_S1;


%band 3: S1 riste time
%does this one need a Kr allowance? possibly not

rt_max_S1 = 13;

cS1rt_first = (dp.aft_t1_samples-dp.aft_t0_samples)<rt_max_S1;

rt_energy = 50;

cS1rt_alt = (dp.pulse_area_phe>rt_energy);

cS1rt = cS1rt_first | cS1rt_alt;


% band 4: skinny box cut
%needs Kr allowance!!!!
skinny_min_1 = 0.7+(-1*exp(-0.3*dp.pulse_area_phe));
skinny_min_2 = 0.75*(1-0.5*sigmf(log10(dp.pulse_area_phe), [3 1.9]));
skinny_energy = 19.8;


cS1skinny_1 = (dp.skinny_pulse_area_phe./dp.pulse_area_phe>skinny_min_1) & (dp.pulse_area_phe<skinny_energy);
cS1skinny_2 = (dp.skinny_pulse_area_phe./dp.pulse_area_phe>skinny_min_2) & (dp.pulse_area_phe>=skinny_energy);


cS1skinny = cS1skinny_1 | cS1skinny_2;



%pmt coincidence
    
singlephe_peak_threshold=0.1;
cut_which_tubes_hit_height = (dp.peak_height_phe_per_sample > singlephe_peak_threshold);
numtubes_height = squeeze(sum(cut_which_tubes_hit_height,2));

cut_two_tubes_height = (numtubes_height>=2);

%if there's only one event in the evt some arrays me be transposed, this
%messes up the combining the cuts
if size(cut_two_tubes_height,2) == 1
    cut_two_tubes_height=cut_two_tubes_height';
end

skinny_area_threshold=0.3;
cut_which_tubes_hit_sarea = (dp.skinny_peak_area_phe > skinny_area_threshold);
numtubes_sarea = squeeze(sum(cut_which_tubes_hit_sarea,2));

cut_two_tubes_sarea = (numtubes_sarea>=2);

%if there's only one event in the evt some arrays me be transposed, this
%messes up the combining the cuts
if size(cut_two_tubes_sarea,2) == 1
cut_two_tubes_sarea=cut_two_tubes_sarea';
end

cut_two_tubes = cut_two_tubes_height & cut_two_tubes_sarea;



% combine all three bands
cS1 = cS1w & cS1rt & cS1skinny & cS1pf;

cS1a = cS1w & cS1rt & cS1skinny & cS1pf & cut_two_tubes;



% single phe
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



cut_one_tube = (numtubes_height == 1);


%if there's only one event in the evt some arrays me be transposed, this
%messes up the combining the cuts
if size(cut_one_tube,2) == 1
    cut_one_tube=cut_one_tube';
end

%skinny box cut for SPEs

skinny_min_spe_1 = 0.9+(-4*exp(-3*dp.pulse_area_phe));
skinny_min_spe_2 = 0.7 + (1.5*log(dp.pulse_area_phe)).^-1;
skinny_energy_spe = 28;

cS1skinny_spe_1 = (dp.skinny_pulse_area_phe./dp.pulse_area_phe>skinny_min_spe_1) & (dp.pulse_area_phe<skinny_energy_spe);
cS1skinny_spe_2 = (dp.skinny_pulse_area_phe./dp.pulse_area_phe>skinny_min_spe_2) & (dp.pulse_area_phe>=skinny_energy_spe);


cS1skinny_spe = cS1skinny_spe_1 | cS1skinny_spe_2;

%collect all cuts

cut_sphe = cut_one_tube & cS1 & cS1skinny_spe;





% S2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     
% band 1: prompt fraction


pf_max_S2 = 0.7*(1-0.7*sigmf(log10(dp.pulse_area_phe), [3.5 2]));


cS2pf = (dp.prompt_fraction <  pf_max_S2) ;

    
% band 2: top-bottom asymmetry
    

S2mean = 0.03-0.50*sigmf(log10(dp.pulse_area_phe), [2 4.8]); z_min_S2 = S2mean - (0.4*log(0.56*dp.pulse_area_phe)).^-1.5; z_min_S2(z_min_S2 < -1.1) = -1.1; z_min_S2(dp.pulse_area_phe < 4) = -1.1 ;
S2mean = 0.15*(1-sigmf(log10(dp.pulse_area_phe), [2 5]));  z_max_S2 = S2mean + (0.4*log(dp.pulse_area_phe)).^-0.5;  z_max_S2(z_max_S2 > 1.1) = 1.1; z_max_S2(dp.pulse_area_phe < 1) = 1.1;



cS2z  = (dp.top_bottom_asymmetry > z_min_S2) & ...
        (dp.top_bottom_asymmetry < z_max_S2) ;


% band 3: skinny box cut

skinny_max = 0.82*(1-0.5*sigmf(log10(dp.pulse_area_phe), [3 2.3]));

cS2skinny = ((dp.skinny_pulse_area_phe./dp.pulse_area_phe)<skinny_max);

    

% in addition to SE, LR
h_min = 1;
cS2h = (dp.pulse_height_phe_per_sample > h_min);


S2_min_e = 100;
cS2e = (dp.pulse_area_phe > S2_min_e);

% combine all three bands
cS2 = cS2pf & cS2z & cS2skinny;

cS2a = cS2pf & cS2z & cS2h & cS2e & cS2skinny;




% single electron %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% that doesn't seem right to require a minimum threshold - at least then we
% to define a minimum threshold for S2s as well!!!!

Single_e_min = 5;
Single_e_max = 100;

cut_se_e = (dp.pulse_area_phe > Single_e_min) & ...
(dp.pulse_area_phe < Single_e_max) ;


% skinny box cut for SE

skinny_max_se = 0.80+(-1*exp(-0.52*dp.pulse_area_phe));

cS2skinny_se = ((dp.skinny_pulse_area_phe./dp.pulse_area_phe)<skinny_max_se);


%collect cuts

cut_se = cS2 & cut_se_e & cS2skinny_se;



            

% neither S1 nor S2 (nor no pulse) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cut_nota = ~cS1a & ~cut_sphe & ~cS2a & ~cut_se & ~cut_no_pulse;
%think about cut_nota again!! does is need cS1 a and cs2a?






%% Create output object %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assign each pulse a numerical category, using this system:

% 0 means no pulse (empty entry in matrix)
% 1 means S1
% 2 means S2
% 3 means single phe
% 4 means single electron S2
% 5 means none of the above (nota)

dp.pulse_classification           = NaN(pulse_event_size); % initialize

dp.pulse_classification(cS2a     ) = 2;   % first, assign both categories of S2
dp.pulse_classification(cut_se  ) = 4;   % 

dp.pulse_classification(cS1a     ) = 1;   % then, assign S1 (since low-energy events are in both cuts,
dp.pulse_classification(cut_sphe) = 3;   %                  are indistinguishible, and are more properly
                                         %                  called S1 than S2)

dp.pulse_classification(cut_nota) = 5;   % (order of this one shouldn't matter)

dp.pulse_classification(cut_no_pulse) = 0; %last, just in case the S1 and S2 definitions include zero-area pulses.

%% Write Output File

status = LUXBinaryRQWriter_framework(settings, dp, filename_rq, data_path_rq, dp.event_number, livetime); % Should add livetime input at the end, remove event_number as settings field

end
