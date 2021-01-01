function status = PulseClassifier_MultiDimensional(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)
%
% status = PulseClassifier_MultiDimensional(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)
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
% top_bottom_asymmetry
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
%   20140404 LR  - Re-incorporate top_bottom_asymmetry cut for S1s (previously replaced by width ratio cut) - but opened it up a bit - still needs to be verified
%   20140435 LR  - Loosen S1 width ratio cut & the s2filter boxcar filter for S1s, increase S1 coincidence threhsold for skinny area to 0.3 (from 0.25) and clean out obsulete cuts and commments
%   20140507 LR  - Loosen prompt_fraction_tlx cut to compenate for gain and baseline changes
%   20140515 LR  - Initialise pulse_classification RQ to zeros
%   20140520 LR  - Loosen top_bottom_asymmetry cut for S2s, prompt_Fraction_tlx and boxcar filter cut for S1s too accomodate for new gains (3.4) and baseline changes
%   20140529 AL  - Added explicit cast to double around operations with aft_t*_samples to avoid situations
%                  with 0/0 now that they are ints (using doubles this just ends up in a NaN).
%   20140429 AL  - Explicitly initialised the pulse_classification RQ as uint32
%   20140704 LR  - lower threshold for S2s to 33 phe (~1.5 x SE area)
%   20190227 PAT - removed 1 event fix. I have edited event loader to fix
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

% band 1: double-boxcar

bc_min_S1 = -0.01+(-0.5*exp(-1.2*(dp.pulse_area_phe)));

% modify to allow for 83mKr S1's
%bc_max_S1 = 0.055 + 9e-2*(log(dp.pulse_area_phe).^-1);  bc_max_S1(dp.pulse_area_phe<1) = 1;  bc_max_S1(bc_max_S1>1) = 1;
%bc_max_S1 = 0.075 + 6.5e-2*(log(dp.pulse_area_phe).^-1.4);  bc_max_S1(dp.pulse_area_phe<1) = 1;  bc_max_S1(bc_max_S1>1) = 1;
bc_max_S1 = 0.07 + 9e-2*(log(dp.pulse_area_phe).^-0.9);  bc_max_S1(dp.pulse_area_phe<1) = 1;  bc_max_S1(bc_max_S1>1) = 1;
bc_max_S1(inrange(dp.pulse_area_phe,[100,500])) = .4;


cS1bc =  (dp.s2filter_max_area_diff./dp.pulse_area_phe > bc_min_S1) & ...
         (dp.s2filter_max_area_diff./dp.pulse_area_phe < bc_max_S1) ;



% band 2: prompt fraction at 10% (tlx)

%pf_min_S1_1 = 0.7+(-1*exp(-0.27*(dp.pulse_area_phe-1)));
%pf_min_S1_1 = 0.62+(-1*exp(-0.33*(dp.pulse_area_phe-1)));
pf_min_S1_1 = 0.56+(-1.2*exp(-0.26*(dp.pulse_area_phe+0.2)));
%pf_min_S1_2 = 0.75*(1-0.8*sigmf(log10(dp.pulse_area_phe), [3 2]));
%pf_min_S1_2 = 0.68*(1-0.8*sigmf(log10(dp.pulse_area_phe), [3 2]));
pf_min_S1_2 = 0.68*(1-0.8*sigmf(log10(dp.pulse_area_phe), [2.6 2]));
%pf_energy = 18.3;
%pf_energy = 20.5;
pf_energy = 32.8;

cS1pf_1 = ((dp.prompt_fraction_tlx >  pf_min_S1_1) & (dp.pulse_area_phe<pf_energy));
cS1pf_2 = ((dp.prompt_fraction_tlx >  pf_min_S1_2) & (dp.pulse_area_phe>=pf_energy));


cS1pf = cS1pf_1 | cS1pf_2;

% band 3: top-bottom asymmetry 
z_min_S1 = -0.55 - (0.5*log(dp.pulse_area_phe)).^-0.7;  z_min_S1(z_min_S1 < -1.1) = -1.1;  z_min_S1(dp.pulse_area_phe < 1) = -1.1;
z_max_S1 = -0.35 + (0.3*log(dp.pulse_area_phe)).^-2.2;  z_max_S1(z_max_S1 > 1.1) = 1.1; z_max_S1(dp.pulse_area_phe < 1) = 1.1;


cS1z  = (dp.top_bottom_asymmetry > z_min_S1) & ...
        (dp.top_bottom_asymmetry < z_max_S1) ;


%band 4: ratios of width LR

%w_max_S1 = exp(-0.035*dp.pulse_area_phe)+0.28-8e-1*(1./(exp(-1*(log10(dp.pulse_area_phe)-6))));
w_max_S1 = exp(-0.025*dp.pulse_area_phe)+0.29-8e-1*(1./(exp(-1*(log10(dp.pulse_area_phe)-6))));
w_max_S1(inrange(dp.pulse_area_phe, [150, 500])) = .34;

cS1w = double(dp.aft_t1_samples-dp.aft_t0_samples)./double(dp.aft_t2_samples-dp.aft_t0_samples)<w_max_S1;


%additional cuts LR

singlephe_peak_threshold=0.09;
cut_which_tubes_hit_height = (dp.peak_height_phe_per_sample > singlephe_peak_threshold);
numtubes_height = squeeze(sum(cut_which_tubes_hit_height,2));

cut_two_tubes_height = (numtubes_height>=2);

%{
Removed this due to global fix PAT 190226 - this old fix was undoing it. 
%if there's only one event in the evt some arrays me be transposed, this
%messes up the combining the cuts
if size(cut_two_tubes_height,2) == 1
    cut_two_tubes_height=cut_two_tubes_height';
end
%}

skinny_area_threshold=0.3;
cut_which_tubes_hit_sarea = (dp.skinny_peak_area_phe > skinny_area_threshold);
numtubes_sarea = squeeze(sum(cut_which_tubes_hit_sarea,2));

cut_two_tubes_sarea = (numtubes_sarea>=2);
%{
Removed this due to global fix PAT 190226 - this old fix was undoing it. 
%if there's only one event in the evt some arrays me be transposed, this
%messes up the combining the cuts
if size(cut_two_tubes_sarea,2) == 1
cut_two_tubes_sarea=cut_two_tubes_sarea';
end
%}

cut_two_tubes = cut_two_tubes_height & cut_two_tubes_sarea;


% combine all four bands
cS1 = cS1bc & cS1pf & cS1w & cS1z;

cS1a = cS1bc & cS1pf & cS1w & cS1z & cut_two_tubes;



% single phe
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


cut_one_tube = (numtubes_height == 1);

%{
Removed this due to global fix PAT 190226 - this old fix was undoing it. 
%if there's only one event in the evt some arrays me be transposed, this
%messes up the combining the cuts
if size(cut_one_tube,2) == 1
    cut_one_tube=cut_one_tube';
end
%}
cut_sphe = cut_one_tube & cS1;




% S2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%band 1: double boxcar

bc_min_S2 = 0.12 - (1.35*log(dp.pulse_area_phe).^-1.5);  bc_min_S2(bc_min_S2 < 2e-4) = 2e-4;  bc_min_S2(dp.pulse_area_phe < 1) = 2e-4;
bc_max_S2 = ones(size(dp.pulse_area_phe));

cS2bc =  (dp.s2filter_max_area_diff./dp.pulse_area_phe > bc_min_S2) & ...
         (dp.s2filter_max_area_diff./dp.pulse_area_phe < bc_max_S2) ;

     
% band 2: prompt fraction


pf_min_S2 = -0.015+(-0.6*exp(-0.3*(dp.pulse_area_phe)));
pf_max_S2 = 0.7*(1-0.9*sigmf(log10(dp.pulse_area_phe), [3 2]));


cS2pf_first = (dp.prompt_fraction >  pf_min_S2) & ...
        (dp.prompt_fraction <  pf_max_S2) ;

    
% alterante cut from Jeremy, modified LR - includes oultier population of
% S2s, i.e. most of the here selected S2s are multiples scatters clustered together
%- should be re-tuned post new pulse finding algorithm
pf_max_S2_alt = 0.2;
pf_max_S2_alt_energy = 3;

cS2pf_alt = (dp.prompt_fraction > pf_min_S2) & (log10(dp.pulse_area_phe) >= pf_max_S2_alt_energy) & (dp.prompt_fraction < pf_max_S2_alt);

    
cS2pf = cS2pf_first | cS2pf_alt;    
    
% band 3: top-bottom asymmetry
    

%S2mean = 0.1-0.50*sigmf(log10(dp.pulse_area_phe), [2 4.8]); z_min_S2 = S2mean - (0.4*log(0.5*dp.pulse_area_phe)).^-1.5; z_min_S2(z_min_S2 < -1.1) = -1.1; z_min_S2(dp.pulse_area_phe < 4) = -1.1 ;
S2mean = 0.05-0.50*sigmf(log10(dp.pulse_area_phe), [2 4.8]); z_min_S2 = S2mean - (0.4*log(0.5*dp.pulse_area_phe)).^-1.5; z_min_S2(z_min_S2 < -1.1) = -1.1; z_min_S2(dp.pulse_area_phe < 4) = -1.1;
S2mean = 0.15*(1-sigmf(log10(dp.pulse_area_phe), [2 5]));  z_max_S2 = S2mean + (0.5*log(dp.pulse_area_phe)).^-0.5;  z_max_S2(z_max_S2 > 1.1) = 1.1; z_max_S2(dp.pulse_area_phe < 1) = 1.1;


cS2z  = (dp.top_bottom_asymmetry > z_min_S2) & ...
        (dp.top_bottom_asymmetry < z_max_S2) ;
    
    
%additional cuts LR

%possible width cut?
%w_min = 10;

%cs2w = (dp.gaus_fit_sigma_samples > w_min);



% in addition to SE, LR
h_min = 1;
cS2h = (dp.pulse_height_phe_per_sample > h_min);

%S2_min_e = 50;
%S2_min_e = 100;
S2_min_e = 33;

cS2e = (dp.pulse_area_phe >= S2_min_e);

% combine all three bands
cS2 = cS2bc & cS2pf & cS2z;

cS2a = cS2bc & cS2pf & cS2z & cS2h & cS2e;




% single electron %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Single_e_min = 5;
%Single_e_max = 100;
Single_e_max = 33;

cut_se =  cS2 & (dp.pulse_area_phe >= Single_e_min) & ...
                (dp.pulse_area_phe < Single_e_max) ;   


            

% neither S1 nor S2 (nor no pulse) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cut_nota = ~cS1a & ~cut_sphe & ~cS2a & ~cut_se & ~cut_no_pulse;





%% Create output object %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assign each pulse a numerical category, using this system:

% 0 means no pulse (empty entry in matrix)
% 1 means S1
% 2 means S2
% 3 means single phe
% 4 means single electron S2
% 5 means none of the above (nota)

dp.pulse_classification           = zeros(pulse_event_size, 'uint32'); % initialize

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
