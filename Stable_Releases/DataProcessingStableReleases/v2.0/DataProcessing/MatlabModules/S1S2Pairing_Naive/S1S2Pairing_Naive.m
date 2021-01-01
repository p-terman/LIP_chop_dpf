function status = S1S2Pairing_Naive(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)
%
%     ,-.   ,  ,-.  ,-.  ;-.                          .  .               
%    (   ` '| (   `    ) |  )     o     o             |\ |     o         
%     `-.   |  `-.    /  |-'  ,-: . ;-. . ;-. ,-:     | \| ,-: . . , ,-. 
%    .   )  | .   )  /   |    | | | |   | | | | |     |  | | | | |/  |-' 
%     `-'   '  `-'  '--' '    `-` ' '   ' ' ' `-|     '  ' `-` ' '   `-' 
%                                             `-' ---                    
% status = S1S2Pairing_Naive(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)
%
% DESCRIPTION:
%    This module will perform a simple logic to find S1-S2 pairs and drift time:
%       1) It will find all the S2s
%       2) It will find the first S1 before the first S2
%       3) It will record z_drift_samples for each of those S2s (relative
%           to the S1 found in step 2)
%
%
%
% Descriptions of the output rqs:
%
% * s1s2_pairing (num_pulses,num_evts). A value of 1 will be placed where
%   paired S1 and the S2 pulses are found.
%   For example, if the pulse classifier has:
%   pulse_classification(:,ievt) = [1 1 2 1 2]
%   , then
%   s1s2_pairing(:,ievt) = [1 0 1 0 1].
%   In the distant future, the pairing algorithm may be improved to allow
%   for multiple pairs per event. The subsequent pairings will have index
%   2, etc.  This Naive module only finds one index.
%
% * z_drift_samples (num_pulses,num_evts). Following the example above,
%   then
%   z_drift_samples(:,ievt) = [0 0 25454 0 29609]
%   where each entry is in the S2 position, and is the time between that S2 and the S1
%
%
%
%
% Required RQs:
%
% Versioning:
%   v1.0  20130418 CHF - Created
%   v1.01 20130507 BE - corrected minor bug on line 118
%   v1.1  20130730 SAH - switched from first S1 to `biggest S1 before first
%                         S2'. Swapped zeros for NaNs in the s1s2_pairing. And
%                         added friendly two helper rqs:  s1_before_s2 and golden.
%   v1.11 20130912 SAH - corrected definition of golden, so that s1 found after
%                         s2 don't disqualify it from being golden.
%   v2.0  20140704 LR  - initialise new RQs to uint32 explicitly
%                        remove golden RQ
%                        remove s1_before_s2 RQ
%                        initialise z_drift_samples to zeros
%                        pair all S2s after the first S1 with the first S1, ignore S2s before the first S1
% 20171231 PAT changed s1 to allow for class five and then have this take the largest s1 as the 'first' s1. 
%20180205 PAT changed to allow for s1 like class 5 and s2 like class 5, dependent on use of new subprogram which_class5_type.m 
%
% RQ versions:
%   v2.0 s1s2_pairing
%   v2.0 z_drift_samples
% PAT modified for LIPs
%% THIS STUFF YOU ALWAYS NEED
% Load .rq file

dp = LUXLoadRQ1s_framework(filename_rq, data_path_rq);  %standard loading

settings.evt_settings = dp.admin.evt_settings; %grabs the settings arrays and makes a new array out of that
settings.daq_settings = dp.admin.daq_settings;
settings.filename_prefix = dp.admin.filename_prefix;

event_number = dp.event_number; %event and livetime infomation collected
livetime = dp.admin.livetime;

%% Bookkeeping - access parameters

myname = 'S1S2Pairing_Naive';
fprintf('\n\n *** Starting module %s\n',myname);

dp_settings_xml = XMLReader_framework(data_processing_xml_path);  %read dp settings xml into program
lug_iqs_xml = XMLReader_framework(iq_xml_path); %read lug settings into program

module_names = {dp_settings_xml.data_processing_settings.module.module_name}; %list all possible module names
index_temp = strfind(module_names,myname); %find the location in the module_names array of the S1S2pairing_naive module (this line with next line)
index_module = find(not(cellfun('isempty', index_temp)));

mymodule_settings = dp_settings_xml.data_processing_settings.module(index_module).parameters; %any particular module settings, I believe that there are none

max_num_pulses = dp_settings_xml.data_processing_settings.global.max_num_pulses;

%% Initialize variables


N = length(dp.event_number);
pulse_event_size = [max_num_pulses N];

dp.s1s2_pairing    = zeros(pulse_event_size, 'uint32');
dp.z_drift_samples =   zeros(pulse_event_size, 'uint32');

dp = which_class5_type(dp);  %180205 PAT added 

%% Loop per event and assign S1 S2 pairs

for ii = 1:N
    
    % find the indices for S1 and S2 pulses
    s1_inds = find(dp.pulse_classification(:,ii) == 1); %180205 only s1 here | dp.pulse_classification(:,ii) == 5); %treating all class 5 as an s1 -  PAT
    if isempty(s1_inds) % no normal s1
        s1_inds = find(dp.s1_like_class5(:, ii)==1); %use class5 as s1 if there
    end
    s2_inds = find(dp.pulse_classification(:,ii) == 2 | dp.pulse_classification(:,ii) == 4 | dp.s2_like_class5(:,ii)==1); % treating all SE as s2. 180205 added s2_like_class5 as well

    
    % how many pulses of each type are there?
    num_s1 = length(s1_inds);
    num_s2 = length(s2_inds);
   
    % are there one or more pulses of each type?
    has_s1 = num_s1 > 0; 
    has_s2 = num_s2 > 0;
   
    
    
    if has_s1 & has_s2  % if an s1/other and s2/SE is extant, pair the s1 to the s2s
        
        % find the first s1
        first_s1_ind = find(dp.pulse_classification(:,ii) == 1 | dp.s1_like_class5(:,ii)==1,1); %180205 brought back and added s1 like class5
        
        %instead, first largest s1 or other pulse
        %[~, max_ind] = max(dp.pulse_area_phe(s1_inds, ii)) ; %the ~ is
        %just a dummy var, we want the index - added by PAT to take the
        %largest s1 as the 'first s1' removed 180205
 
       
        first_s1_ind = find(dp.pulse_classification(:,ii) == 1 | dp.s1_like_class5(:,ii)==1,1); %180205 brought back and added s1 like class5
             s2_after_s1 = find(dp.pulse_classification(:,ii) == 2 & (dp.aft_t0_samples(:,ii) > dp.aft_t0_samples(first_s1_ind,ii)) ...   % after the first s1  - this is using aft_t0 as the 'start' time of the pulse - thus this just means that the s1 is timed before the s2           
            | dp.pulse_classification(:,ii) == 4 & (dp.aft_t0_samples(:,ii) > dp.aft_t0_samples(first_s1_ind,ii)) ...  % now allows for SE
            | dp.s2_like_class5(:,ii)==1 &   (dp.aft_t0_samples(:,ii) > dp.aft_t0_samples(first_s1_ind,ii)) ); %180205 now includes s2 like class5 pulses (logical structure may be redundent but I preserved for consistancy)       
        %takes all s2s following 'THE s1' 
        num_s2_after_s1 = length(s2_after_s1);

        if num_s2_after_s1 > 0 %are there s2s that came AFTER the s1

            first_s2_ind = s2_after_s1(1); %this is the first s2 immediately following the s1 (even if other types of pulses inbetween)


            % we choose a single s1 to pair with all s2 pulses (here, the first s1)
            dp.s1s2_pairing(first_s1_ind,ii) = 1;
            
            % loop over all s2 after the first s1, pairing them with the s1, and calculating drift time
            for ss = 1:num_s2_after_s1
                
                % label as a member of the pair - that this s2 has a
                % corresponding s1
                dp.s1s2_pairing(s2_after_s1(ss),ii) = 1;
                
                % calculate a drift time with respect to the s1
                delta_t = dp.aft_t0_samples(s2_after_s1(ss),ii) - dp.aft_t0_samples(first_s1_ind,ii);
                dp.z_drift_samples(s2_after_s1(ss),ii) = delta_t; %save drift time info
                  
            end
            
            
        end
        
    end
    
end


%% Write Output File
dp.s1_like_class5 = double(dp.s1_like_class5);
dp.s2_like_class5=double(dp.s2_like_class5);

status = LUXBinaryRQWriter_framework(settings, dp, filename_rq, data_path_rq, event_number, livetime); % Should add livetime input at the end, remove event_number as settings field





