function status = PulseFinder_TransparentRubiksCube_with_mod(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)
% status = PulseFinder_TransparentRubiksCube(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)
%
% Analysis is like a Rubik's cube, meaning it was all invented in Hungary in the 70's. No,
% no, meaning there are many aspects, but they are all just parts (faces) of the same
% thing. So if you really want to get the essence of the thing, you need a transparent
% Rubik's cube. Note that the analogy scales gracefuly to higher dimensions, but immediately 
% loses its tang and its poetry.
%
% For the literal-minded, this is a PulseFinder module and also a PulseTiming module. It
% even has designs on being a PulseClassifier module. 
%
% It has no outputs, other than writing a bunch of rqs to the dp struct
%
% Finally, despite the name, it should be noted that this module is only semi-transparent
% at present.
%
% 20130321 pfs - first sketch, based on existing module framework
% 20130325 pfs - checking backward and forward module compatibility, adding minimum
%                requisite rqs.
% 20130326 pfs - getting somewhere.
% 20130402 pfs - updating module algorithm to optimize S1 pulse finding
% 20130315 JRV - Now skips files with 0 events
% 20130625 JD & pfs - pulse edges *always* defined from box adges (not from hft)
% 20130709 pfs, KO, LR, FN, CG - if done searching region before two largest pulses, but
%                not done filling up max_num_pulses buffer, carry on searching for pulses
%                after the two largest (can be e trains, but also valid pulses)
% 20130731 sah - added new monte carlo truth value rqs to the
%                       `fields_to_ignore' cell array.  It would be nice to
%                       not have to edit this list everytime new rqs are
%                       added before pulse finder...  hmm...
% 20130801 JD, pfs, KO - fixed double use of 'range' variable bug that was 
%          stopping execution of region_a/region_b logic implemented in 20130709
%          commit.
% 20130811 pfs & JD - additional aft_tlx_samples and aft_trx_samples timing variables
% 20130911 JD & KO - dynamically decide which variables are pulse-level when sorting 
%          as a function of pulse start time, fixes crash when new event-level variables 
%          introduced upstream but not included in old explicit varibales to sort list 
% 20131011 JD - made tlx and trx timing fractions configurable using xml settings
% 20140204 SS - all rqs other than pulse_start_samples, pulse_end_samples and index_kept_sumpods
%          moved to PulseTiming_BasicSet and PulseQuantities_MinimumSet for the improvement of
%          module compaitility within the DPF.
% 20140425 SS - Added 20 sample buffer to end of pulses for testing (to
%          prevent missing SPE tails
% 20140723 AL - Initialzed RQs properly and with the correct type

%20180606 PAT - changed to allow for area calculation to find big pulse. This version allows for both chopping of pulses and not 
%		chopping if the pulses are not the req size
%20180830 PAT - added chopped flag to let someone know if a pulse has been chopped without having to look elsewhere
%% Load .rq file

status = [];

dp = LUXLoadRQ1s_framework(filename_rq, data_path_rq);
settings.evt_settings = dp.admin.evt_settings;
settings.daq_settings = dp.admin.daq_settings;
settings.filename_prefix = dp.admin.filename_prefix;
 
event_number = dp.event_number;
livetime = dp.admin.livetime;

%% Bookkeeping

myname = 'PulseFinder_TransparentRubiksCube';
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
max_num_pulses = int32(dp_settings_xml.data_processing_settings.global.max_num_pulses);% moded by PAT, mod removed
%max_num_pulses = int32(dp_settings_xml.data_processing_settings.global.max_num_pulses)
parts = double(max_num_pulses)%  added by PAT 151208
%parts=10
% note, this is a serendipitous re-use of noiseThre
sum_thresh = mymodule_settings.noiseThre; % sum pulse_area_phe at which the algorithm declares itself done

%% Load .cvt file or calculate it
 
filename_cvt = strrep(filename_evt,'evt','cvt');
 
if ~exist([data_path_evt filesep filename_cvt],'file')
    fprintf('*** WARNING: Did not find .cvt file. Running Summer Module\n');
    status = PulseCalibration_BaselineZen(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path);
    status = PODSummer_LUXSumPOD(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path);    
end
 
fprintf('Loading sum from %s\n',filename_cvt);
[cvt_struct settings] = LUXCVTLoader_framework(data_path_evt,strrep(filename_evt,'evt','cvt'));


%% Initialize variables

N = length(cvt_struct);
%max_num_pulses=3
pulse_event_size = [max_num_pulses N];
%optional_extra_iterations = 3; % 3 is the magic number (yes it is)

% 5 lines copied from PulseFinder_Naive.m
	dp.pulse_start_samples          = zeros(pulse_event_size, 'int32');
	dp.pulse_end_samples            = zeros(pulse_event_size, 'int32');
	dp.index_kept_sumpods           = zeros(pulse_event_size, 'uint8');
        dp.num_pulses_found             = zeros(N, 'uint32');
	settings.evt_settings           = dp.admin.evt_settings;
	settings.daq_settings           = dp.admin.daq_settings;

% initialise values to the minimum values of int32
   dp.pulse_start_samples(:)       = 999999;
   dp.pulse_end_samples(:)         = 999999;

                                                                    
    
%% Compute timing RQs
for evt = 1:N
    clear range pulse_time_samples* pulse_data_phe*;
    fprintf('.');    
    % Check if fields exist and has data
    if isfield(cvt_struct(evt),'sumpod_data_phe_per_sample') && cvt_struct(evt).empty == 0
                pulse_time_samples = cvt_struct(evt).sumpod_time_samples(cvt_struct(evt).sumpod_time_samples<1e5);
                pulse_time_samples_full = cvt_struct(evt).sumpod_time_samples(cvt_struct(evt).sumpod_time_samples<1e5);
                pulse_data_phe = cvt_struct(evt).sumpod_data_phe_per_sample(cvt_struct(evt).sumpod_time_samples<1e5);
              


		done_searching_region_a = false;
                pp = 1;
		fullBoxSamples = mymodule_settings.fullBoxSamples;
		while pp <= max_num_pulses % loop over PULSES
%			if pp>1 & flags(1)<mymodule_settings.skinnyBoxSamples % the pulse just found was determined to have a width smaller than the skinny box
%				fullBoxSamples = mymodule_settings.skinnyBoxSamples % so search again, using the skinnyBox as the fullBox
%			end
			if pp==3 & ~exist('range','var') % after finding the first two pulses, define search regions a and b
				%disp('** beginning search of region_a')
				range = pulse_time_samples<max(dp.pulse_start_samples(1,evt),dp.pulse_start_samples(2,evt)); % min or max, a minor quandry
				% divide up the pulse so you only look before the 2nd (in time) of the two largest pulses found
				pulse_time_samples_a = pulse_time_samples(range);
				pulse_data_phe_a = pulse_data_phe(range);
				pulse_time_samples_b = pulse_time_samples(~range);
				pulse_data_phe_b = pulse_data_phe(~range);

				pulse_time_samples = pulse_time_samples_a;
				pulse_data_phe = pulse_data_phe_a;
			end

			% DO NOT REMOVE FOLLOWING LINE -pfs
			if length(pulse_time_samples)<10; break; end; % need a cutoff, this one is a bit arbitrary, should check it...
			% perusePeeks looks for an S2 or S2-like pulse	
			% perusePeeks returns the actual sample times, not the indices to the times!
			[afTiming, pp_areas] = ...
				perusePeeks( pulse_time_samples ,pulse_data_phe ,[...
							mymodule_settings.preBoxSamples ...
							fullBoxSamples ...
							mymodule_settings.postBoxSamples ...
							mymodule_settings.edgeFraction ...
							mymodule_settings.txFraction ...
							mymodule_settings.skinnyBoxSamples ...
							mymodule_settings.maximumGap ...
							mymodule_settings.nLookAhead ...
							mymodule_settings.nLookBehind ...
							mymodule_settings.noiseThre ...
							] );

			% next check if the pulse should be thrown out
			% any pulse with only 3 samples are the result of fluctuations following S2's
                                                                     
			if (afTiming(2) - afTiming(1)) > 3 ...
			   || (afTiming(2) ~= afTiming(1)+1) && (pp_areas(2)/(afTiming(2) - afTiming(1))) > mymodule_settings.noiseThre
				keep_pulse = 1;
			else
				keep_pulse = 0;
			end
			if keep_pulse
				% next line needed for PulseQuantities (note it is NOT related to pods)
				dp.index_kept_sumpods(pp,evt) = 1; 
	
				% next, we need the "true" pulse edges in time. this is the region that will be 
				% integrated to give the pulse area
					dp.pulse_start_samples(pp,evt) 			= double((afTiming(1)));
					dp.pulse_end_samples(pp,evt) 			= double((afTiming(2)));

		
					
				
	
			else
				%disp('going to toss pulse...');%keyboard;
			end
			
			% zero part of waveform where the pulse was found
			% we HAVE to use the afTiming array for this info, as dp.pulse_start_samples and dp.pulse_end_samples
			% might not be filled
			zeroMe = (pulse_time_samples>=afTiming(1)) ...
					&(pulse_time_samples<=afTiming(2));
			pulse_data_phe(zeroMe) = 0;		

			% if we have not exhausted max_num_pulses, go and search the region *after* the 2nd (in time) largest pulse
			if (pp_areas(2)<sum_thresh) & ~done_searching_region_a
			%if (~keep_pulse) & ~done_searching_region_a
				done_searching_region_a = true;
				if exist('pulse_time_samples_b','var')
					pulse_time_samples = pulse_time_samples_b;
					pulse_data_phe = pulse_data_phe_b;
					%disp('** beginning search of region_b')
				else
					break; % nothing left to find
				end
			elseif ~(keep_pulse) & done_searching_region_a
				break; % nothing of consequence left to find
			end

			%if (pp==max_num_pulses) 
			%	break; % call it a day
			%end
					% make sure to increment the pulse number if appropriate
			if keep_pulse
				pp = pp + 1;
			end
		end % for pulse pp

		dp.num_pulses_found = uint32(sum(dp.index_kept_sumpods==1));% convenience RQ
		% Sort all pulse-level quantities by time (previously variables to be 
		% sorted were listed explicitly, this caused a crash if event-level 
		% varibales were added but not included in the list. Instead we now
		% only sort dp varibales with identical dimensions to pulse_event_size)                
		fields_all = fieldnames(dp);
		[Y,I]=sort(dp.pulse_start_samples(:,evt),1,'ascend');
		for f=1:length(fieldnames(dp))
			if isequal(size(dp.(fields_all{f})), pulse_event_size)
				dp.(fields_all{f})(:,evt) = dp.(fields_all{f})(I,evt);
			end 
		end % for fieldnames

          
    % Buffer for SPE tails (set in xml settings)
    pMax = dp.num_pulses_found(evt);
     if pMax > 0 
         
        for p = 1:pMax-1
            
            if (dp.pulse_end_samples(p,evt)+ mymodule_settings.extendPulse < dp.pulse_start_samples(p+1,evt))
                dp.pulse_end_samples(p,evt) = dp.pulse_end_samples(p,evt) + mymodule_settings.extendPulse;
            end
            if (dp.pulse_end_samples(p,evt)+ mymodule_settings.extendPulse >= dp.pulse_start_samples(p+1,evt))
                dp.pulse_end_samples(p,evt) = dp.pulse_start_samples(p+1,evt)-1;
            end
            
        end
         
        if (dp.pulse_end_samples(pMax,evt) + mymodule_settings.extendPulse <= max(pulse_time_samples_full))
            dp.pulse_end_samples(pMax,evt) = dp.pulse_end_samples(pMax,evt) + mymodule_settings.extendPulse;
        end
        if (dp.pulse_end_samples(pMax,evt) + mymodule_settings.extendPulse > max(pulse_time_samples_full))
            dp.pulse_end_samples(pMax,evt) = max(pulse_time_samples_full);
        end
       
     end
    
 %}       
        
    end % fi
   
end % for event evt
 
 
 
 
%% Start of PAT main mods 151208

detector_or_found = 1; % 1 = detector length imposed, 0 = detector length is the pulse length

minimium_sample_length =4000; %samples  changed from 1000 on 180618 - we were seeing problem where there was a fake s1 that was the end of a chopped pulse
				%this stopped at 0.2 charge. Assuming need to chop at 0.2 and not less than that. Need s2 pulses below 500k phd
				% this is achieved assuming equal spacing - with 8 allowed pulses in the 32000 samples, thus chop if length greater than 4000 samples. 
big_pulse_min = 1000; % 1000 phd this translates to 1/4th (although realy expect 1/2) of phd expected from 0.1 LIP - where we don't expect needing to chop
s1_length = 50; %samples
detector_sample_length = 32000; %samples

pmt_chs = 1:122;
pulse_event_size = [max_num_pulses N];
per_channel_event_size = [max_num_pulses length(pmt_chs) N];

new_pulse_start = zeros(parts, N, 'int32');
new_pulse_end = zeros(parts, N,  'int32');
new_index_kept_sumpods = zeros(parts, N, 'uint8');
dp.pulse_area_phe = zeros(pulse_event_size,'single');
dp.peak_area_phe = zeros(per_channel_event_size,'single');
dp.chopped_flag = uint8(zeros(1, size(dp.pulse_start_samples, 2))); %PAT added 180830

pulse_length = dp.pulse_end_samples - dp.pulse_start_samples;
delay_buffer = dp.admin.daq_settings.sis3301.global.delay_buffer;
% Figure out which iq has the pmt_gains - assuming only ONE iq was returned
% for each type
pmt_gains_mVns_per_phe = [];

% We still need the mV data...
event_struct = LUXEventLoader_framework(data_path_evt, filename_evt);

for ii = 1:length(lug_iqs_xml.iq)
    if isfield(lug_iqs_xml,'iq') && isfield(lug_iqs_xml.iq(ii),'global') && isfield(lug_iqs_xml.iq(ii).global,'iq_type')
        if strcmp(lug_iqs_xml.iq(ii).global.iq_type,'pmt_gains') == 1
            pmt_gains_mVns_per_phe = [lug_iqs_xml.iq(ii).fit.channel.mVns_per_phe];
            break
        end
    end
end
% Figure out which iq has the per_pod_spurious_area - assuming 1 iq was returned
% for each type
per_pod_spurious_area_phe = [];

for ii = 1:length(lug_iqs_xml.iq)
    if isfield(lug_iqs_xml,'iq') && isfield(lug_iqs_xml.iq(ii),'global') && isfield(lug_iqs_xml.iq(ii).global,'iq_type')
        if strcmp(lug_iqs_xml.iq(ii).global.iq_type,'per_pod_spurious_area') == 1
            per_pod_spurious_area_phe = [lug_iqs_xml.iq(ii).fit.channel.per_pod_spurious_area_mVns]./pmt_gains_mVns_per_phe;
            break
        end
    end
end

% actually calculate the areas of the 'naturally found' pulses
for evt = 1:N
    fprintf('.');
    
    % Check if fields exist and has data
    if isfield(cvt_struct(evt),'sumpod_data_phe_per_sample') && cvt_struct(evt).empty == 0
        loop_max = min([sum(isfinite(dp.index_kept_sumpods(:,evt))) max_num_pulses]);

        % for every pulse (NOT pod)
        for pp = 1:loop_max
            
            if dp.index_kept_sumpods(pp,evt) == 0
                continue
            end
            
            ch_map = find(~[cvt_struct(evt).ch(:).empty]);
            ch_map = intersect(ch_map,find(~[event_struct(evt).ch(:).empty])); %need the phe and the mV scales
            
            if ~isempty(ch_map)
                for ch = ch_map
                    peak_cut = inrange(cvt_struct(evt).ch(ch).pod_time_samples,dp.pulse_start_samples(pp,evt)-1,dp.pulse_end_samples(pp,evt)+1);

                    if ~isempty(peak_cut) && sum(peak_cut) > 0
                        peak_data_phe = cvt_struct(evt).ch(ch).pod_data_phe_per_sample(peak_cut);
                        dp.peak_area_phe(pp,ch,evt) = sum(peak_data_phe);
                        if(ch<=length(per_pod_spurious_area_phe)&&delay_buffer==0);
                            dp.peak_area_phe(pp,ch,evt) = dp.peak_area_phe(pp,ch,evt)-per_pod_spurious_area_phe(ch);
                        end
                    end
                end
            end % ch not empty
            
            dp.pulse_area_phe(pp,evt) = squeeze(sum(dp.peak_area_phe(pp,:,evt),2));
   
        end % for pulse pp

    end % fi sumpod
    
end % for event evt


%% now chop pulses as needed

for evt = 1 : N % note that N is the total number of events, requested earlier in pulsefinder
   
   long_pulses = find( pulse_length(:,evt) > minimium_sample_length, 1); % is there a pulse above the req length?
   
   if ~isempty(long_pulses)
	dp.chopped_flag(evt) = 1;
   	large_pulses = find(dp.pulse_area_phe(:,evt) > big_pulse_min);
	
        if isempty(large_pulses)
            first_large_pulse = 1 ; %if there is not a very large pulse then take first pulse
        else       	
            first_large_pulse = large_pulses(1);  %we want to take the first large pulse to be the start of the event
        end

        for k = 1 : parts 
                
             %this was for the 'equal break' version
             %PAT added condition to segregate the s1 180523
             if detector_or_found == 0
                 detector_sample_length = pulse_length(first_large_pulse, evt);
             end
                
             if k == 1 %imposing an s1 size
                new_pulse_start(k,evt) = dp.pulse_start_samples(first_large_pulse, evt);
                new_pulse_end(k, evt) = dp.pulse_start_samples(first_large_pulse, evt) + s1_length; 
             else
                new_pulse_start(k, evt) = dp.pulse_start_samples(first_large_pulse, evt) + s1_length + ((k-2)/(parts-1)) * (detector_sample_length - s1_length) + 1;
                new_pulse_end(k, evt) = new_pulse_start(k, evt) + 1/(parts-1) * (detector_sample_length - s1_length);
                %this should divide the detector equally, into the number of 'parts'
                %NOTE: there might be rounding errors, giving you
                %slighly higher or lower number of total samples
             end % if k==1
                
             if dp.index_kept_sumpods(first_large_pulse,evt)==1
                new_index_kept_sumpods(k, evt)=1;
             end
 
        end % for k
   else % there is no long pulse found
        new_pulse_start(:,evt) = dp.pulse_start_samples(:,evt);
        new_pulse_end(:,evt) = dp.pulse_end_samples(:,evt);
        new_index_kept_sumpods(:,evt) = dp.index_kept_sumpods(:,evt);
   end % if there is a long pulse found

end % for evt
	

 %}

%version 160226 
%{
%
new_pulse_start = zeros(parts, N, 'int32');
new_pulse_end = zeros(parts, N,  'int32');
new_index_kept_sumpods = zeros(parts, N, 'uint8');
 
pulse_length = dp.pulse_end_samples - dp.pulse_start_samples;
minimium_sample_length= 1000;
parts = pulse_length/1000;
parts = max(parts,[],1) %take longest pulse in event
parts = double(parts)

for evt = 1 : N; % note that N is the total number of events, requested earlier in pulsefinder
    j=1; %j is the running tally in the new 'expanded' pulse count
    
         for i=1 : max_num_pulses;  %i is the current original pulse count
     
        if pulse_length(i, evt)>= minimium_sample_length
            for k = 0 : parts(evt) - 1
                
                    new_pulse_start(k+1, evt) = dp.pulse_start_samples(i, evt) + ((k)/(parts(evt)))*pulse_length(i,evt); %detector_sample_length;
                    new_pulse_end(k+1, evt) =  dp.pulse_start_samples(i, evt) +  ((k+1)/(parts(evt)))*pulse_length(i,evt); %detector_sample_length;
                if dp.index_kept_sumpods(i,evt)==1
                       new_index_kept_sumpods(:, evt)=1;
                end
               
            end
            break
        end
        j=j+parts(evt);
        
       
    
    end
end


 %}

 %{

%modified to have 100 sample initial cutoff. but it will break up the pulses equally over the pulselenght

new_pulse_start = zeros(parts*max_num_pulses, N, 'int32');
new_pulse_end = zeros(parts*max_num_pulses, N,  'int32');
new_index_kept_sumpods = zeros(parts*max_num_pulses, N, 'uint8');
 
pulse_length = dp.pulse_end_samples - dp.pulse_start_samples;
detector_sample_length= 100%20000; % 20000 samples= 200ms allowed drift time, which is the whole of the 'real length' of the detector
 
for evt = 1 : N; % note that N is the total number of events, requested earlier in pulsefinder
    j=1; %j is the running tally in the new 'expanded' pulse count
    for i=1 : max_num_pulses;  %i is the current original pulse count
        
        detector_sample_length = pulse_length(i, evt); % new version modification
        if pulse_length(i, evt)>= 100 %detector_sample_length;
          % detector_sample_length = pulse_length(i, evt);
            for k = 0 : parts - 1; 
                %{
                this was for the 'equal break' version
                new_pulse_start(j + k, evt) = dp.pulse_start_samples(i, evt) + (k/parts)*(dp.pulse_end_samples(i, evt) - dp.pulse_start_samples(i, evt));
                new_pulse_end(j + k, evt) = dp.pulse_end_samples(i, evt) - (((parts-1)-k)/parts)*(dp.pulse_end_samples(i, evt) - dp.pulse_start_samples(i, evt));
                %this should divide the found start and end times equally, into the number of 'parts' 
                %}
                new_pulse_start(j + k, evt) = dp.pulse_start_samples(i, evt) +  (k/parts)*detector_sample_length;
                
                if k == parts-1 %last pulse, thus etrain/ not real length of detector only
                    new_pulse_end(j + k, evt) = dp.pulse_end_samples(i, evt);
                else
                    new_pulse_end(j + k, evt) = dp.pulse_start_samples(i, evt) + ((k+1)/parts)*detector_sample_length + - 1; %thus the end time -1 sample is the start of the next
                    
                end
                %this should divide the 'real detector length' into parts-1
                %equal parts, and the last new 'pulse' will be the
                %resulting e train and last part of 'real length'
                if dp.index_kept_sumpods(i,evt)==1
                       new_index_kept_sumpods(j+k, evt)=1;
                end
        
            end
        end
        j=j+parts;
    
    end
end
 
 %}

 
%{
%v3 discarded
 new_pulse_start = zeros(parts, N, 'int32');
new_pulse_end = zeros(parts, N,  'int32');
new_index_kept_sumpods = zeros(parts, N, 'uint8');
 
pulse_length = dp.pulse_end_samples - dp.pulse_start_samples;
minimium_sample_length= 1000
s1_length = 350;

for evt = 1 : N; % note that N is the total number of events, requested earlier in pulsefinder
    j=1; %j is the running tally in the new 'expanded' pulse count
    for i=1 : max_num_pulses;  %i is the current original pulse count
        
        if pulse_length(i, evt)>= minimium_sample_length;
            detector_sample_length = pulse_length(i,evt);
            
            for k = 0 : parts - 1; %k is the current 'part'
                
                %{
                this was for the 'equal break' version
                new_pulse_start(j + k, evt) = dp.pulse_start_samples(i, evt) + (k/parts)*(dp.pulse_end_samples(i, evt) - dp.pulse_start_samples(i, evt));
                new_pulse_end(j + k, evt) = dp.pulse_end_samples(i, evt) - (((parts-1)-k)/parts)*(dp.pulse_end_samples(i, evt) - dp.pulse_start_samples(i, evt));
                %this should divide the found start and end times equally, into the number of 'parts' 
                %}
                if k == 0 % first part we want is s1 only 
                    
                    new_pulse_start(j, evt) =  dp.pulse_start_samples(i, evt);
                    new_pulse_end(j, evt) =  dp.pulse_start_samples(i, evt) + s1_length;   
               
                elseif k == parts-1 %last pulse, thus etrain/ not real length of detector only
                    new_pulse_start(j + k, evt) = s1_length + dp.pulse_start_samples(i, evt) +  ((k-1)/(parts-1))*detector_sample_length;
                    new_pulse_end(j + k, evt) = dp.pulse_end_samples(i, evt);
                else
                    new_pulse_end(j + k, evt) = dp.pulse_start_samples(i, evt) + ((k)/(parts-1))*detector_sample_length + s1_length - 1; %thus the end time -1 sample is the start of the next
                    new_pulse_start(j + k, evt) = s1_length + dp.pulse_start_samples(i, evt) +  ((k-1)/(parts-1))*detector_sample_length;
                end
                %this should divide the 'real detector length' into parts-1
                %equal parts, and the last new 'pulse' will be the
                %resulting e train and last part of 'real length'
                if dp.index_kept_sumpods(i,evt)==1
                       new_index_kept_sumpods(j+k, evt)=1;
                end
        
            end
            break %this should mean that the treatment would only be on one event (the first one)
        end
        j=j+parts;
    
    end
end
%}

%
dp.pulse_end_samples= new_pulse_end;
dp.pulse_start_samples = new_pulse_start;
dp.index_kept_sumpods=new_index_kept_sumpods;
dp.num_pulses_found = uint32(sum(dp.index_kept_sumpods==1));% this will help visualux
dp.chopped = double(dp.chopped); %it likes to save things this way PAT added 180830
%{
for i =1:parts
    for j=1:N
        
        if dp.pulse_start_samples(i,j)==0 && dp.pulse_end_samples(i,j) ==0
            dp.pulse_start_samples(i,j) = NaN;
            dp.pulse_end_samples(i,j) = NaN;
            dp.index_kept_sumpods(i,j) = 0;
        end
    end
end
%dp.num_pulses_found= dp.num_pulses_found *parts; % do I need to have this line? Does it matter if I have originial or expanded num of pulses?
%end of main mods
%%
%} 

fprintf('Done!\n')

%% Write Output File

status = LUXBinaryRQWriter_framework(settings, dp, filename_rq, data_path_rq, event_number, livetime); % Should add livetime input at the end, remove event_number as settings field


%mex /Users/peter/matlab/library/LUX/Trunk/DataProcessing/MatlabModules/PulseFinder_TransparentRubiksCube/perusePeeks.c
%status = PulseFinder_TransparentRubiksCube(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path);
