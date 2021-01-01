function status = PulseFinder_TransparentRubiksCube(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)
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
% 20131011 JD - made tlx and trx timing fractions configurable using xml settings

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
max_num_pulses = dp_settings_xml.data_processing_settings.global.max_num_pulses;

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
	dp.pulse_start_samples          = nan(pulse_event_size);
	dp.pulse_end_samples            = nan(pulse_event_size);
	dp.index_kept_sumpods           = nan(pulse_event_size);	
	settings.evt_settings           = dp.admin.evt_settings;
	settings.daq_settings           = dp.admin.daq_settings;

% rqs returned by this module
	dp.aft_box_start_samples        = zeros(pulse_event_size);
	dp.aft_t0_samples               = zeros(pulse_event_size);
	dp.aft_tlx_samples              = zeros(pulse_event_size);
	dp.aft_t1_samples               = zeros(pulse_event_size);
	dp.aft_trx_samples              = zeros(pulse_event_size);
	dp.aft_t2_samples               = zeros(pulse_event_size);
	dp.aft_box_end_samples          = zeros(pulse_event_size);

	dp.hft_t0_samples               = zeros(pulse_event_size);
	dp.hft_t10l_samples             = zeros(pulse_event_size);
	dp.hft_t50l_samples             = zeros(pulse_event_size);
	dp.hft_t1_samples               = zeros(pulse_event_size);
	dp.hft_t50r_samples             = zeros(pulse_event_size);
	dp.hft_t10r_samples             = zeros(pulse_event_size);
	dp.hft_t2_samples               = zeros(pulse_event_size);
	
	dp.pre_pulse_area_positive_phe  = zeros(pulse_event_size);
	dp.pre_pulse_area_negative_phe  = zeros(pulse_event_size);
	dp.post_pulse_area_positive_phe = zeros(pulse_event_size);
	dp.post_pulse_area_negative_phe = zeros(pulse_event_size);

	dp.full_box_samples             = zeros(pulse_event_size);
	dp.amis1_fraction               = zeros(pulse_event_size); % am i s1? es juan? no se...

	dp.n_samples_in_evt             = zeros(1,N);
	dp.full_evt_area_phe            = zeros(1,N);

	dp.skinny_pulse_start_samples	= nan(pulse_event_size);
	dp.skinny_pulse_end_samples     = nan(pulse_event_size);

        % JIMHACK: calculate this here so don't have to call any other modules as part of the LED data hack
        dp.pulse_area_phe               = zeros(pulse_event_size);
        dp.pulse_length_samples         = zeros(pulse_event_size);
	
%% Compute timing RQs
for evt = 1:N
    clear range pulse_time_samples* pulse_data_phe*;
    fprintf('.');    
    % Check if fields exist and has data
    if isfield(cvt_struct(evt),'sumpod_data_phe_per_sample') && cvt_struct(evt).empty == 0
                pulse_time_samples = cvt_struct(evt).sumpod_time_samples(cvt_struct(evt).sumpod_time_samples<1e5);
                pulse_data_phe = cvt_struct(evt).sumpod_data_phe_per_sample(cvt_struct(evt).sumpod_time_samples<1e5);

		dp.n_samples_in_evt(evt)      = length(pulse_time_samples);
		dp.full_evt_area_phe(evt)     = sum(pulse_data_phe);

		done_searching_region_a = false;
                pp = 1;
		fullBoxSamples = mymodule_settings.fullBoxSamples;
		while pp <= max_num_pulses % loop over PULSES
%			if pp>1 & flags(1)<mymodule_settings.skinnyBoxSamples % the pulse just found was determined to have a width smaller than the skinny box
%				fullBoxSamples = mymodule_settings.skinnyBoxSamples % so search again, using the skinnyBox as the fullBox
%			end
			if pp==3 & ~exist('range','var') % after finding the first two pulses, define search regions a and b
				%disp('** beginning search of region_a')
				range = pulse_time_samples<max(dp.aft_box_start_samples(:,evt)); % min or max, a minor quandry
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
			[afTiming hfTiming pp_areas flags] = ...
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

                        % JIM suggestion: changed for LED
                        % Override returned values from perusePeeks so that tricks code into just integrating 
                        % between 30 and 60 samples. N.B. to use this you must set the max number of pulses to 1 
                        % in the .xml config and with this method many of the timing and other RQ 
                        % variables will be junk. To do this properly would need to edit perusePeeks.c to 
                        % calculate proper timing variables for a fixed interval (pulse) between 30 and 60 samples.
                        afTiming = [30, 30, 30, 45, 60, 60, 60];
                        hfTiming = [30, 30, 30, 45, 60, 60, 60];
                        pp_areas = [1000.0, 1000.0, 1000.0];
                        flags = [31, 0, 0];

			% next check if the pulse should be thrown out
			% any pulse with only 3 samples are the result of fluctuations following S2's
			if (afTiming(7) - afTiming(1)) > 3 ...
			   || (pp_areas(2)/(afTiming(7) - afTiming(1))) > mymodule_settings.noiseThre
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
					dp.pulse_end_samples(pp,evt) 			= double((afTiming(7)));

                                % JIMHACK - fill in pulse area here as part of LED hack. N.B. we are not using 
                                % sum of peaks as is usual when filled in using PulseQuantities_MinimumSet 
                                found_pulse = (pulse_time_samples>=afTiming(1)) ...
                                             &(pulse_time_samples<=afTiming(7));
                                found_pulse_data_phe = pulse_data_phe(found_pulse);
                                dp.pulse_area_phe(pp,evt) = sum(found_pulse_data_phe);
                                dp.pulse_length_samples(pp,evt) = dp.pulse_end_samples(pp,evt)-dp.pulse_start_samples(pp,evt)+1 
	
				% error-trap and record the timing breakpoints
					% error-trap out of range values (can happen on very small, narrow pulses)
					for qq=2:6
						if afTiming(qq)<afTiming(1); afTiming(qq)=afTiming(1); end;
						if afTiming(qq)>afTiming(7); afTiming(qq)=afTiming(7); end;
					end
					for qq=1:7
						if hfTiming(qq)<afTiming(1); hfTiming(qq)=afTiming(1); end;
						if hfTiming(qq)>afTiming(7); hfTiming(qq)=afTiming(7); end;
					end
					
					% record area fraction timing
					dp.aft_box_start_samples(pp,evt)	    = double((afTiming(1)));
					dp.aft_t0_samples(pp,evt) 			    = double((afTiming(2)));
					dp.aft_tlx_samples(pp,evt) 			    = double((afTiming(3)));
					dp.aft_t1_samples(pp,evt) 			    = double((afTiming(4)));
					dp.aft_trx_samples(pp,evt) 			    = double((afTiming(5)));
					dp.aft_t2_samples(pp,evt) 			    = double((afTiming(6)));
					dp.aft_box_end_samples(pp,evt)   	    = double((afTiming(7)));
	 
					% record height fraction timing (legacy RQ)
					dp.hft_t0_samples(pp,evt)               = double(hfTiming(1));
					dp.hft_t10l_samples(pp,evt)             = double(hfTiming(2));
					dp.hft_t50l_samples(pp,evt)             = double(hfTiming(3));
					dp.hft_t1_samples(pp,evt)               = double(hfTiming(4));
					dp.hft_t50r_samples(pp,evt)             = double(hfTiming(5));
					dp.hft_t10r_samples(pp,evt)             = double(hfTiming(6));
					dp.hft_t2_samples(pp,evt)               = double(hfTiming(7));
					
				% aaand, we need to record the pre- and post- areas:
				% note we are NOT using the perusePeeks values 
				% (since it does not return negative and positive areas, separately)
					tmp_pre_t0 = max(pulse_time_samples(1),dp.pulse_start_samples(pp,evt)-mymodule_settings.preBoxSamples);
					tmp_pre_range = find((pulse_time_samples >= tmp_pre_t0) & (pulse_time_samples<=(dp.pulse_start_samples(pp,evt)-1)));
					tmp_pre_pulse_data_phe = cvt_struct(evt).sumpod_data_phe_per_sample(tmp_pre_range);
	
					dp.pre_pulse_area_positive_phe(pp,evt)  = sum(tmp_pre_pulse_data_phe(tmp_pre_pulse_data_phe>=0));
					dp.pre_pulse_area_negative_phe(pp,evt)  = sum(tmp_pre_pulse_data_phe(tmp_pre_pulse_data_phe<0));
					
					tmp_post_t2 = min(pulse_time_samples(end),dp.pulse_end_samples(pp,evt)+mymodule_settings.preBoxSamples);
					tmp_post_range = find((pulse_time_samples >= (dp.pulse_end_samples(pp,evt)+1)) & (pulse_time_samples<=tmp_post_t2));
					tmp_post_pulse_data_phe = cvt_struct(evt).sumpod_data_phe_per_sample(tmp_post_range);
					
					dp.post_pulse_area_positive_phe(pp,evt) = sum(tmp_post_pulse_data_phe(tmp_post_pulse_data_phe>=0));
					dp.post_pulse_area_negative_phe(pp,evt) = sum(tmp_post_pulse_data_phe(tmp_post_pulse_data_phe<0));
				
				% flags(1) indicates the actual fullBoxSamples used, if less than the requested value (default 200)
					dp.full_box_samples(pp,evt) = flags(1);
					
				% flags(2) is a measure of S1-ness: 0=>S2, 1=>S1 (ish!)
					dp.amis1_fraction(pp,evt) = flags(2);
					if isnan(dp.amis1_fraction(pp,evt)); dp.amis1_fraction(pp,evt) = 0; end; % reasonable asymptote
					if isinf(dp.amis1_fraction(pp,evt)); dp.amis1_fraction(pp,evt) = -1; end; % reasonable asymptote
					if dp.amis1_fraction(pp,evt)>10 | dp.amis1_fraction(pp,evt)<-10; dp.amis1_fraction(pp,evt) = -2; end; % error code
	
				% flags(3) currently being used to track the starting point of the max area skinny box, within the fullBox
					dp.skinny_pulse_start_samples(pp,evt)    = flags(3); % possible S1 start
					dp.skinny_pulse_end_samples(pp,evt)      = flags(3)+mymodule_settings.skinnyBoxSamples;
				
				% look at it. does it all hang together??
				if 0 %cvt_struct(evt).event_number == 103942
				   figure(32193); clf
				   %plot(cvt_struct(evt).sumpod_time_samples,cvt_struct(evt).sumpod_data_phe_per_sample,'k.-'); hold on
				   h=plot(pulse_time_samples,pulse_data_phe,'k.-'); set(h,'Color',[1 1 1]*0.7);
				   hold on;
	
				   plot(dp.pulse_start_samples(pp,evt),0,'go','markers',10,'markerFaceColor','g');
				   plot(dp.pulse_end_samples(pp,evt),0,'go','markers',10,'markerFaceColor','g');
	
				   plot(dp.aft_t0_samples(pp,evt),0,'bp','markers',10,'markerFaceColor','b');
				   plot(dp.aft_tlx_samples(pp,evt),0,'bp','markers',10,'markerFaceColor','b');
				   plot(dp.aft_t1_samples(pp,evt),0,'bp','markers',10,'markerFaceColor','b');
				   plot(dp.aft_trx_samples(pp,evt),0,'bp','markers',10,'markerFaceColor','b');
				   plot(dp.aft_t2_samples(pp,evt),0,'bp','markers',10,'markerFaceColor','b');
	
				   %xlim([dp.pulse_start_samples(pp,evt)-30  dp.pulse_end_samples(pp,evt)+30])
				   qut = (pulse_time_samples>=dp.pulse_start_samples(pp,evt)) & (pulse_time_samples<=dp.pulse_end_samples(pp,evt));
				   title(dis('evt %d pp %d (%1.2f phe)',evt,pp,sum(pulse_data_phe(qut))));
				   %dis('pre+: %1.1f phe ... post+: %1.1f phe',dp.pre_pulse_area_positive_phe(pp,evt),dp.post_pulse_area_positive_phe(pp,evt));
				   %dis('pre-: %1.1f phe ... post-: %1.1f phe',dp.pre_pulse_area_negative_phe(pp,evt),dp.post_pulse_area_negative_phe(pp,evt));
				   dis('am i s1: %3.2f',dp.amis1_fraction(pp,evt));
				   
				   keyboard
				end
			else
				%disp('going to toss pulse...');%keyboard;
			end
			
			% zero part of waveform where the pulse was found
			% we HAVE to use the afTiming array for this info, as dp.pulse_start_samples and dp.pulse_end_samples
			% might not be filled
			zeroMe = (pulse_time_samples>=afTiming(1)) ...
					&(pulse_time_samples<=afTiming(7));
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

		dp.num_pulses_found = sum(dp.index_kept_sumpods==1); % convenience RQ
		% sort everything by time, as requested by the troops
		[Y,I]=sort(dp.pulse_start_samples(:,evt),1,'ascend');
		fields_all = fieldnames(dp);
		fields_to_ignore = {'event_number' 'file_number' 'adc_ppe' 'adc_sds' 'zen_applied' 'admin' 'n_samples_in_evt' 'full_evt_area_phe' 'num_pulses_found', 'event_timestamp', 'event_timestamp_samples',...
                            'mc_key_index','mc_luxsim_random_seed','mc_event_number_luxsim','mc_prim_par_type','mc_prim_par_E_keV','mc_prim_par_pos_x_cm','mc_prim_par_pos_y_cm','mc_prim_par_pos_z_cm',...
                            'mc_prim_par_dir_x','mc_prim_par_dir_y','mc_prim_par_dir_z','mc_par_type','mc_NEST_num_gammas','mc_NEST_num_electrons','mc_E_keV','mc_x_cm','mc_y_cm','mc_z_cm','mc_missing_E_keV',...
                            'mc_num_scats','mc_msd_cm','mc_timestamp','mc_luxsim2evt_random_seed'};%'full_evt_area_thr_phe'};
        fields_to_sort = setdiff(fields_all,fields_to_ignore); % WARNING: the order of the fields is all shuffled here, thanks to setdiff

		for f=1:length(fields_to_sort)
			dp.(fields_to_sort{f})(:,evt) = dp.(fields_to_sort{f})(I,evt);
		end

    end % fi
end % for event evt

fprintf('Done!\n')

%% Write Output File

status = LUXBinaryRQWriter_framework(settings, dp, filename_rq, data_path_rq, event_number, livetime); % Should add livetime input at the end, remove event_number as settings field


%mex /Users/peter/matlab/library/LUX/Trunk/DataProcessing/MatlabModules/PulseFinder_TransparentRubiksCube/perusePeeks.c
%status = PulseFinder_TransparentRubiksCube(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path);
