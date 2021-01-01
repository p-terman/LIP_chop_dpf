function status = select_reprocessing_dpf_brown(LUXcode_path, data_path_evt, data_path_rq, data_processing_xml_path, iq_xml_path, ArrayID)
tic
%Inputs:
%   data_path_evt = evt files location (unzipped) for specific dataset
%   data_path_rq =  location of the rq.mat file should have alredy been processed with 100 pulses
%   LUXcode_path = path to the dpf modules location with the special modules
%   data_processing_xml_path = this is actually the dp xml settings file with location (full path) 
%   iq_xml_path = the iq xml file with location
%
% output
%    New folder in the matfiles folder called 'chopped events' which will contain the rq files
%    of just the chopped events, no combined files (to save space). Named XXXX_chop.rq
%    Subfolder of 'chopped events' called 'matfiles' for the rq.mat files
%
%  intermediate files:
%       program will make a _chop.evt.mat file and a _chop.cvt file but
%       these are deleted at the end of the program loop
%
% Note: this program will loop over all rq files in a dataset.


addpath(genpath(LUXcode_path));

lug_iqs_xml_file=iq_xml_path;

data_processing_settings_path=data_processing_xml_path;
pathbase = extractBefore(data_path_evt, 'lux'); %be careful here - this might need to be changed

if ~contains(data_path_rq, 'matfiles')
    data_path_rq = strcat(data_path_rq, '/matfiles/');
end

filelist = dir([data_path_rq filesep 'lux*.rq.mat']);
filelist_length = length(filelist);
ArrayID = str2num(ArrayID);

if ArrayID * 100 < filelist_length
	end_point = ArrayID * 100;
	start_point = end_point - 99;
	disp('choice one')
else
	end_point = filelist_length;
	start_point = (ArrayID - 1) * 100 + 1;
	disp('choice else') 
end
	
disp(ArrayID)

disp(length(filelist))
for ii= start_point : end_point
   
    filename_rq = filelist(ii).name;
    disp(filename_rq)	
    namebase = filename_rq(1:27);

    filename_evt = dir([data_path_evt filesep namebase '*.evt']); % previous script should have made sure that this is not zip
    disp(filename_evt)
    filename_evt = filename_evt.name;
    
    dp = load([data_path_rq filesep filename_rq ]);
    
    if contains(filename_rq , '.mat')
        filename_rq = erase(filename_rq, '.mat');
    end
    settings.evt_settings = dp.admin.evt_settings;
    settings.daq_settings = dp.admin.daq_settings;
    settings.filename_prefix = dp.admin.filename_prefix;


    livetime = dp.admin.livetime;

    %% Bookkeeping

    dp_settings_xml = XMLReader_framework(data_processing_xml_path);
    lug_iqs_xml = XMLReader_framework(iq_xml_path);
    data_processing_settings_path = data_processing_xml_path;
    lug_iqs_xml_file = iq_xml_path;

    %%


    if ~isfield(dp, 'pulse_end_samples')
        fprintf('Should have run Pulse Height Timing and Pulse Quantities Minimum Set !\n')
        system([strcat(LUXcode_path, '/CppModules/bin/PulseTiming_HeightTiming ') filename_evt ' ' data_path_evt ' ' filename_rq ' ' pathbase ' 5 ' data_processing_settings_path ' ' lug_iqs_xml_file ])
        status = PulseQuantities_MinimumSet(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path);
        dp = LUXLoadRQ1s_framework(filename_rq, data_path_rq);
    end


    %set your constants
    s1_length = 50; %samples
    minimum_sample_length = 4000; %samples required to chop a pulse
    max_num_pulses = dp_settings_xml.data_processing_settings.global.max_num_pulses;
    parts = double(max_num_pulses); % number of parts that we will chop into
    detector_sample_length= 32000;% 320ms allowed drift time, which is the whole of the 'real length' of the detector
    min_pulse_area = 1000; %phd this is the minimum area required to say that this is the 'first major pulse' in event
                         %this min area will help to screen noise but keep s1 -
                         %remember, looking for large pulses, normal dpf will
                         %get the small s1s


    is_long_pulse = dp.pulse_length_samples > minimum_sample_length;
    re_process = find(sum(is_long_pulse));

    if ~isempty(re_process) % an event has a pulse that is long

        event_struct = LUXEventLoader_framework(data_path_evt, filename_evt);

        chop_filename_rq = strrep(filename_rq , '.rq' , '_chop.rq');

        %% chop down the evt file 
        evt_list = fieldnames(event_struct);
        clear new_evt_file       
        N_new = length(re_process);
	evt_list_length = length(evt_list);
        for i = 1:evt_list_length
            current_evt = evt_list{i};

            for j = 1:N_new % go through each event
                event = re_process(j);

                if j==1
                    if strcmp(current_evt, 'ch_map') || strcmp(current_evt, 'filename_prefix') ||...
                        strcmp(current_evt, 'posttrigger') || strcmp(current_evt, 'pretrigger')

                        new_evt_file(j).(current_evt) = event_struct(:,1).(current_evt);

                    elseif strcmp(current_evt, 'filename')
                        filename_evt_new = [filename_evt(1:end-4) '_chop.evt'];
                        new_evt_file(j).(current_evt) = filename_evt_new;
                    else
                        new_evt_file(j).(current_evt) = event_struct(:, event).(current_evt); 
                    end
                elseif j == N_new
                    if  strcmp(current_evt, 'filename_prefix') ||...
                        strcmp(current_evt, 'posttrigger') || strcmp(current_evt, 'pretrigger')

                        new_evt_file(j).(current_evt) = event_struct(:,end).(current_evt);
                    else
                        new_evt_file(j).(current_evt) = event_struct(:,event).(current_evt); 
                    end

                else
                    % j not 1 or last 
                new_evt_file(j).(current_evt) = event_struct(:, event).(current_evt); 
                end
            end
        end % evt_list

        %% now we want to run another pass at the evt info to see what else we
        % can eliminate. First generate a cvt like object
        %this is like the baseline zen
        amp_gain = dp.admin.daq_settings.global.preamp .* dp.admin.daq_settings.global.postamp;

        % Figure out which iq has the pmt_gains - assuming only ONE iq was returned
        % for each type
        pmt_gains_mVns_per_phe = [];
	iqs_length = length(lug_iqs_xml.iq);
        for rr = 1:iqs_length
            if isfield(lug_iqs_xml,'iq') && isfield(lug_iqs_xml.iq(rr),'global') && isfield(lug_iqs_xml.iq(rr).global,'iq_type')
                if strcmp(lug_iqs_xml.iq(rr).global.iq_type,'pmt_gains') == 1
                    pmt_gains_mVns_per_phe = [lug_iqs_xml.iq(rr).fit.channel.mVns_per_phe];
                 break
                end
            end
        end

        myname = 'PulseCalibration_BaselineZen';

        module_names = {dp_settings_xml.data_processing_settings.module.module_name};
        index_temp = strfind(module_names,myname);
        index_module = find(not(cellfun('isempty', index_temp)));

        if ~isempty(index_module)
            dp_settings_xml.data_processing_settings.module(index_module);
            mymodule_settings = dp_settings_xml.data_processing_settings.module(index_module).parameters;
        else
        %     error(sprintf('*** ERROR: Module was not found in settings file:\n%s\n',data_processing_xml_path));
        end

        holder = new_evt_file;% do this to the reduced evt struct.
        holder(1).thr = mymodule_settings.flatten_thr; % simple way to port this value into the function
        holder =  LUXCalibratePulses_framework(holder,pmt_gains_mVns_per_phe,amp_gain);
        cvt_struct = LUXSumPOD_framework(holder); %this is instead of the pod summing module
        clear holder;

        final_ind = 1;
        clear final_evt_new new_cvt_file; %make sure on a clean slate

	re_process_length = length(re_process)
        for p = 1:re_process_length % still need these numbers because we are also looking at the rq file
            event = re_process(p);
            long_pulses = find(dp.pulse_length_samples(:,event)>4000);
            first_long_pulse = long_pulses(1);
            pulse_start = dp.pulse_start_samples(first_long_pulse, event);
            pulse_end = dp.pulse_end_samples(first_long_pulse, event);

            if dp.pulse_length_samples(first_long_pulse, event) > 32000
                pps_test_mid = pulse_start + 16000;
                pps_test_twothirds = pulse_start + 21000;
              %  pps_test_threequarters = pulse_start + 24000;
                pps_test_end = pulse_start + 31500;
            else
                pps_test_mid = pulse_start + (pulse_end - pulse_start)/2;
                pps_test_twothirds = pulse_start + (pulse_end - pulse_start)*2/3;
             %   pps_test_threequarters = pulse_start + (pulse_end - pulse_start)*3/4;
                pps_test_end = pulse_end - 500; % back off a bit from the exact end since this will be going toward zero
            end


            [~ , mid_location_start ] = min( abs(cvt_struct(p).sumpod_time_samples - double(pps_test_mid - 100)  ) );
            [~ , mid_location_end ] = min( abs(cvt_struct(p).sumpod_time_samples - double(pps_test_mid + 100)  ) );

            [~ , range_max_loc] = max(cvt_struct(p).sumpod_data_phe_per_sample(mid_location_start : mid_location_end));
            mid_max_loc = mid_location_start + range_max_loc; % this should be the location of the max in terms of cvt samples                                      

            [~ , twothirds_location_start ] = min( abs(cvt_struct(p).sumpod_time_samples - double(pps_test_twothirds - 100)  ) );
            [~ , twothirds_location_end ] = min( abs(cvt_struct(p).sumpod_time_samples - double(pps_test_twothirds + 100)  ) );

            [~ , range_max_loc] = max(cvt_struct(p).sumpod_data_phe_per_sample(twothirds_location_start : twothirds_location_end));
            twothirds_max_loc = twothirds_location_start + range_max_loc;

            %[~ , threequarters_location_start ] = min( abs(cvt_struct(p).sumpod_time_samples - double(pps_test_threequarters - 100)  ) );
            %[~ , threequarters_location_end ] = min( abs(cvt_struct(p).sumpod_time_samples - double(pps_test_threequarters + 100)  ) );

            %[~ , range_max_loc] = max(cvt_struct(p).sumpod_data_phe_per_sample(threequarters_location_start : threequarters_location_end));
            %threequarters_max_loc = threequarters_location_start + range_max_loc;

            [~ , end_location_start ] = min( abs(cvt_struct(p).sumpod_time_samples - double(pps_test_end - 100)  ) );
            [~ , end_location_end ] = min( abs(cvt_struct(p).sumpod_time_samples - double(pps_test_end + 100)  ) );

            [~ , range_max_loc] = max(cvt_struct(p).sumpod_data_phe_per_sample(end_location_start : end_location_end));
            end_max_loc = end_location_start + range_max_loc;

            pps_mid_ave = sum(cvt_struct(p).sumpod_data_phe_per_sample(mid_max_loc - 25 : mid_max_loc + 25))/50;

            pps_twothirds_ave = sum(cvt_struct(p).sumpod_data_phe_per_sample(twothirds_max_loc - 25 : twothirds_max_loc + 25))/50; 

            %pps_threequarters_ave = sum(cvt_struct(p).sumpod_data_phe_per_sample(threequarters_max_loc - 25 : threequarters_max_loc + 25))/50; 

            pps_end_ave = sum(cvt_struct(p).sumpod_data_phe_per_sample(end_max_loc - 25 : end_max_loc + 25)) / 50;

            pps_threshold = 4; %phd per sample required to keep the pulse

            %threshold tests
            tm = pps_mid_ave > pps_threshold;
            t23 = pps_twothirds_ave > pps_threshold;
            %t34 = pps_threequarters_ave > pps_threshold;
            te = pps_end_ave > pps_threshold;

            if (tm && t23) ||  (tm && te)||...(tm && t34) ||
                     (t23 && te) %|| (t23 && t34) ||(t34 && te)
                     % if two out of four are above threshold then we keep it
                final_evt_new(final_ind) = new_evt_file(p);
                new_cvt_file(final_ind) = cvt_struct(p);
                final_re_process(final_ind) = re_process(p);
                final_ind = final_ind + 1;
            end       
        end
        clear new_evt_file; %don't need this, we have final_evt_new and event_struct

        if final_ind == 1 % there has been no entry to the final structs.
            continue;
        end

        evt_list = fieldnames(event_struct);          
	evt_list_length = length(evt_list);
        for j = 1: evt_list_length
            current_evt = evt_list{j};

            if strcmp(current_evt, 'ch_map')
                final_evt_new(1).(current_evt) = event_struct(:,1).(current_evt);
            elseif strcmp(current_evt, 'filename')
                filename_evt_new = [filename_evt(1:end-4) '_chop.evt'];
                final_evt_new(1).(current_evt) = filename_evt_new;
            elseif  strcmp(current_evt, 'filename_prefix') ||...
                    strcmp(current_evt, 'posttrigger') || strcmp(current_evt, 'pretrigger')

                final_evt_new(1).(current_evt) = event_struct(:,1).(current_evt);
                final_evt_new(end).(current_evt) = event_struct(:,end).(current_evt);
            end
        end

        N_new = length(final_evt_new);

        %save the new file
        fprintf('Saving EVT file... ');
        filename_evt_fullpath_new = [data_path_evt filesep filename_evt_new '.mat'];
        save(filename_evt_fullpath_new, 'final_evt_new' , '-v7.3');
        chop_filename_evt = filename_evt_new;
        clear final_evt_new event_struct; % we are done with these, save the RAM
      %  clear new_cvt_file; % make sure we are on a clean slate
            
        cvt_list = fieldnames(cvt_struct);
	cvt_list_length = length(cvt_list);
        for jj = 1: cvt_list_length
            current_cvt = cvt_list{jj};

            if strcmp(current_cvt, 'ch_map')
                new_cvt_file(1).(current_cvt) = cvt_struct(:,1).(current_cvt);

            elseif  contains(current_cvt, 'filename') || strcmp(current_cvt, 'thr') ...
                   || strcmp(current_cvt, 'posttrigger') || strcmp(current_cvt, 'pretrigger')
 
                new_cvt_file(end).(current_cvt) = cvt_struct(:,end).(current_cvt);
                new_cvt_file(1).(current_cvt) = cvt_struct(:,1).(current_cvt);

            end
        end


        %save the new file
        fprintf('Saving CVT file... ');
        chop_filename_cvt = strrep(filename_evt , '.evt' , '_chop.cvt');
        filename_cvt_fullpath_new = [data_path_evt filesep chop_filename_cvt];

        status = LUXCVTWriter_framework( new_cvt_file, settings, livetime, filename_cvt_fullpath_new );
        clear new_cvt_file cvt_struct dp_chop;

        dp_chop.admin = dp.admin;
        dp_chop.file_number =  dp.file_number(final_re_process);
        dp_chop.event_number = dp.event_number(final_re_process);
        dp_chop.source_filename = dp.source_filename;
        dp_chop.event_timestamp_samples = dp.event_timestamp_samples(final_re_process);
        dp_chop.luxstamp_samples = dp.luxstamp_samples(final_re_process);
        dp_chop.time_since_livetime_start_samples = dp.time_since_livetime_start_samples(final_re_process);
        dp_chop.time_until_livetime_end_samples = dp.time_until_livetime_end_samples(final_re_process);

        %% Section in place of Transparent Rubix Cube
        % look at the finalized non-chop result
        % then find where to chop

        new_pulse_start = zeros(parts, N_new, 'int32');
        new_pulse_end = zeros(parts, N_new,  'int32');
        new_index_kept_sumpods = zeros(parts, N_new, 'uint8');


        for new_evt = 1: N_new %this is to count the event number in the reduced evt
            evt = final_re_process(new_evt); % only go through events that need reprocessing

            large_pulse = find(dp.pulse_area_phe(:, evt) > min_pulse_area);   %find first large pulse 
            first_large_pulse = large_pulse(1);  
                    %this can only work if the first large pulse is a currently
                    %saved pulse
            for k = 1 : parts 

                 if k == 1 %imposing an s1 size
                    new_pulse_start(k, new_evt) = dp.pulse_start_samples(first_large_pulse, evt);
                    new_pulse_end(k, new_evt) = dp.pulse_start_samples(first_large_pulse, evt) + s1_length; 
                 else
                    new_pulse_start(k, new_evt) = dp.pulse_start_samples(first_large_pulse, evt) + s1_length + ((k-2)/(parts-1)) * (detector_sample_length - s1_length) + 1; % +1 to avoid start/stop overlap 
                    new_pulse_end(k, new_evt) = new_pulse_start(k, new_evt) + 1/(parts-1) * (detector_sample_length - s1_length);
                    %this should divide the detector equally, into the number of 'parts'
                    %NOTE: there might be rounding errors, giving you
                    %slighly higher or lower number of total samples
                 end % if k==1

                 if dp.index_kept_sumpods(first_large_pulse, evt)==1
                    new_index_kept_sumpods(k, new_evt)=1;
                 end

            end % for k

        end %going through all events that have long pulse


        new_num_pulses_found = uint32(sum(new_index_kept_sumpods==1));% this will help visualux


        dp_chop.num_pulses_found = new_num_pulses_found;
        dp_chop.pulse_start_samples = new_pulse_start;
        dp_chop.pulse_end_samples = new_pulse_end;
        dp_chop.index_kept_sumpods = new_index_kept_sumpods;
        %now that we have this info let's save it
        dp_chop.chopped_flag =  uint8(ones(1, size(dp_chop.pulse_start_samples, 2)));

        LUXBinaryRQWriter_framework(settings, dp_chop, chop_filename_rq, data_path_rq, final_re_process, livetime);
        clear dp_chop dp; 
        codepath = strcat(LUXcode_path , '/CppModules/bin/PulseTiming_HeightTiming');
        system([ codepath ' ' chop_filename_evt ' ' data_path_evt ' ' chop_filename_rq ' ' data_path_rq ' 5 ' data_processing_settings_path ' ' lug_iqs_xml_file ]);

        status = PulseQuantities_MinimumSet_split_sim(chop_filename_evt,data_path_evt,chop_filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path);
        status = PulseQuantities_PhotonCounting_split_sim(chop_filename_evt,data_path_evt,chop_filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path);
        status = PulseClassifier_MultiDimensional(chop_filename_evt,data_path_evt,chop_filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path);
        status = S1S2Pairing_Naive(chop_filename_evt,data_path_evt,chop_filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path);
        status = Event_Classification(chop_filename_evt,data_path_evt,chop_filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path);
        status = PositionReconstruction_MercuryI (chop_filename_evt,data_path_evt,chop_filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path);
        status = Corrections_PositionCorrection_split_sim(chop_filename_evt,data_path_evt,chop_filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path);
        status = Corrections_ApplyCorrections(chop_filename_evt,data_path_evt,chop_filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path);
        status = PulseQuantities_TimeSince(chop_filename_evt,data_path_evt,chop_filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path); %deals with etrains, and looks at the time since the last 'big' event
        status = AdditionalFileFormat_SaveMatFile(chop_filename_evt,data_path_evt,chop_filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path);
   
%%%%%%%%%%%%%%%%%%%%%
%{
        %now the dp_chop is saved as a .mat file
        dp_chop = load([data_path_rq '/matfiles/' chop_filename_rq '.mat']);
        %dp_reg = LUXLoadRQ1s_framework(filename_rq, data_path_rq);
        dp = post_dpf_corrections(filename_evt , dp);

    dp_combined = dp;
    clear dp;
    fields = fieldnames(dp_combined);
    fields_length = length(fields);
    for ii = 1: fields_length
        current_field = fields{ii};
        if strcmp(current_field, 'admin') || strcmp(current_field, 'source_filename') ...
                || strcmp(current_field, 'livetime_latch_samples') || strcmp(current_field, 'livetime_end_samples')...
                || ~isfield(dp_chop, current_field)
            continue
        elseif ndims(dp_combined.(current_field))>2
            dp_combined.(current_field)(:, :, final_re_process) = dp_chop.(current_field);
        else
            dp_combined.(current_field)(:, final_re_process) = dp_chop.(current_field);
        end
    end
    
    dp_combined.chopped_flag = zeros(1, length(dp_combined.file_number));
    dp_combined.chopped_flag(final_re_process) = 1;

    %save it
    new_filename_rqmat = strrep(filename_rq , '.rq' , '_comb_chop.rq.mat');

    matfiles_dir = dir([data_path_rq '/matfiles_comb_chop/']);
    if isempty(matfiles_dir)
        mkdir(data_path_rq, '/matfiles_comb_chop/');
    end

    save([data_path_rq '/matfiles_comb_chop/' new_filename_rqmat] , 'dp_combined');
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('Done!\n') 
    end % if re_process is empty


    %now something to clean up the intermediary files.
    if exist ('chop_filename_cvt' , 'var') 
        new_folder = [extractBefore(data_path_rq, 'matfiles') filesep 'chopped_events'];
        mkdir (new_folder)
        mkdir ([new_folder '/matfiles/'])
        movefile (([data_path_rq filesep  chop_filename_rq]) , (new_folder))
        movefile(([data_path_rq  '/matfiles/' chop_filename_rq '.mat']), ([new_folder '/matfiles/']));
        rmdir ([data_path_rq  '/matfiles/'])
        delete ([data_path_evt filesep chop_filename_evt '.mat']);
        delete ([data_path_evt filesep chop_filename_cvt]);
        clear chop_filename_cvt chop_filename_rq chop_filename_evt
    end
    
    clear is_long_pulse final_re_process livetime
    clear new_index_kept_sumpods new_pulse_end new_pulse start
    clear new_num_pulses_found
end % ii over all the rq files
toc
'end'

exit
end
