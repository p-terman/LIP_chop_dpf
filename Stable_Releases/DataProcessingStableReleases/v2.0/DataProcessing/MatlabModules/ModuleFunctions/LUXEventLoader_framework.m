function [event_struct settings] = LUXEventLoader_framework(data_path, arg2, options)
% [event_struct settings] = LUXEventLoader_framework(data_path, arg2, options)
%
%     if arg2 is filename, then load that file. if it is a number or range,
%     look in the data_path for those numbers
%
%     options is .load_xenon_chs (default 1)
%                .load_water_chs (default 0)
%                .load_xlm_info (default 1)
%                .load_tte_ch (default 0)
%
% Example: [event_struct settings] = LUXEventLoader_framework('~/path/to/data/lux10_20120402T1010/', 'lux10_20120402T1010_f000000001.evt');
%
%
% VERSIONING:
%
% 2011-08-16 JJC re-written to work with index in file header. See
%                LUXEventLoader_old for previous history.
% 2011-08-22 JJC preallocated structs (+ other mods) to be slightly (~25%)
% faster.
% 2011-08-24 JJC if channel is empty, then event_struct.ch.pulse = [];
% 2011-09-14 JJC fixed bug that was loading wrong events, fixed empty flag
%                counting trigger channel.
% 2012-01-04 JJC now loads .evt.gz files with ease.
% 2012-02-14 JJC above now works if you don't have permissions to write to
%                data_path. File is copied to current working directory and deleted
%                afterwards.
% 2012-02-17 JJC modified for .pod instead of .pulse
% 2012-03-06 JJC fread for 2009a works - cast size as double
% 2012-05-29 JJC uses NewGetFiles2Load which "intelligently" searches files
%                rather than looking in every damn file starting at 1. Also sorts output
%                based on event number
% 2012-09-12 JJC now adds analog sum channels to the structure.
% 2012-09-18 JJC looks for .ind file, loads if it exists and you are
%                looking for a specific event list (rather than a file or all events).
% 2012-09-18 JRV Now uses /tmp/ instead of ./ for unzipping files
% 2012-09-19 JJC Fixed indexing evt_list in files to load.
% 2012-09-20 JJC Added options input with fields .load_only_analog_sum_chs
%                and .load_only_xlm_info.
% 2012-11-05 JJC If index file is empty, will read the headers of files for
%                index - (why the hell would the index file be empty? still
%                investigating...)
% 2012-12-12 JJC - replaced pmt_chs with xenon_daq_chs
% 2012-12-13 JJC - If arg2 is a negative number then load that many of the
%                  latest events based on the index file.
% 2012-12-18 JJC - XLM update: trigseqnum, ccheck. backwards compatible.
% 2013-01-03 JJC - oops - don't use 'char' as a datatype, use 'uint8'
%                  instead.
% 2013-01-10 JJC - improving options structure above to be able to load
%                  xenon_chs and/or water_chs and/or analog_sum and/or tte_ch and/or xlm_info
% 2013-02-04 CHF - renamed baseline_mV to pod_baseline_mV
%                  Now mV is a positive pulse on a positive scale, and
%                  pod_baseline_mV is a positive quantity.
%                  Formatted this header to be more readable.
%                  Replaced (.122) for mV_per_adc, defined as 2000/2^14
% 2013-03-18 CHF - Minor fix with pulse_start_samples (was off by 1 sample)
% 2013-04-15 CHF - Now can load zero event files. It simply loads an empty
%                  event_struct, but the settings are still given as
%                  output.
% 2013-06-17 pfs - added ability to read in a single event (previously required read of entire file)
%
%
%% Constants

mV_per_adc = 2000./2^14;

%% Initialize and check inputs

event_struct = [];
event = [];
settings = [];
if nargin<2
    arg2 = [];
end

if ~exist('options','var')
    options.load_analog_sum_chs=0;
    options.load_xlm_info=0;
    options.load_xenon_chs=1;
    options.load_water_chs=0;
    options.load_tte_ch=0;
else
    if ~isfield(options, 'load_analog_sum_chs')
        options.load_analog_sum_chs=0;
    end
    if ~isfield(options, 'load_xlm_info')
        options.load_xlm_info=0;
    end
    if ~isfield(options, 'load_xenon_chs')
        options.load_xenon_chs=1;
    end
    if ~isfield(options, 'load_water_chs')
        options.load_water_chs=0;
    end
    if ~isfield(options, 'load_tte_ch')
        options.load_tte_ch=0;
    end
end

gunzip_loc = '/tmp/';
[direxst, fa] = fileattrib(gunzip_loc);

try
    % check if /tmp exists and has proper permissions
    if ~(direxst && ...
            fa.UserRead == 1 && ...
            fa.UserWrite == 1 && ...
            fa.GroupRead == 1 && ...
            fa.GroupWrite == 1 && ...
            fa.OtherRead == 1 && ...
            fa.OtherWrite == 1)
        
        % /tmp doesn't exist or isn't properly writable
        fprintf('Warning: /tmp is not writable. Using ./ for temp files.\n');
        gunzip_loc = './';
    end
catch
    % don't use Windows...
    fprintf('WARNING: Can not use /tmp directory. Using ./ for temp files.\n');
    gunzip_loc = './';
end


%% Decide what to load

if isnumeric(arg2) && ~isempty(arg2)
%    flag_search_events = 1;
    arg2=sort(arg2);
    index_file = dir([data_path '/*.ind']);
    if ~isempty(index_file)
        files_to_load = GetFiles2LoadFromIndex(data_path, arg2);
    	filename = files_to_load(1).name;
        evt_list = []; % need this, don't know why
    else
		disp('*** missing .ind file. please go get it.');return
    end
    
end
%event_struct.file_name = files_to_load.name; % works for 1 file

if ischar(arg2)
    filename = arg2;
    evt_list = [];
end

%%

ii = 0;

options.load_settings=1;

temp_loc=0;

% Gunzip if needed
if strcmp(filename((end-1):end),'gz')
    
    fprintf('gunzipping file %d...\n',str2num(filename(22:30)));
    
    if ~exist([gunzip_loc '/' filename],'file')
        copyfile([data_path '/' filename], gunzip_loc);
    end
    
    temp_loc=1;
    gunzip([gunzip_loc filename]);
    filename = filename(1:(end-3));
    
end

% Load settings
if temp_loc==1
    [settings, eventfid] = LUXSuperLoader_framework(filename,gunzip_loc,options);
else
    [settings, eventfid] = LUXSuperLoader_framework(filename,data_path,options);
end

if settings.evt_settings.event_builder_version<7
    fprintf('You need to rebuild this data with the new event builder, or use LUXEventLoader_old.\n');
    return;
end

vrange_bot = settings.daq_settings.sis3301.global.vrange_bot * 1000;
%data_start = ftell(eventfid);
file_header = fread(eventfid,4,'uint32');
nb_evts_in_file = file_header(4);

% SKIP THE REST IF NO EVENTS IN FILE. 20130415 CHF
if nb_evts_in_file > 0
    
    evt_index_block = fread(eventfid,2*nb_evts_in_file,'uint32');
    events_in_file_numbers = evt_index_block(logical(mod(1:(2*nb_evts_in_file),2)));
    events_in_file_positions = evt_index_block(~logical(mod(1:(2*nb_evts_in_file),2)));
    
    if isempty(evt_list)
        evt_list = events_in_file_numbers';
    else
        evt_list = intersect(evt_list,events_in_file_numbers);
    end
    
    % Initialize size of event structure
    event(1).ch_map = [];
    event(length(evt_list)).ch_map = [];
    event(1).filename_prefix = filename(1:19);
    event(length(evt_list)).filename_prefix = filename(1:19);
    event(1).pretrigger = settings.evt_settings.pretrigger;
    event(length(evt_list)).pretrigger = settings.evt_settings.pretrigger;
    event(1).posttrigger = settings.evt_settings.posttrigger;
    event(length(evt_list)).posttrigger = settings.evt_settings.posttrigger;
    
    % Initialize settings
    if isfield(settings.daq_settings.sis3301.global,'pmt_chs')
        event(1).ch_map = settings.daq_settings.sis3301.global.pmt_chs;
    elseif isfield(settings.daq_settings.sis3301.global,'xenon_daq_chs') % 2012-12-12 JJC
        if options.load_xenon_chs
            event(1).ch_map = [event(1).ch_map settings.daq_settings.sis3301.global.xenon_daq_chs];
        end
        if options.load_water_chs
            event(1).ch_map = [event(1).ch_map settings.daq_settings.sis3301.global.water_daq_chs];
        end
        if options.load_tte_ch
            if isfield(settings.daq_settings.sis3301.global,'tte_ch')
                event(1).ch_map = [event(1).ch_map settings.daq_settings.sis3301.global.tte_ch];
            else
                event(1).ch_map = [event(1).ch_map 126];
            end
        end
    else
        fprintf('no xenon_daq_chs in settings, so just assuming ch_map is 1:122\n');
        event(1).ch_map = 1:122;
    end
    
    % Loop for the event list
    for ee = evt_list
        
        if isempty(ee), continue; end
        
        nb_pulses_tot=0;
        
        ii=ii+1;
        
        event(ii).event_number = ee;
        ei = events_in_file_numbers == ee;
        fseek(eventfid,events_in_file_positions(ei),'bof');
        
        if ~isfield(settings.daq_settings.sis3301.global, 'read_xlm')
            settings.daq_settings.sis3301.global.read_xlm=0;
        end
        
        if settings.daq_settings.sis3301.global.read_xlm==0
            event_gid = fread(eventfid,9,'*uint32');
            %record_format = fread(eventfid,1,'uint32');
            %record_size = fread(eventfid,1,'uint32');
            %event(ii).timestamp = fread(eventfid,1,'uint64');
            event(ii).timestamp = typecast(event_gid(8:9),'uint64');
            nb_chs = event_gid(4);
        end
        
        % ---- TPB_140401: Determine just how we'll compute the ADC baseline.
        % Initialize to true a flag that tells the program how the baseline will be computed. 
        % If set to 1, compute the baseline from the first few samples of pod_data. If set to
        % 0, compute the baseline from the DAQ provided value.
        do_baseline_from_pod_data = 1;
       % We also require at least 
       % 5 (arbitrary choice by TPB) pre-trigger samples of which 3 will be used to compute the
       % mean baseline.
        if ~isfield(settings.daq_settings.sis3301.global, 'pulse_detect_pretrigger')
            % pre-trigger sample number is not given in the settings file
            fprintf('\n\n\n\npulse_detect_pretrigger not found in daq_settings. pod_baseline will be used to compute the POD baseline instead of the pre-trigger samples\n\n\n\n');
            % Use the POD baseline as returned by the DAQ
            do_baseline_from_pod_data = 0; % set the flag to use the DAQ values
        else
            % The pre-trigger length is given. Now make sure that it contains enough samples
            if settings.daq_settings.sis3301.global.pulse_detect_pretrigger < 5
                % if it doesn't contain enough samples
                fprintf('\n\n\n\npulse_detect_pretrigger does not contain enough samples to average. pod_baseline will be used to compute the POD baseline instead of the pre-trigger samples\n\n\n\n');
                % Use the POD baseline as returned by the DAQ
                do_baseline_from_pod_data = 0; % set the flag to use the DAQ values
            % else we'll get te baseline from the data
            end
        end
        % ----
        
        if settings.daq_settings.sis3301.global.read_xlm==1
            event_gid = fread(eventfid,5,'*uint32');
            
            event(ii).trigger_timestamp = fread(eventfid,1,'uint64');
            
            if settings.daq_settings.global.daq_version>7.0
                event(ii).trigseqnum = fread(eventfid,1,'uint32');
            end
            
            event(ii).max_filter_response = fread(eventfid,1,'uint32');
            event(ii).max_ch_ID = fread(eventfid,1,'int8');
            xlm_nb_ddcs = fread(eventfid,1,'uint16');
            
            for d_xlm=1:xlm_nb_ddcs
                event(ii).S1_hit_vector(d_xlm) = fread(eventfid,1,'uint8');
                event(ii).S2_hit_vector(d_xlm) = fread(eventfid,1,'uint8');
            end
            
            if settings.daq_settings.global.daq_version>7.0
                event(ii).xlm_ccheck = fread(eventfid,1,'uint8');
            end
            
            if ~options.load_xlm_info
                rmfield(event(ii),'trigger_timestamp');
                if settings.daq_settings.global.daq_version>7.0; rmfield(event(ii),'trigseqnum'); end
                rmfield(event(ii),'max_filter_response');
                rmfield(event(ii),'max_ch_ID');
                rmfield(event(ii),'S1_hit_vector');
                rmfield(event(ii),'S2_hit_vector');
                if settings.daq_settings.global.daq_version>7.0; rmfield(event(ii),'xlm_ccheck'); end;
            end
            
            record_info = fread(eventfid,2,'uint32');
            event(ii).timestamp = fread(eventfid,1,'uint64');
            nb_chs = event_gid(4);
        end
        
        ch_list= sort(1:max(event(1).ch_map));
        ch_seek=ftell(eventfid);
        
        %event(ii).ch(1:nb_chs) = struct('pod',struct('start',0,'length',0,'baseline',0,'pod_data',[],'baseline_mV',0,'pod_data_mV',[]));
        if options.load_xenon_chs || options.load_water_chs || options.load_analog_sum_chs
            
            for ch = 1:length(ch_list)
                
                if ch==1
                    fseek(eventfid,ch_seek,'bof');
                end
                
                if options.load_analog_sum_chs
                    fseek(eventfid,event_gid(2),'bof');
                end
                
                channel_header_in_bytes = fread(eventfid,45,'*int8');
                nb_pulses = typecast(channel_header_in_bytes(42:45),'uint32');
                
                if any(ch==event(1).ch_map)
                    nb_pulses_tot = nb_pulses_tot+nb_pulses;
                end
                
                if nb_pulses > 0
                    event(ii).ch(ch).empty = 0;
                    event(ii).ch(ch).pod_start_samples = fread(eventfid,double(nb_pulses),'int32');
                    event(ii).ch(ch).pod_length_samples = fread(eventfid,double(nb_pulses),'uint32');
                    event(ii).ch(ch).pod_baseline = fread(eventfid,double(nb_pulses),'uint32');
                    % TPB_140401: We technically no longer need the above pod_baseline since we'll
                    % compute it from the begining of the POD data itself but we must advance the
                    % file pointer anyway so we read it and don't do anything with it. Also, if for some
                    % unholy reason the pre-trigger length isn't recorder or isn't long enough then
                    % we'll actually use event(ii).ch(ch).pod_baseline as the baseline.
                    %
                    %event(ii).ch(ch).pod_baseline_mV = -(double(event(ii).ch(ch).pod_baseline) * mV_per_adc + vrange_bot); % pod_baseline_mV is a POSITIVE quantity.
                    % TPB_140401: The above baseline estimate is biased due to integer division truncation.
                    % Below we use the first few samples of each POD to estimate the baseline instead if the
                    % daq_settings tell us how man pre-trigger samples are available. We also require at least 
                    % 5 (arbitrary choice by TPB) pre-trigger samples of which 3 will be used to compute the
                    % mean baseline.
                    %
                    % read in the data for all pulses in this channel in this event
                    event(ii).ch(ch).pod_data = fread(eventfid,double(sum(event(ii).ch(ch).pod_length_samples)),'*uint16');
                    %
                    % Define the baseline in mV if there aren't enough pre-trigger samples or initialize the
                    % array that will store the baseline if there are enough pre-trigger samples
                    if do_baseline_from_pod_data == 0
                        % Use the POD baseline as returned by the DAQ
                        event(ii).ch(ch).pod_baseline_mV = -(double(event(ii).ch(ch).pod_baseline) * mV_per_adc + vrange_bot); % pod_baseline_mV is a POSITIVE quantity.
                    else
                        % Initialize the array that will hold the computed baselines for each POD
                        event(ii).ch(ch).pod_baseline_mV = zeros(1,length(event(ii).ch(ch).pod_baseline));
                    end
                    
                    % Initialize the array that stores the "time" arrays for each POD
                    event(ii).ch(ch).pod_time_samples = [];
                    
                    % Initialize the array that will store the baslines expanded to the shape of the 
                    % data array. Essentially, a single baseline value will be replicated into a number
                    % of samples to match the length of the pulse it belongs to.
                    baselines_mV_vector = zeros(1,length(event(ii).ch(ch).pod_data));
                    
                    % Loop over all of the PODs in this channel
                    for dd = 1:nb_pulses;
                        % get the length of this POD
                        M = event(ii).ch(ch).pod_length_samples(dd);
                        % define the "time" array for the POD
                        tt = (0:M-1) + event(ii).ch(ch).pod_start_samples(dd);
                        % store the above "time" array by appending it to an existing array
                        event(ii).ch(ch).pod_time_samples = [event(ii).ch(ch).pod_time_samples tt];
                        % get the slice of the data that corresponds to the current time array
                        slicecut = inrange(event(ii).ch(ch).pod_time_samples,event(ii).ch(ch).pod_start_samples(dd),event(ii).ch(ch).pod_start_samples(dd)+M+1);
                        % compute the mean baseline from the data, if possible
                        if do_baseline_from_pod_data == 1
                        	% get the piece of data that corresponds to the POD
                            tempdata = event(ii).ch(ch).pod_data(slicecut);
                            % and compute the baseline from that
                            event(ii).ch(ch).pod_baseline_mV(dd) = -(mean(double(tempdata(1:settings.daq_settings.sis3301.global.pulse_detect_pretrigger-2))) * mV_per_adc + vrange_bot); % pod_baseline_mV is a POSITIVE quantity.
                        end
                        % populate the vector of a shape corresponding to the data shape
                        baselines_mV_vector(slicecut) = event(ii).ch(ch).pod_baseline_mV(dd);
                    end
                    
                    % pod_data_mV is a POSITIVE quantity. We
                    % convert from ADC to mV, offset and flip
                    % the sign. Then we subtract the POSITIVE
                    % pod_baseline_mV
                    event(ii).ch(ch).pod_data_mV = -(double(event(ii).ch(ch).pod_data) * mV_per_adc + vrange_bot) - baselines_mV_vector';
                    
                    % Visualize
                    if 0
                        figure(123); clf
                        plot(event(ii).ch(ch).pod_time_samples,event(ii).ch(ch).pod_data_mV,'k.-');
                        hold on
                        plot(event(ii).ch(ch).pod_time_samples,-(double(event(ii).ch(ch).pod_data) * mV_per_adc + vrange_bot),'b.-');
                        plot(event(ii).ch(ch).pod_time_samples,baselines_mV_vector,'r.-')
                        keyboard
                    end
                    
                    %                 if nb_pulses > 1; keyboard; end
                    
                else
                    event(ii).ch(ch).empty = 1;
                end
                
            end % for ch
        end % fi
        
        if nb_pulses_tot==0
            event(ii).empty=1;
        else
            event(ii).empty=0;
        end
        
        
        %     fprintf('.');
        %
        %     % Print to screen loader status
        %     if mod(ee,100) == 0
        %         fprintf('\n')
        %     end
        
        
    end % for event list
    
    if temp_loc==1
        delete([gunzip_loc filename]);
        delete([gunzip_loc filename '.gz']);
    end
    
    
    % Sort events
    if ~isempty(event)&&~issorted([event(:).event_number])
        [unused order] = sort([event(:).event_number]);
        event_struct = event(order);
    else
        event_struct = event;
    end
    
end

event_struct(1).filename = filename; % wannnnt this (and assume, possibly incorrectly, that all events in event_struct come from same filename)
fclose(eventfid);

fprintf('\n');




