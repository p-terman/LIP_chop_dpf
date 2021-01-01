function [event_struct settings] = LUXCVTLoader_framework(data_path, arg2, options)
% [event_struct settings] = LUXCVTLoader_framework(data_path, arg2, options)
%
%     if arg2 is filename, then load that file. if it is a number or range,
%     look in the data_path for those numbers
%
%     options is .load_xenon_chs (default 1)
%                .load_water_chs (default 0)
%                .load_xlm_info (default 1)
%                .load_tte_ch (default 0)
%
% Example: [event_struct settings] = LUXCVTLoader_framework('~/path/to/data/lux10_20120402T1010/', 'lux10_20120402T1010_f000000001.evt');
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
% 2013-02-17 CHF - Major changes to make it work as CVT loader
% 2013-03-11 pfs - added min_time_samples and max_time_samples, so we can load evt files as
%                  well as cvt files (probably only useful for event viewing)
% 2013-03-18 CHF - Minor fix with pulse_start_samples (was off by 1 sample)
% 2013-03-22 CHF - Fixed sumpod offset readout problem that was plaging CVT files. 
%                  Cleaned up the code in a major way. Added comments.
% 2013-03-22 pfs - modified logic on line 204 so that nb_chs gets defined
%                  regardless (code was bonking due to lack of this variable)
% 2013-04-05 AC  - added filename_prefix from settings.daq_settings 
%                  to root of settings, as expected by LUXBinaryWriter_framework
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

settings.filename_prefix=settings.daq_settings.global.filename_prefix;
vrange_bot = settings.daq_settings.sis3301.global.vrange_bot * 1000;
%data_start = ftell(eventfid);
file_header = fread(eventfid,4,'uint32');
nb_evts_in_file = file_header(4);

evt_index_block = fread(eventfid,2*nb_evts_in_file,'uint32');
events_in_file_numbers = evt_index_block(1:2:end);
events_in_file_positions = evt_index_block(2:2:end);

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

% Loop for the event list
for ee = evt_list
    
    if isempty(ee), continue; end
    
    nb_pulses_total=0;
    
    ii=ii+1;
    
    event(ii).event_number = ee;
    ei = events_in_file_numbers == ee;
    fseek(eventfid,events_in_file_positions(ei),'bof');
    
    if ~isfield(settings.daq_settings.sis3301.global, 'read_xlm')
        settings.daq_settings.sis3301.global.read_xlm=0;
    end
    
    if 1%settings.daq_settings.sis3301.global.read_xlm==0
        event_gid = fread(eventfid,9,'*uint32');
        %record_format = fread(eventfid,1,'uint32');
        %record_size = fread(eventfid,1,'uint32');
        %event(ii).timestamp = fread(eventfid,1,'uint64');
        event(ii).timestamp = typecast(event_gid(8:9),'uint64');
        nb_chs = event_gid(4);
    end
    
    % Ignore all this POO
    %
    %     if settings.daq_settings.sis3301.global.read_xlm==1
    %         event_gid = fread(eventfid,5,'*uint32');
    %
    %         event(ii).trigger_timestamp = fread(eventfid,1,'uint64');
    %
    %         if settings.daq_settings.global.daq_version>7.0
    %             event(ii).trigseqnum = fread(eventfid,1,'uint32');
    %         end
    %
    %         event(ii).max_filter_response = fread(eventfid,1,'uint32');
    %         event(ii).max_ch_ID = fread(eventfid,1,'int8');
    %         xlm_nb_ddcs = fread(eventfid,1,'uint16');
    %
    %         for d_xlm=1:xlm_nb_ddcs
    %             event(ii).S1_hit_vector(d_xlm) = fread(eventfid,1,'uint8');
    %             event(ii).S2_hit_vector(d_xlm) = fread(eventfid,1,'uint8');
    %         end
    %
    %         if settings.daq_settings.global.daq_version>7.0
    %             event(ii).xlm_ccheck = fread(eventfid,1,'uint8');
    %         end
    %
    %         if ~options.load_xlm_info
    %             rmfield(event(ii),'trigger_timestamp');
    %             if settings.daq_settings.global.daq_version>7.0; rmfield(event(ii),'trigseqnum'); end
    %             rmfield(event(ii),'max_filter_response');
    %             rmfield(event(ii),'max_ch_ID');
    %             rmfield(event(ii),'S1_hit_vector');
    %             rmfield(event(ii),'S2_hit_vector');
    %             if settings.daq_settings.global.daq_version>7.0; rmfield(event(ii),'xlm_ccheck'); end;
    %         end
    %
    %         record_info = fread(eventfid,2,'uint32');
    %         event(ii).timestamp = fread(eventfid,1,'uint64');
    %         nb_chs = event_gid(4);
    %     end
    
    event(1).ch_map = 1:nb_chs;
    ch_list = sort(1:max(event(1).ch_map));
    ch_seek=ftell(eventfid);
    
    event(ii).ch(1:nb_chs) = struct('pod',struct('start',0,'length',0,'baseline',0,'pod_data',[],'baseline_mV',0,'pod_data_mV',[]));
    
    % Loop per channel
    for ch = 1:length(ch_list)
        
        if ch==1
            fseek(eventfid,ch_seek,'bof');
        end
        
        % Get header information
        bit_precision = fread(eventfid, 1, 'int8');
        v_range_V = fread(eventfid, 1, 'double');
        v_offset = fread(eventfid, 1, 'double');
        time_step_s = fread(eventfid, 1, 'double');
        pretrigger_samples = fread(eventfid, 1, 'int32');
        temp1 = fread(eventfid, 1, 'uint32');
        n_preceding_samples = fread(eventfid, 1, 'uint32');
        n_following_samples = fread(eventfid, 1, 'uint32');
        
        % Get number of pulses
        nb_pulses = fread(eventfid, 1, 'uint32');
        
        % If we are in the sumpod channel
        if ch == 137
                        
            % ...and it has stuff
            if nb_pulses > 0
                
                event(ii).sumpod_start_samples = fread(eventfid,nb_pulses,'int32');
                event(ii).sumpod_length_samples = fread(eventfid,nb_pulses,'uint32');
                pod_baseline_mV = fread(eventfid,nb_pulses,'uint32'); % not for use
                
                % Catch a read offset before it screws up your memory
                if any(pod_baseline_mV ~= 0)
                    error('*** ERROR: found nonzero pod_baseline_mV (written always as 0 in CVT file). This means there is an offset in the readout, which can cause major problems');
                end
                
                % Load summed data                    
                event(ii).sumpod_data_phe_per_sample = fread(eventfid,sum(event(ii).sumpod_length_samples),'single');
                
                % Construct timing vector
                event(ii).sumpod_time_samples = [];
                for dd = 1:nb_pulses;
                    M = event(ii).sumpod_length_samples(dd);
                    tt = (0:M-1) + event(ii).sumpod_start_samples(dd);
                    event(ii).sumpod_time_samples = horzcat(event(ii).sumpod_time_samples,tt);
                end
                
                event(ii).sumpod_time_samples = event(ii).sumpod_time_samples';
                
                
            end
            
        else
            
            % If it's a channel 1:136
            % ...and it has stuff
            if nb_pulses > 0
                event(ii).ch(ch).empty = 0;
                
                % Increment global counter for this event
                nb_pulses_total = nb_pulses_total + nb_pulses;
                
                event(ii).ch(ch).pod_start_samples = fread(eventfid,nb_pulses,'int32');
                event(ii).ch(ch).pod_length_samples = fread(eventfid,nb_pulses,'uint32');
                event(ii).ch(ch).pod_baseline_mV = fread(eventfid,nb_pulses,'uint32');
                
                % Catch a read offset before it screws up your memory
                if any(event(ii).ch(ch).pod_baseline_mV ~= 0)
                    error('*** ERROR: found nonzero pod_baseline_mV (written always as 0 in CVT file). This means there is an offset in the readout, which can cause major problems');
                end
                
                % Read the data
                event(ii).ch(ch).pod_data_phe_per_sample = fread(eventfid,sum(event(ii).ch(ch).pod_length_samples),'single');
                
                % Construct timing vector
                event(ii).ch(ch).pod_time_samples = [];
                for dd = 1:nb_pulses;
                    M = event(ii).ch(ch).pod_length_samples(dd);
                    tt = (0:M-1) + event(ii).ch(ch).pod_start_samples(dd);                    
                    event(ii).ch(ch).pod_time_samples = [event(ii).ch(ch).pod_time_samples tt];
                end

                % Visualize for troubleshooting
                if 0
                    figure(123); clf
                    plot(event(ii).ch(ch).pod_time_samples,event(ii).ch(ch).pod_data_mV,'k.-');
                    hold on
                    plot(event(ii).ch(ch).pod_time_samples,-(double(event(ii).ch(ch).pod_data) * mV_per_adc + vrange_bot),'b.-');
                    plot(event(ii).ch(ch).pod_time_samples,baselines_mV_vector,'r.-')
                    keyboard
                end
            else
                event(ii).ch(ch).empty = 1;
            end
        end
    end % for ch
    
    % Check if global counter had any events in it. Label event as empty if not
    if nb_pulses_total==0
        event(ii).empty=1;
    else
        event(ii).empty=0;
    end
        
end % for event list

% Close file
fclose(eventfid);

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

fprintf('\n');


