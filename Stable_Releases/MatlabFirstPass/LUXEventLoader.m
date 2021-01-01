function [event_struct settings] = LUXEventLoader(filename_prefix, data_path, event_list, settings)
% [event_struct settings] = LUXEventLoader(filename_prefix, data_path, event_list, settings)
%
% inputs: filename_prefix - the dataset you are loading, or the filename that you want to load.
%         data_path - the path to the data containing the file(s) you want
%         to load, or the dataset folder
%         event_list - list of events to load, if there is an index file.
%
%
% outputs: event_struct - structure containing the data
%          settings - the settings structure from the head of the file.
%
% example usage:
% [event_struct settings] = LUXEventLoader_catxml_v4('lux01_20090804T1253,'./',1:100)
% [event_struct settings] = LUXEventLoader_catxml_v4([],./',[]) % will
% load from latest dataset in ./, or if it can't find any, load the .evt files
% [event_struct settings] = LUXEventLoader_catxml_v4('lux01_20090804T1253_f000000001.evt,'./',[])
% 
% 
% 2010-04-15 JJC v5.0
% Load everything if event_list isn't specified
% 2010-11-08 JJC v5.1 replaced all dual_XXX with XXX because dual was
% crashing.


load_file = 0;
if nargin < 3
    event_list = [];
end

if ~exist('data_path','var') || isempty(data_path)
    fprintf('you did not specify a data path. assuming ./\n');
    data_path = './';
end

if ~isempty(filename_prefix) && (length(filename_prefix)>19) % data file
    if (length(filename_prefix)==34) && strcmp(filename_prefix(31:34),'.evt')
        load_file=1;
        filename_list = dir([data_path,filesep,filename_prefix]);
        if isempty(filename_list)
            fprintf('cannot find file %s in %s\naborting...\n',filename_prefix,data_path);
            return;
        end
    else
        fprintf('%s is not the name of a dataset or the name of a .evt file.\n',filename_prefix);
        return;
    end
end

if isempty(filename_prefix)
    fprintf('Looking for the latest dataset in %s ...\n', data_path);
    latest_filename_list = {};
    latest_filename_n = [];
    
    % Making assumption here -- all files / datasets are to start with 'lux' prefix
    % Note -- inclusive with LUX 0.1 files, so this is OK
    filename_list = dir([data_path,filesep,'lux*']);
    if ~isempty(filename_list)
        % Check whether the folder listed in filename_list fits the right format
        ii_file = 1;
        while ii_file < length(filename_list)
            if filename_list(ii_file).isdir
                % test to see if bit after underscore is a well-formed date
                % string
                [prefix, rem] = strtok(filename_list(ii_file).name, '_');
                datestring = strtok(rem, '_');
                try
                    datenum = datestr(datestring, 'yyyymmddTHHMM');
                    ii_file = ii_file+1;
                catch
                    filename_list(ii_file) = [];
                end
            else
                filename_list(ii_file) = [];
            end
        end
    end
    
    if isempty(filename_list) || ~any([filename_list.isdir]) % no results found or left after name checking
        fprintf('There are no datasets in %s\nlooking for .evt files instead\n', data_path);
        filename_list = dir([data_path,filesep,'*.evt']);
        if isempty(filename_list)
            fprintf('There are also no .evt files in %s\naborting...\n',data_path);
            return;
            
        else
            load_file=1;
        end
    else
        if load_file==0
            dates = [filename_list.datenum];
            [tmp, ii_latest] = max(dates);
            filename_prefix = filename_list(ii_latest).name;
            fprintf('Loading settings from dataset %s\n', filename_prefix);
        end
    end
    
end

if load_file==0
    %%%% Check for Index File %%%%%
    filename_ind = [data_path, '/', filename_prefix, '.ind' ];
    if ~exist(filename_ind)
        fprintf('no index file found. you need to find the index file,\nor run the index file generator.');
        return;
    end
    
    indexfid = fopen(filename_ind, 'rb', 'l');
    
    if isempty(event_list)
        offset = -7*4; % skip to beginning of last event entry
        fseek(indexfid, double(offset), 'eof');
        event_gid = fread(indexfid, 4, 'uint32');
        event_list = 1:(event_gid(3)+1);
    end
    filename_open = '';
    for evt=event_list
        
        offset = 4*(1+(evt-1)*7);
        fseek(indexfid,offset,'bof');
        index_entry = fread(indexfid, 7, 'ulong');
        event_gid = index_entry(1:4);
        record_format = index_entry(5);
        file_number = index_entry(6);
        file_pos = index_entry(7);
        filename_evt = [filename_prefix sprintf('_f%09d.evt',file_number)];
        if ~strcmp(filename_evt,filename_open)
            if exist('eventfid','var')
            fclose(eventfid);
            end
        
        options.lug_query=0;
        
        if evt==event_list(1) && ~exist('settings','var');
            options.load_settings=1;
            [settings, eventfid] = LUXSuperLoader(filename_evt,data_path,options);
            filename_open = filename_evt; 
        else
            options.load_settings=0;
            [settings_tmp, eventfid] = LUXSuperLoader(filename_evt,data_path,options);
            filename_open = filename_evt;
        end
        event(evt).filename_prefix = filename_prefix; % added 2010-02-03 JJC
            event(evt).pretrigger = settings.evt_settings.pretrigger;
            event(evt).posttrigger = settings.evt_settings.posttrigger;
        
        cur_pos = ftell(eventfid);
        file_header = fread(eventfid,5,'uint32');
        vrange_bot = settings.daq_settings.sis3301.global.vrange_bot * 1000;
        end
        empty=0;
        fseek(eventfid, double(file_pos)+double(cur_pos), 'bof');
        if (settings.daq_settings.global.event_builder_version==4) % added to fix different index file position in v4.0 20090321 jjc
        fread(eventfid, 1, 'uint32');%1
        fread(eventfid, 1, 'uint32'); %2
        fread(eventfid, 1, 'uint32'); %evt
        fread(eventfid, 1, 'uint32'); %nb_chs
        fread(eventfid, 1, 'uint32'); %eventsize
        fread(eventfid, 1, 'uint32'); %recordformat
        fread(eventfid, 1, 'uint32'); %recordsize
    end
    %xml_settings.global.event_builder_version = 1;
    %xml_settings.global.daq_version = 6.5 ;
    %if ((settings.daq_settings.global.event_builder_version>=3) && (settings.daq_settings.global.daq_version>=6.2));
        event(evt).timestamp = fread(eventfid, 1, 'uint64');
    %elseif ((settings.daq_settings.global.event_builder_version>=2) && (settings.daq_settings.global.daq_version>=6.2));
    %    event(evt).timestamp = fread(eventfid, 1, 'uint32');
    %elseif ((settings.daq_settings.global.event_builder_version>=1) || (settings.daq_settings.global.daq_version<6.2));
    %    event(evt).timestamp = 0;
    %end

    %timestamp = event(evt).timestamp
    evt;
    nb_chs = event_gid(4);
    
    for ch=1:nb_chs
        ch;
        bin_datatype = fread(eventfid, 1, 'int8');
        v_res = fread(eventfid, 1, 'double');
        v_off = fread(eventfid, 1, 'double');
        t_res = fread(eventfid, 1, 'double');
        pretrigger = fread(eventfid, 1, 'int32');
        total_samples = fread(eventfid, 1, 'uint32');
        pulse_detect_pretrigger = fread(eventfid, 1, 'uint32');
        pulse_end_posttrigger = fread(eventfid, 1, 'uint32');
        nb_pulses = fread(eventfid, 1, 'uint32');

        if nb_pulses ~= 0
            pulse_start = fread(eventfid, double(nb_pulses), '*int32');
            pulse_length = fread(eventfid, double(nb_pulses), '*uint32');
            pulse_baseline = fread(eventfid, double(nb_pulses), '*uint32');

            for ps=1:nb_pulses
                event(evt).ch(ch).pulse(ps).start = pulse_start(ps);
                event(evt).ch(ch).pulse(ps).length = pulse_length(ps);
                event(evt).ch(ch).pulse(ps).baseline = pulse_baseline(ps);
                pulse_data = fread(eventfid, double(pulse_length(ps)), '*uint16');
                event(evt).ch(ch).pulse(ps).pulse_data = pulse_data;
                event(evt).ch(ch).pulse(ps).pulse_data_mV = double(pulse_data) .* (2000/(2^14)) + vrange_bot; % dynamic range 0.1 -> -1.9
                if ~isempty(event(evt).ch(ch).pulse(ps).pulse_data)
                    event(evt).ch(ch).pulse(ps).baseline_mV = double(event(evt).ch(ch).pulse(ps).baseline) * (2000/(2^14)) + vrange_bot; % 2009-07-15 jjc - now use baseline value from struck
                %event(evt).ch(ch).pulse(ps).baseline_mV = double(max(event(evt).ch(ch).pulse(ps).pulse_data(2:4))) * (2000/(2^14)) + vrange_bot;
                else
                    event(evt).ch(ch).pulse(ps).baseline_mV = 0;
                end
            end

        end
        if nb_pulses == 0
            empty=empty+1;
        end
    end
    event(evt).empty=0;
    if empty == nb_chs % if all channels are empty, then event is empty
        event(evt).empty=1;
    end
    end
    event_struct = event;
end





if load_file==1 % load files
    event_list_input = event_list;
    for ii=1:length(filename_list)
        %evt = evt+1;
        
        filename_evt = filename_list(ii).name;
        %eventfid = OpenBinaryFile(filename_evt, 1);
        %filename_open = filename_evt;
        %endiancheck = fread(eventfid,1,'uint32');
        %settings_length = fread(eventfid,1,'uint32');
        %settings_string = fread(eventfid,settings_length,'char');
        %settings = XMLParser(char(settings_string'));
        
        options.lug_query=0;
        if ii==1 && ~exist('settings','var');
            options.load_settings=1;
            [settings, eventfid] = LUXSuperLoader(filename_evt,data_path,options);
        else
            options.load_settings=0;
            [settings_tmp, eventfid] = LUXSuperLoader(filename_evt,data_path,options);
        end
        file_header = fread(eventfid,5,'uint32');
        vrange_bot = settings.daq_settings.sis3301.global.vrange_bot * 1000;
        
        
        event_list = (file_header(4)+1):((file_header(4))+file_header(5));
        
        
        %if (settings.daq_settings.global.event_builder_version>=4.3 && settings.daq_settings.global.daq_version>=6.8)
            nb_seqs_in_file = fread(eventfid,1,'uint16');
            for i_seq=1:nb_seqs_in_file
                timestamp_latch(i_seq) = fread(eventfid,1,'uint64');
                timestamp_end(i_seq) = fread(eventfid,1,'uint64');
                
            end
        %end
        
        for evt=event_list
            empty=0;
            event(evt).filename_prefix = filename_evt(1:19); % added 2010-02-03 JJC
            event(evt).pretrigger = settings.evt_settings.pretrigger;
            event(evt).posttrigger = settings.evt_settings.posttrigger;
            
            event_gid = fread(eventfid,7,'uint32');
            nb_chs = event_gid(4);
            %if ((settings.daq_settings.global.event_builder_version>=3) && (settings.daq_settings.global.daq_version>=6.2));
                event(evt).timestamp = fread(eventfid, 1, 'uint64');
            %elseif ((settings.daq_settings.global.event_builder_version>=2) && (settings.daq_settings.global.daq_version>=6.2));
            %    event(evt).timestamp = fread(eventfid, 1, 'uint32');
            %elseif ((settings.daq_settings.global.event_builder_version>=1) || (settings.daq_settings.global.daq_version<6.2));
            %    event(evt).timestamp = 0;
            %end
            
            for ch=1:nb_chs
                ch;
                bin_datatype = fread(eventfid, 1, 'int8');
                v_res = fread(eventfid, 1, 'double');
                v_off = fread(eventfid, 1, 'double');
                t_res = fread(eventfid, 1, 'double');
                pretrigger = fread(eventfid, 1, 'int32');
                total_samples = fread(eventfid, 1, 'uint32');
                pulse_detect_pretrigger = fread(eventfid, 1, 'uint32');
                pulse_end_posttrigger = fread(eventfid, 1, 'uint32');
                nb_pulses = fread(eventfid, 1, 'uint32');
                if nb_pulses ~= 0
                    pulse_start = fread(eventfid, double(nb_pulses), '*int32');
                    pulse_length = fread(eventfid, double(nb_pulses), '*uint32');
                    pulse_baseline = fread(eventfid, double(nb_pulses), '*uint32');
                    %keyboard
                    for ps=1:nb_pulses
                        event(evt).ch(ch).pulse(ps).start = pulse_start(ps);
                        event(evt).ch(ch).pulse(ps).length = pulse_length(ps);
                        event(evt).ch(ch).pulse(ps).baseline = pulse_baseline(ps);
                        pulse_data = fread(eventfid, double(pulse_length(ps)), '*uint16');
                        event(evt).ch(ch).pulse(ps).pulse_data = pulse_data;
                        event(evt).ch(ch).pulse(ps).pulse_data_mV = double(pulse_data) .* (2000/(2^14)) + vrange_bot; % dynamic range 0.1 -> -1.9
                        if ~isempty(event(evt).ch(ch).pulse(ps).pulse_data)
                            event(evt).ch(ch).pulse(ps).baseline_mV = double(event(evt).ch(ch).pulse(ps).baseline) * (2000/(2^14)) + vrange_bot; % 2009-07-15 jjc - now use baseline value from struck
                            %event(evt).ch(ch).pulse(ps).baseline_mV = double(max(event(evt).ch(ch).pulse(ps).pulse_data(2:4))) * (2000/(2^14)) + vrange_bot;
                        else
                            event(evt).ch(ch).pulse(ps).baseline_mV = 0;
                        end
                    end
                    
                end
                if nb_pulses == 0
                    empty=empty+1;
                end
            end
            event(evt).empty=0;
            if empty == nb_chs % if all channels are empty, then event is empty
                event(evt).empty=1;
            end
            
        end
        fclose(eventfid);
        
    end
    event_struct = event;
end







































