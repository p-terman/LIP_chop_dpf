function [event_struct] = LUXPMTEventLoader(xml_settings, event_list)
% [event_struct] = LUXPMTEventLoader(xml_settings, event_list)
%
% Inputs:
%	* xml_settings:	structure with the dataset settings, created with
%                   LUXLoadSettings
%   * event_list: list of event numbers you want to load.
%                 for example, 1:100, 1, 2, 3, 67.
%
% Outputs:
%   * event_struct: struct array of info for all events
%                 .pretrigger: number of pretrigger (samples? us?)
%                 .posttrigger: number of posttrigger samples (samples?
%                   us?)
%                 .timestamp: timestamp for this event
%                 .ch: struct array of all channel info for given event
%                    .pulse: struct array of pulses for given channel
%                          .start: start location of pulse
%                          .length: length of pulse
%                          .baseline: baseline of pulse (ADC counts)
%                          .baseline_mV: baseline of pulse (mV)
%                          .pulse_data: pulse data (ADC counts)
%                          .pulse_data_mV: pulse data (mV)
%
%
% 2008-06-25 jjc
% 2008-07-15 ldv - renamed it from event_loader to LUXEventLoader
% 2008-08-04 ldv - modified it to use xml_settings as an input
% 2009-02-11 dcm - added some help text for the output
% 2009-05-13 dcm - changed .evt loading to use OpenBinaryFile, all file 
%  reads from f* functions now use dual_f*. Event loop only calls
%  OpenBinaryFile when necessary; previously called fopen for every event
% 2009-05-31 dcm - loads all events if event_list isn't specified (i.e.
%  function is called as [event_struct] = LUXEventLoader(xml_settings))
% 2009-06-16 jjc - can now read events with livetime information in the header
% 2009-07-15 jjc - fixed baseline calculation, now explicitly use value from struck
% 20091030 CHF - changed name to LUXPMTEventLoader, using averaging for baseline instead of STRUCK baseline


% Load everything if event_list isn't specified
if nargin < 2
    event_list = [];
end


%%%% Load Settings from xml_settings %%%%%
filename_prefix = xml_settings.filename_prefix;
data_path = xml_settings.data_path;
mode = isfield(xml_settings.sis3301.global, 'nb_evts_per_file');
vrange_bot = xml_settings.sis3301.global.vrange_bot * 1000; % should be 1900 mV


%%%% Check for Index File %%%%%
filename_ind = [data_path, '/', filename_prefix, '/', filename_prefix, '.ind' ];
while ~exist(filename_ind)
    disp(sprintf('Index File Not Found'));
    disp(sprintf('You need to run the Event Builder'));
    disp(sprintf('Do you want to run the Event Builder now?'));
    disp(sprintf('(Note that it will only run in the DAQ computer)'));
    evt_build_answer = input('[Yes]/No','s');
    if strcmp(evt_build_answer,'Y') || strcmp(evt_build_answer,'Yes') || strcmp(evt_build_answer,'y') || strcmp(evt_build_answer,'yes') || isempty(evt_build_answer)
        event_builder_file = '/home/luxdaq/daq/LUX_Event_Builder/LUXEventBuilder';
        if ~exist(event_builder_file,'file')
            disp(sprintf('The event builder could not be located in this machine.  Please run it manually.'));
            disp(sprintf('%s',event_builder_file))
            return;
        else
            unix(sprintf('%s %s %s',event_builder_file,data_path,filename_prefix));
        end
    else
        return
    end
end



%%%% Load the Index File %%%%%
indexfid = fopen(filename_ind, 'rb', 'l');

if isempty(event_list)
    offset = -7*4; % skip to beginning of last event entry
    fseek(indexfid, double(offset), 'eof');
    event_gid = fread(indexfid, 4, 'uint32');
    event_list = 1:(event_gid(3)+1);
end

%tic
filename_open = '';
for evt=event_list
    % for evt=1:1
    
    if mod(evt,100)==0;
        percentage_done = (evt-event_list(1))/length(event_list)*100;
        disp(sprintf('%0.1f%% loaded',percentage_done));
    end
    event(evt).pretrigger = 0;
    event(evt).posttrigger = 200;
    if isfield(xml_settings.sis3301.global, 'pretrigger');
        event(evt).pretrigger = xml_settings.sis3301.global.pretrigger;
        event(evt).posttrigger = xml_settings.sis3301.global.posttrigger;
    end

    offset = 4*(1+(evt-1)*7); %skip to beginning of event in index file
    fseek(indexfid, offset, 'bof');
    event_gid = fread(indexfid, 4, 'uint32');

    if (event_gid(3) ~= (evt-1))
        disp('wrong event number!');
    end
    %nb_chs = event_gid(4)
    nb_chs = 8;
    record_format = fread(indexfid, 1, '*uint32');
    file_number = fread(indexfid, 1, '*uint32');

    if (xml_settings.global.event_builder_version==4) % added to fix C file number starting at 0 in version 4.0 20090321 jjc
        file_number = file_number +1;
    end

    file_pos = fread(indexfid, 1, '*uint32');
    filename_evt = [data_path,'/',filename_prefix,'/',filename_prefix,sprintf('_f%09d',file_number),'.evt'];

    %disp(filename_evt);
    if ~strcmp(filename_evt, filename_open)
        if exist('eventfid','var')
            dual_fclose(eventfid);
        end
        eventfid = OpenBinaryFile(filename_evt, 1);
        filename_open = filename_evt;
    end
    dual_fread(eventfid, 1, 'uint32'); %0
    dual_fread(eventfid, 1, 'uint32'); %1
    dual_fread(eventfid, 1, 'uint32'); %2
    dual_fread(eventfid, 1, 'uint32'); %3
    dual_fread(eventfid, 1, 'uint32'); %4
    
    
    if (xml_settings.global.event_builder_version>=4.3 && xml_settings.global.daq_version>=6.8)
    	nb_seqs_in_file = dual_fread(eventfid,1,'uint16');
    	for i_seq=1:nb_seqs_in_file
    		timestamp_latch(i_seq) = dual_fread(eventfid,1,'uint64');
    		timestamp_end(i_seq) = dual_fread(eventfid,1,'uint64');
    		
    	end
    end
    empty = 0;
    file_pos;
    dual_fseek(eventfid, double(file_pos), 'bof');
    if (xml_settings.global.event_builder_version==4) % added to fix different index file position in v4.0 20090321 jjc
        dual_fread(eventfid, 1, 'uint32');%1
        dual_fread(eventfid, 1, 'uint32'); %2
        dual_fread(eventfid, 1, 'uint32'); %evt
        dual_fread(eventfid, 1, 'uint32'); %nb_chs
        dual_fread(eventfid, 1, 'uint32'); %eventsize
        dual_fread(eventfid, 1, 'uint32'); %recordformat
        dual_fread(eventfid, 1, 'uint32'); %recordsize
    end
    %xml_settings.global.event_builder_version = 1;
    %xml_settings.global.daq_version = 6.5 ;
    if ((xml_settings.global.event_builder_version>=3) && (xml_settings.global.daq_version>=6.2));
        event(evt).timestamp = dual_fread(eventfid, 1, 'uint64');
    elseif ((xml_settings.global.event_builder_version>=2) && (xml_settings.global.daq_version>=6.2));
        event(evt).timestamp = dual_fread(eventfid, 1, 'uint32');
    elseif ((xml_settings.global.event_builder_version>=1) || (xml_settings.global.daq_version<6.2));
        event(evt).timestamp = 0;
    end
    

    %timestamp = event(evt).timestamp
    evt;
    for ch=1:nb_chs
        ch;
        bin_datatype = dual_fread(eventfid, 1, 'int8');
        v_res = dual_fread(eventfid, 1, 'double');
        v_off = dual_fread(eventfid, 1, 'double');
        t_res = dual_fread(eventfid, 1, 'double');
        pretrigger = dual_fread(eventfid, 1, 'int32');
        total_samples = dual_fread(eventfid, 1, 'uint32');
        pulse_detect_pretrigger = dual_fread(eventfid, 1, 'uint32');
        pulse_end_posttrigger = dual_fread(eventfid, 1, 'uint32');
        nb_pulses = dual_fread(eventfid, 1, 'uint32');

        if nb_pulses ~= 0
            pulse_start = dual_fread(eventfid, double(nb_pulses), '*int32');
            pulse_length = dual_fread(eventfid, double(nb_pulses), '*uint32');
            pulse_baseline = dual_fread(eventfid, double(nb_pulses), '*uint32');

            for ps=1:nb_pulses
                event(evt).ch(ch).pulse(ps).start = pulse_start(ps);
                event(evt).ch(ch).pulse(ps).length = pulse_length(ps);
                event(evt).ch(ch).pulse(ps).baseline = pulse_baseline(ps);
                pulse_data = dual_fread(eventfid, double(pulse_length(ps)), '*uint16');
                event(evt).ch(ch).pulse(ps).pulse_data = pulse_data;
                event(evt).ch(ch).pulse(ps).pulse_data_mV = double(pulse_data) .* (2000/(2^14)) + vrange_bot; % dynamic range 0.1 -> -1.9
                if ~isempty(event(evt).ch(ch).pulse(ps).pulse_data)
%                     event(evt).ch(ch).pulse(ps).baseline_mV = double(event(evt).ch(ch).pulse(ps).baseline) * (2000/(2^14)) + vrange_bot; % 2009-07-15 jjc - now use baseline value from struck
                event(evt).ch(ch).pulse(ps).baseline_mV = double(mean(event(evt).ch(ch).pulse(ps).pulse_data(2:10))) * (2000/(2^14)) + vrange_bot;
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
    %dual_fclose(eventfid);
end
dual_fclose(eventfid);
fclose(indexfid);
%toc


event_struct=event;

%for evt=1:1000
%if isfield(event(evt), 'ch')
%if ~isempty(event(evt).ch)

%for ch=1:length(event(evt).ch)
%if length(event(evt).ch(ch))>0
%disp(evt)
%end
%end
%end
%end
%end
