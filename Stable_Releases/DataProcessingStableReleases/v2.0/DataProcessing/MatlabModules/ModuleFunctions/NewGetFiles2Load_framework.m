function [files_to_load] = NewGetFiles2Load_framework(data_path, event_list)
%[files_to_load] = NewGetFiles2Load_framework(data_path, event_list)
% Used by LUXEventLoader to parse a directory of .evt files and find the
% files that have the desired events for loading.
%
% Example: files_to_load = NewGetFiles2Load_framework('~/data/dir/lux10_20120402T1010/',[1:10 57]);
%
% 2011-08-15 JJC created
% 2012-02-13 JJC gunzip .gz files, and puts them in current directory if no
% write privileges in data_path.

files_to_search = dir([data_path,'/','lux*.evt*']);
FF=0;

for hh=1:length(files_to_search)
    file_numbers(hh) = str2num(files_to_search(hh).name(22:30));
    files_to_search(hh).file_number = str2num(files_to_search(hh).name(22:30));
end

if any(diff(file_numbers)>1) && length(files_to_search)>100
    [files_to_load] = GetFiles2Load_framework(data_path, event_list);
    return;
end


%keyboard
% open last file - what is the event range? guess number per file, and max
% number of events
if strcmp(files_to_search(end).name((end-1):end), 'gz')
    copyfile([data_path '/' files_to_search(end).name], './');
    temp_loc=1;
    gunzip(['./' files_to_search(end).name]);
    fid = fopen(['./' files_to_search(end).name(1:(end-3))], 'rb','l');
else
    fid = fopen([data_path, '/', files_to_search(end).name], 'rb','l');
    temp_loc=0;
end
endendian_check = fread(fid,1,'uint32');
settings_length = fread(fid,1,'uint32');
fseek(fid,settings_length,'cof'); % skip settings
file_header = fread(fid,4,'uint32');
nb_evts_in_file = file_header(4);
evt_index_block = fread(fid,2*nb_evts_in_file,'uint32');
events_in_file_numbers = evt_index_block(logical(mod(1:(2*nb_evts_in_file),2)));
last_event_in_set = events_in_file_numbers(end);
events_in_file_positions = evt_index_block(~logical(mod(1:(2*nb_evts_in_file),2)));
guess_nb_evts_per_file = length(events_in_file_numbers); % round(events_in_file_numbers(end) / length(files_to_search));
events_per_file(1) = length(events_in_file_numbers);
% guess which file has the first event
gg = ceil(event_list(1) / guess_nb_evts_per_file);
fclose(fid);
if temp_loc==1
    delete(['./' files_to_search(end).name]);
    %delete(['./' files_to_search(end).name '.gz']);
    delete(['./' files_to_search(end).name(1:(end-3))]);
end
gglt=[];
gggt=[];
files_to_skip=[];
increment_iterator=0;
while ~isempty(event_list) && ~isempty(files_to_search)
    %keyboard
    %fprintf('event list = [%d:%d]\n',event_list(1), event_list(end));
    gg=min(gg,length(files_to_search));
    if gg==0; gg=1; end
    
    % load the guess file
    if strcmp(files_to_search(gg).name((end-1):end), 'gz')
        copyfile([data_path '/' files_to_search(gg).name], './');
        temp_loc=1;
        gunzip(['./' files_to_search(gg).name]);
        fid = fopen(['./' files_to_search(gg).name(1:(end-3))], 'rb','l');
    else
        fid = fopen([data_path, '/', files_to_search(gg).name], 'rb','l');
        temp_loc=0;
    end
    endian_check = fread(fid,1,'uint32');
    settings_length = fread(fid,1,'uint32');
    fseek(fid,settings_length,'cof'); % skip settings
    file_header = fread(fid,4,'uint32');
    nb_evts_in_file = file_header(4);
    evt_index_block = fread(fid,2*nb_evts_in_file,'uint32');
    events_in_file_numbers = evt_index_block(logical(mod(1:(2*nb_evts_in_file),2)));
    events_in_file_positions = evt_index_block(~logical(mod(1:(2*nb_evts_in_file),2)));
    events_to_load = intersect(event_list,events_in_file_numbers);
    event_list = setdiff(event_list,events_to_load); % remove the events already found from the list
    events_per_file(end+1) = length(events_in_file_numbers);
    fclose(fid);
    if temp_loc==1
        delete(['./' files_to_search(gg).name]);
        %delete(['./' files_to_search(gg).name '.gz']);
        delete(['./' files_to_search(gg).name(1:(end-3))]);
    end
    if ~isempty(events_to_load)
        FF=FF+1;
        files_to_load(FF).name = files_to_search(gg).name;
        files_to_load(FF).evt_list = events_to_load;
    end
    %files_to_skip(end+1) = gg;
    %fprintf('gg = %d, %s (%d:%d)\n',gg, files_to_search(gg).name, events_in_file_numbers(1), events_in_file_numbers(end));
    
    % if guess was completely off (no desired events in the guess file)
    if isempty(events_to_load)
        if increment_iterator<5 || length(files_to_search)<10
            increment_iderator = increment_iterator+1;
            if ~isempty(events_in_file_numbers) && event_list(1)<events_in_file_numbers(1)
                
                files_to_search = files_to_search(setdiff(1:length(files_to_search), gg));
                gg=gg-1;
                
            else % event_list(1)>events_in_file_numbers(end)
                
                files_to_search = files_to_search(setdiff(1:length(files_to_search), gg));
                gg=gg;
                
            end
        else
            % incremented enough, let's make a new guess
            increment_iterator=0;
            gg_start = ceil(event_list(1) / mean(events_per_file));
            [unused gg] = min(abs([files_to_search(:).file_number]-gg_start));
            %keyboard
        end
    end
    
    
    % if in the middle of a range of desired events, then keep going back to previous file until out of range
    ggg=gg;
    while ~isempty(events_to_load) && ~isempty(event_list) && event_list(1)<=events_in_file_numbers(1)
        ggg=ggg-1;
        ggg=max(ggg,1);
        
        % load the guess file
        if strcmp(files_to_search(ggg).name((end-1):end), 'gz')
            copyfile([data_path '/' files_to_search(ggg).name], './');
            temp_loc=1;
            gunzip(['./' files_to_search(ggg).name]);
            fid = fopen(['./' files_to_search(ggg).name(1:(end-3))], 'rb','l');
        else
            fid = fopen([data_path, '/', files_to_search(ggg).name], 'rb','l');
            temp_loc=0;
        end
        endian_check = fread(fid,1,'uint32');
        settings_length = fread(fid,1,'uint32');
        fseek(fid,settings_length,'cof'); % skip settings
        file_header = fread(fid,4,'uint32');
        nb_evts_in_file = file_header(4);
        evt_index_block = fread(fid,2*nb_evts_in_file,'uint32');
        events_in_file_numbers = evt_index_block(logical(mod(1:(2*nb_evts_in_file),2)));
        events_in_file_positions = evt_index_block(~logical(mod(1:(2*nb_evts_in_file),2)));
        events_to_load = intersect(event_list,events_in_file_numbers);
        event_list = setdiff(event_list,events_to_load); % remove the events already found from the list
        %fprintf('(%d) ggg=ggg-1 = %d, %s (%d:%d)\n', event_list(1) ,ggg, files_to_search(ggg).name, events_in_file_numbers(1), events_in_file_numbers(end));        %load previous file, look for data
        events_per_file(end+1) = length(events_in_file_numbers);
        if ~isempty(events_to_load)
            FF=FF+1;
            files_to_load(FF).name = files_to_search(ggg).name;
            files_to_load(FF).evt_list = events_to_load;
        end
        if isempty(event_list)
            fclose(fid);
            if temp_loc==1
                delete(['./' files_to_search(ggg).name]);
                %delete(['./' files_to_search(ggg).name '.gz']);
                delete(['./' files_to_search(ggg).name(1:(end-3))]);
            end
            return;
        end
        
        fclose(fid);
        if temp_loc==1
            delete(['./' files_to_search(ggg).name]);
            %delete(['./' files_to_search(ggg).name '.gz']);
            delete(['./' files_to_search(ggg).name(1:(end-3))]);
        end
        files_to_search = files_to_search(setdiff(1:length(files_to_search), ggg));
    end
    
    % if in the middle of a range of desired events, then keep going to next file until out of range
    while ~isempty(events_to_load) && ~isempty(event_list) && event_list(1)>events_in_file_numbers(end)
        ggg=ggg;%+1;
        ggg=min(ggg,length(files_to_search));
        
        %load previous file, look for data
        % load the guess file
        if strcmp(files_to_search(ggg).name((end-1):end), 'gz')
            copyfile([data_path '/' files_to_search(ggg).name], './');
            temp_loc=1;
            gunzip(['./' files_to_search(ggg).name]);
            fid = fopen(['./' files_to_search(ggg).name(1:(end-3))], 'rb','l');
        else
            fid = fopen([data_path, '/', files_to_search(ggg).name], 'rb','l');
            temp_loc=0;
        end
        endian_check = fread(fid,1,'uint32');
        settings_length = fread(fid,1,'uint32');
        fseek(fid,settings_length,'cof'); % skip settings
        file_header = fread(fid,4,'uint32');
        nb_evts_in_file = file_header(4);
        evt_index_block = fread(fid,2*nb_evts_in_file,'uint32');
        events_in_file_numbers = evt_index_block(logical(mod(1:(2*nb_evts_in_file),2)));
        events_in_file_positions = evt_index_block(~logical(mod(1:(2*nb_evts_in_file),2)));
        events_to_load = intersect(event_list,events_in_file_numbers);
        event_list = setdiff(event_list,events_to_load); % remove the events already found from the list
        events_per_file(end+1) = length(events_in_file_numbers);
        %fprintf('(%d) ggg=ggg+1 = %d, %s (%d:%d)\n',event_list(1), ggg, files_to_search(ggg).name, events_in_file_numbers(1), events_in_file_numbers(end));
        if ~isempty(events_to_load)
            FF=FF+1;
            files_to_load(FF).name = files_to_search(ggg).name;
            files_to_load(FF).evt_list = events_to_load;
        end
        if isempty(event_list)
            fclose(fid);
            if temp_loc==1
                delete(['./' files_to_search(ggg).name]);
                %delete(['./' files_to_search(ggg).name '.gz']);
                delete(['./' files_to_search(ggg).name(1:(end-3))]);
            end
            return;
        end
        
        fclose(fid);
        if temp_loc==1
            delete(['./' files_to_search(ggg).name]);
            %delete(['./' files_to_search(ggg).name '.gz']);
            delete(['./' files_to_search(ggg).name(1:(end-3))]);
        end
        files_to_search = files_to_search(setdiff(1:length(files_to_search), gg));
    end
    
    %if isempty(event_list)
    %    fclose(fid);
    %    if temp_loc==1
    %        delete(['./' files_to_search(gg).name]);
    %        %delete(['./' files_to_search(gg).name '.gz']);
    %        delete(['./' files_to_search(gg).name(1:(end-3))]);
    %   end
    %  return;
    %end
    
    %     fclose(fid);
    %     if temp_loc==1
    %         delete(['./' files_to_search(gg).name]);
    %         delete(['./' files_to_search(gg).name '.gz']);
    %     end
    %
    
end






