function [files_to_load] = GetFiles2Load_framework(data_path, event_list)
%[files_to_load] = GetFiles2Load_framework(data_path, event_list)
% Used by LUXEventLoader to parse a directory of .evt files and find the
% files that have the desired events for loading.
%
% Example: files_to_load = GetFiles2Load_framework('~/data/dir/lux10_20120402T1010/',[1:10 57]);
%
% 2011-08-15 JJC created
% 2012-02-13 JJC gunzip .gz files, and puts them in current directory if no
% write privileges in data_path.

files_to_search = dir([data_path,'/','*.evt*']);
FF=0;
for ff = 1:length(files_to_search)
    temp_loc=0;
    if strcmp(files_to_search(ff).name((end-1):end),'gz')
       %[stat mess] = fileattrib([data_path '/' files_to_search(ff).name]);
       %fprintf('gunzipping file %d...\n',ff);
       %if mess.GroupWrite==1
       %    gunzip([data_path '/' files_to_search(ff).name]);
       %    files_to_search(ff).name = files_to_search(ff).name(1:(end-3));
       %    temp_loc=0;
       %else
           %fprintf('I do not have write privileges for gunzipping in the data_path,\nso I will put the gunzipped file in your current directory.\nI will delete it when I am done.\n');
           copyfile([data_path '/' files_to_search(ff).name], './');
           temp_loc=1;
           gunzip(['./' files_to_search(ff).name]);
           files_to_search(ff).name = files_to_search(ff).name(1:(end-3));
       %end
    end
    if temp_loc==1
        fid = fopen(['./' files_to_search(ff).name],'rb','l');
    else
        fid = fopen([data_path, '/' files_to_search(ff).name],'rb','l');
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
    if ~isempty(events_to_load)
        FF=FF+1;
        if temp_loc==1
            files_to_load(FF).name = [files_to_search(ff).name '.gz'];
        else
            files_to_load(FF).name = files_to_search(ff).name;
        end
        files_to_load(FF).evt_list = events_to_load;
    end
    
    if isempty(event_list)
        fclose(fid);
        if temp_loc==1
            delete(['./' files_to_search(ff).name]);
            delete(['./' files_to_search(ff).name '.gz']);
        end
        return;
    end
    
    fclose(fid);
    if temp_loc==1
       delete(['./' files_to_search(ff).name]);
       delete(['./' files_to_search(ff).name '.gz']);
    end
end
if ~exist('files_to_load','var')
    fprintf('did not find any files with these events...\n');
    files_to_load = [];
end
return;
    
    
    