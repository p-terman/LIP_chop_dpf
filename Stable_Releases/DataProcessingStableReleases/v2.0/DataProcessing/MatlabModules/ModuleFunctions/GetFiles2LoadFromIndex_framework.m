function files_to_load = GetFiles2LoadFromIndex_framework(data_path, event_list)
%function files_to_load = GetFiles2LoadFromIndex_framework(data_path, event_list)
%data_path is the path to the .evt files and .ind file.
%event_list is the list of events that you want to get the file numbers
%for.
%files_to_load is a structure with all the files that contain the events
%you want, with a list showing which events are in which file.
%
%2012-09-19 JJC - Created.
%2012-12-13 JJC - If event_list is one negative number, then load that many
%of the latest events.
%%2012-12-14 JJC - only read complete entries (skip last one if it's not 3
%%elements).

index_file = dir([data_path '/*.ind']);
fidind = fopen([data_path '/' index_file(1).name],'rb');
index_vec = fread(fidind,index_file(1).bytes/8,'double');
index_vec = index_vec(1:(length(index_vec)-mod(length(index_vec),3)));
fclose(fidind);
index_block = reshape(index_vec,3,floor(length(index_vec)/3));
if numel(event_list)==1 && event_list<0
    [temp_e event_block_ind] = intersect(index_block(1,:),index_block(1,(end+event_list+1):end));
else
    [temp_e event_block_ind] = intersect(index_block(1,:),event_list);
end
index_block_trim = index_block(:,event_block_ind);
file_numbers_temp = unique(squeeze(index_block_trim(2,:)));
files_in_path = dir([data_path '/*.evt*']);
for fff=1:length(files_in_path)
    files_in_path_numbers(fff) = str2num(files_in_path(fff).name(22:30));
end
[c ia file_numbers_to_load] = intersect(file_numbers_temp, files_in_path_numbers);
files_to_load = files_in_path(file_numbers_to_load);
for ff=1:length(files_to_load)
    ff_number = str2num(files_to_load(ff).name(22:30));
    [r] = find(ff_number==index_block_trim(2,:));
    files_to_load(ff).evt_list = index_block_trim(1,r);
end

%keyboard
end