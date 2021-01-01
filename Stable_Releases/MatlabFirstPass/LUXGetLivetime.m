function [livetime] = LUXGetLivetime(data_path, file_list)
% [livetime] LUXGetLivetime(xml_settings, file_list);
% returns livetime(file_list_number).latch, .end, .file_number
%	latch - is the timestamp of when the acquisition begain
%	end - is the tiemstamp of the last sample that ended the acquisition (multi-event mode)
%	file_number - is the actual file number from the dataset that the timestamps correspond to
%		(file_list_number is the position in the file_list
% 2009-06-24 JJC created

%if ~(xml_settings.global.event_builder_version>=4.3 && xml_settings.global.daq_version>=6.8)
%	livetime=[];
%	return;
%end

%filename_prefix = xml_settings.filename_prefix;
%data_path = xml_settings.data_path;

if nargin<2
	file_list = [];
end

if 1%isempty(file_list) % set file_list to all files. first this just does 1:file containing last event (via index file). it should look through files and make list that way.
%	filename_ind = [data_path,'/',filename_prefix,'/',filename_prefix,'.ind'];
%	indexfid = fopen(filename_ind,'rb','l');
%	fseek(indexfid,double(-7*4),'eof');
%	event_gid = fread(indexfid,4,'uint32');
%	record_format = fread(indexfid,1,'*uint32');
%	nb_files = fread(indexfid,1,'*uint32');
	
%	file_list = 1:nb_files ;
dir_files = dir(data_path);
jj=0;
for ii=1:length(dir_files)
    if dir_files(ii).isdir==0 && strcmp(dir_files(ii).name(end-2:end),'evt') && strcmp(dir_files(ii).name(1:3),'lux')
       jj=jj+1;
       %file_list(jj) = str2num(dir_files(ii).name(22:30));
       files(jj) = dir_files(ii);
    end
end
end
if ~isempty(file_list)
    files = files(file_list);
end
%keyboard
%filename_prefix = xml_settings.filename_prefix;
%data_path = xml_settings.data_path;

for iifile=1:length(file_list)
    options.load_settings=0;
    
    [s fid] = LUXSuperLoader(files(iifile).name,data_path,options);
    fseek(fid,20,'cof');
	%filename = [data_path,'/',filename_prefix,'/',filename_prefix,sprintf('_f%09d',file_list(iifile)),'.evt'];
	%eventfid = fopen(filename,'rb','l');
	%fread(eventfid, 1, 'uint32'); %0
    %fread(eventfid, 1, 'uint32'); %1
    %fread(eventfid, 1, 'uint32'); %2
    %fread(eventfid, 1, 'uint32'); %3
    %fread(eventfid, 1, 'uint32'); %4
    
    

    	nb_seqs_in_file = fread(fid,1,'uint16');
    	for i_seq=1:nb_seqs_in_file
    		livetime(iifile).latch(i_seq) = fread(fid,1,'uint64');
    		livetime(iifile).end(i_seq) = fread(fid,1,'uint64');
    	end
    	livetime(iifile).file_number = file_list(iifile);
       %keyboard
       fclose(fid);

end