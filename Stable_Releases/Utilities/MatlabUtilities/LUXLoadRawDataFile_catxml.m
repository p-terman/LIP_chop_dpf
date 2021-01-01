function [raw livetime] = LUXLoadRawDataFile_catxml( filename_prefix, file_numbers_list, data_path )
%{
%function [raw livetime] = LUXLoadRawDataFile( filename_prefix, file_numbers_list, data_path )
%
%This function loads the raw .dat file into memory. This is useful if your acquisition is in multi-event
%mode and you want to look at everything the Struck saves, ignoring event triggers
%
%inputs
%


%}

if ~exist('filename_prefix', 'var') || isempty(filename_prefix)
	disp('you must enter a filename_prefix');
	return
end

if ~exist('file_numbers_list','var')
	file_numbers_list = []; % default load all files
end

if ~exist('data_path', 'var')
	data_path = [];
end

[unix_status,hostname] = system('hostname');
if strncmp(hostname,'gsk',3) && ~strncmp(hostname(8),'.',1) % Brown Cluster
	flag_gsk10_cluster = 1;
else
	flag_gsk10_cluster = 0;
end

if ~isempty(data_path)
	data_path = [data_path '/' filename_prefix '/'];
end

if isempty(data_path)
	if flag_gsk10_cluster
		data_path = ['/disk/paRAID02/data/ParticleAstro/LUX01/Run009/' filename_prefix '/']; % currently looks for Run009 by default, should search as LUXLoadSettings searches
 	else
		data_path = ['Temporary/LUX01/' filename_prefix '/'];
	end
end

file_list = dir([data_path '*.dat']); % look for data files in data_path
if isempty(file_list)
	fprintf('No .dat files were found on the data path:\n%s', data_path);
	return;
end
tic
for iifile=1:length(file_list)
	file_list_numbers(iifile) = str2double(file_list(iifile).name(22:30));
	iseq_tot=0;
	ps_tot=0;
end
if isempty(file_numbers_list)
	file_numbers_list = file_list_numbers;
end
for file_number = file_numbers_list % loop over all files you want

	file_list_index = find(file_number==file_list_numbers); % find the index in the list of files of the file number that we're looking for
	if isempty(file_list_index)
		fprintf('Could not find a file with file number %d.\n',file_number);
		run_it = 0;
	else
		run_it = 1;
		file_to_load = file_list(file_list_index).name;
	end

	

	% Begin Reading File %
	if run_it
		fid = fopen([data_path file_to_load],'rb','l');
        
        % Read dat xml file
        endianness = fread(fid,1,'ulong');
        xml_length = fread(fid,1,'ulong');
        xml_string = fread(fid,xml_length,'char');
        
        

        
        
        
        % / Read dat xml file
        
		%endianness = fread(fid, 1, 'uint32');
		% read file header
		file_header = fread(fid,8,'uint32');
		nb_chs = file_header(6);
		%if ~exist('nb_seqs','var') % if first file, initialize variables
		%	for ch=1:nb_chs
		%		pulse_data = [];
		%		baseline = [];
		%		timestamp = [];
		%	end
		%end
		
		% read livetime header
		nb_seqs = fread(fid,1,'uint16');
		%livetime_latch = zeros(1,nb_seqs);
		%livetime_end = zeros(1,nb_seqs);
		for iseq=1:nb_seqs
			iseq_tot = iseq_tot + 1;
			livetime(iseq_tot).latch = fread(fid,1,'uint64');
			livetime(iseq_tot).end = fread(fid,1,'uint64');
		end
		
		% loop over channels
		for ch=1:nb_chs
			if file_number == file_numbers_list(1)
				ps_tot(ch) = 0;
			end
			channel_header = fread(fid,4,'uint16');
			nb_pulses = channel_header(2);
			
			% loop over pulses
			for ps=1:nb_pulses
			
				pulse_header = fread(fid,9,'uint32');
				pulse_length = pulse_header(7);
				pulse_baseline = pulse_header(9);
				pulse_timestamp = pulse_header(6) + pulse_header(5)*2^32;
				pulse_samples = fread(fid,pulse_length,'uint16')';
				pulse_samples=pulse_samples(2:end); %kludge to fix first sample from acq
				pulse_length=pulse_length-1;
				
				%{
				pulse_baseline_vec = ones(1,pulse_length)*pulse_baseline;
				pulse_timestamp_vec = ones(1,pulse_length)*(pulse_timestamp-1) + (1:1:pulse_length);
				
				if file_number==file_numbers_list(1) && ps==1 % first file
					raw.ch(ch).pulse_data = pulse_samples; % add new pulse on to end of previous pulses
					raw.ch(ch).baseline = pulse_baseline_vec;
					raw.ch(ch).timestamp = pulse_timestamp_vec;
				else
					raw.ch(ch).pulse_data = [raw.ch(ch).pulse_data pulse_samples]; % add new pulse on to end of previous pulses
					raw.ch(ch).baseline = [raw.ch(ch).baseline pulse_baseline_vec];
					raw.ch(ch).timestamp = [raw.ch(ch).timestamp pulse_timestamp_vec];
				end
				%}
				ps_tot(ch) = ps_tot(ch) + 1;
				raw.ch(ch).pulse(ps_tot(ch)).pulse_data = pulse_samples';
				raw.ch(ch).pulse(ps_tot(ch)).pulse_data_mV = pulse_samples'.*(2000/2^14);
				raw.ch(ch).pulse(ps_tot(ch)).baseline = pulse_baseline;
				raw.ch(ch).pulse(ps_tot(ch)).baseline_mV = pulse_baseline.*(2000/2^14);
				raw.ch(ch).pulse(ps_tot(ch)).timestamp = pulse_timestamp;
				raw.ch(ch).pulse(ps_tot(ch)).length = pulse_length;
			end
			
		end
		
		fclose(fid);
	end
end
toc