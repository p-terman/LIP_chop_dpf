% function dp=LUX01LoadRQ1s(filename_prefix,RQ1s_path,evt_list,RQ1s_list,filenumber)
%
%   Reads the contents of a simple binary data file.
%
% Inputs:
%       filename_prefix - dataset name (e.g. lux01_20090101T1234)
%       evt_list        - event list to load. Can be a number, or an array to load (Optional)
%       RQ1s_list       - cell array with list of RQs to load (Optional)
%       RQ1s_path       - path where the RQ1 files are located (Optional)
%
% Outputs:
%       dp              - 1st pass RQs, structure (i.e. dp{1}.fieldname)
%
%   2008-05-12 CED
%   2009-04-23 LdV Renamed it from ReadNtupleFile to LUX01LoadRQ1s; minor
%                  modifications, limits the number of events read
%   2009-04-24 ChF Using filename_prefix and RQ1s_path separately. Check if
%                  system is gsk-10 to use default path. Outputs dp{1} with
%                  correct dimensions. Looks for .mrg files in folder.
%   2009-04-27 JJC Changed nb_evts to evt_list so that it can be an array not starting at 1.
%   2009-04-27 JJC Implemented RQ1s_list input. Also with evt_list starting above event 1, remove empty events.
%   2010-02-04 JJC Fixed endianness so it always reads accurately
%   regardless of what machine you are on, or the file was written on.
%   2010-05-27 JJC v5, fopen_endian_safe now reads the settings as well and
%   so you start at the beginning of the block header string size, as
%   normal.
function dp=LUXLoadRQ1s(filename_prefix,RQ1s_path,evt_list,RQ1s_list,filenumber)
tic
%% set defaults
dp =[];
data = cell(1,1);
names = {{}};
ioerror = 0;
if ~exist('filename_prefix','var')
    filename_prefix = [];
end;
if ~exist('evt_list','var')
    evt_list = [];
end;
if ~exist('RQ1s_path','var')
    RQ1s_path = [];
end;
if ~exist('RQ1s_list','var')
    RQ1s_list = [];
end;

if ~exist('filenumber','var')
	filenumber = [];
end

if isempty(filename_prefix)
    fprintf('Looking for the latest dataset in %s ...\n', RQ1s_path);
    latest_filename_list = {};
    latest_filename_n = [];
    
    % Making assumption here -- all files / datasets are to start with 'lux' prefix
    % Note -- inclusive with LUX 0.1 files, so this is OK
    filename_list = dir([RQ1s_path,filesep,'lux*']);
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
    
    if isempty(filename_list) % no results found or left after name checking
        fprintf('There are no datasets in %s\n', RQ1s_path);
        return;
    else
        dates = [filename_list.datenum];
        [tmp, ii_latest] = max(dates);
        filename_prefix = filename_list(ii_latest).name;
        fprintf('Loading RQ1s from dataset %s\n', filename_prefix);
    end
end



%% check input
if iscell(filename_prefix) && length(filename_prefix)==1
    filename_prefix = filename_prefix{1};
end;
%if ~ischar(filename_prefix)
%    disp('requires filename_prefix (string or cell containing string) as input argument');
%    return;
%end;

%% set reference constants
variable_type_strings={ ...
    'char', ...
    'int8', ...
    'int16', ...
    'int32', ...
    'int64', ...
    'uint8', ...
    'uint16', ...
    'uint32', ...
    'uint64', ...
    'single', ...
    'double', ...
    'float' ...
    'float128' ...
    'logical' ...
    };

variable_type_bytes=[1 1 2 4 8 1 2 4 8 4 8 4 16 1];


%% Find files
%RQ1s_path = [RQ1s_path filesep filename_prefix];
if isempty(filenumber)
  disp(RQ1s_path);
	file_list = dir([RQ1s_path filesep '*.rq1.mrg']); % Look for merged files
    if isempty(file_list)
        fprintf('Could not find any merged files. Looking for filenumber 1...\n');
        filenumber = 1;
        file_list = dir([RQ1s_path filesep sprintf('*_f%09d.rq1',filenumber)]);
    end
else
	file_list = dir([RQ1s_path filesep sprintf('*_f%09d.rq1',filenumber)]);
end
if isempty(file_list)
	fprintf('Could not find file\n');
	return;
end
filename = [RQ1s_path filesep file_list(end).name]; % For now, loads the last one - assumes that's the one with the most files

% file_list = dir([RQ1s_path filename_prefix '*.rq1']); % Look for merged files
% if isempty(file_list)
%     dis(['There are no merged files on the RQ1s path: ',sprintf('%s',[RQ1s_path filename_prefix])])
%     return;
% end
% filename = [RQ1s_path file_list(2).name]; % For now, loads the last one - assumes that's the one with the most files

% %% Find files
% file_list = dir([RQ1s_path '*.rq1.mrg']); % Look for merged files
% if isempty(file_list)
% 	dis(['There are no merged files on the RQ1s path: ',sprintf('%s',RQ1s_path)])
% 	return;
% end
% filename = [RQ1s_path file_list(end).name]; % For now, loads the last one - assumes that's the one with the most files

% file_list = dir([RQ1s_path filename_prefix '*.rq1']); % Look for merged files
% if isempty(file_list)
%     dis(['There are no merged files on the RQ1s path: ',sprintf('%s',[RQ1s_path filename_prefix])])
%     return;
% end
% filename = [RQ1s_path file_list(2).name]; % For now, loads the last one - assumes that's the one with the most files

%% open data file, load to memory
%filename = '~/RQ1s/lux01_20090729T2216/lux01_20090729T2216_f000000001.rq1.tmp2';
fprintf('Loading %s\n',filename);
%keyboard
[fid switch_endian] = fopen_endian_safe(filename); % load at most 100 MB at once
% if switch_endian
%    [fn,pm,mf]=fopen(1);
%    if strcmp(mf,'ieee-le')
%       fids = fopen(filename,'rb','b');
%       filename = [filename '.se'];
%       fidw = fopen(filename,'wb','l');
%       fseek(fids,0,'eof');
%       swfilesize = ftell(fids);
%       fseek(fids,0,'bof');
%       swfile = fread(fids,swfilesize,'int8');
%       fwrite(fidw,swfile,'int8',0,'l');
%       fprintf('switched endianness, now reading filename %s\n',filename);
%    end
%    
%    if strcmp(mf,'ieee-be')
%       fids = fopen(filename,'rb','l');
%       filename = [filename '.se'];
%       fidw = fopen(filename,'rb','b');
%       fseek(fids,0,'eof');
%       swfilesize = ftell(fids);
%       fseek(fids,0,'bof');
%       swfile = fread(fids,swfilesize,'int8');
%       fwrite(fidw,swfile,'int8');
%       fprintf('switched endianness, now reading filename %s\n',filename);
%    end
%    fclose(fidw);
%    fclose(fids);
%    [status, result] = unix(['chmod 755 ' filename]);
%     if status ~= 0
%         disp(['Warning: unable to set file permissions for ' filename]);
%         dis('Message: %s', result);
%     end
%    fid = fopen_endian_safe(filename);
% end

if fid<1
    ioerror=-1;
    disp('Failed to open file');
    return;
end;
%endiantest = hex2dec('01020304');
%reverse_endiantest = hex2dec('04030201');

%endianness = fread(fid,1,'uint32')
%if isempty(endianness) % reached end of file!
%    return;
%end
%if endianness == reverse_endiantest
%    fprintf('switching endianness\n');
%    switch_endian = 1;
%    MEM_SwitchEndianness(fid.fid);
%elseif endianness ~= endiantest
%    fseek(fid,-4,0);
%    fprintf('endian test inconclusive\n');
%elseif endianness == endiantest
%    fprintf('endianness good\n');
%    switch_endian = 0;
%end

%% set some useful variables


%% loop over blocks

for n_block = 1:2 % loop over only 2 blocks, assuming that we're reading a single file w/ 2 blocks (i.e. a merged file)
    %keyboard
    % read and parse header
    if n_block == 2 % skip to header of second block
        [fid switch_endian] = fopen_endian_safe(filename); % load at most 100 MB at once
        %if switch_endian
        %    MEM_SwitchEndianness(fid.fid);
        %    fprintf('switching endianness\n');
        %end
        %keyboard
        fseek(fid, block2_offset, 'bof');
    end

    
    
    headerlength=fread(fid,1,'uint16');
    header=char(fread(fid,headerlength,'*uint8')');
    
    delims=[0 regexp(header,'[\.:; ]')];

    linewords=(length(delims)-1)/3;

    names{n_block}=cell(1,linewords);
    vartype=cell(1,linewords);
    varsize=cell(1,linewords);
    varbytes=zeros(1,linewords);
    for nw=1:linewords
        names{n_block}{nw}=header((delims(3*nw-2)+1):(delims(3*nw-1)-1));
        vartype{nw}=header((delims(3*nw-1)+1):(delims(3*nw)-1));
        varsize{nw}=sscanf(header((delims(3*nw)+1):(delims(3*nw+1)-1)),'%d,')';
    end

    % determine variable types
    for nw=1:linewords
        found_vartype = find(strcmp(variable_type_strings,vartype{nw}),1,'first');
        if isempty(found_vartype)
            disp(sprintf( ...
                'did not understand vartype %s in header, while reading file %s in ReadNtupleFile', ...
                vartype{nw},filename_prefix));
            ioerror = 1;
            fclose(fid);
            return;
        end
        varbytes(nw) = variable_type_bytes(found_vartype);
    end

    % read number of lines in this block
    numlines_in_file = fread(fid,1,'int32');

    % If the user specified the number of events to load
    % then limit the number of lines to load from the block
    if ~isempty(evt_list)
        if n_block > 1
            numlines = min(evt_list):min(max(evt_list),numlines_in_file);
        else
            numlines = 1:abs(numlines_in_file);
        end
    else
        numlines = 1:abs(numlines_in_file);
    end

    numlines_global = numlines ; % store numlines requested total

    chunks = ceil(length(numlines_global) / 100e3) ;



    for chk=1:chunks

        starter = (chk-1)*100e3+1;
        ender = (chk-1)*100e3+100e3;
        if ender>length(numlines_global)
            ender = length(numlines_global);
        end
        numlines = numlines_global(starter:ender) ;


        % pre-allocate space for this block
        % There is possibly room for improvement here, definitely room for
        % change as needed.
        %
        % Currently, allocation is based on abs(numlines), so if numlines==0
        % there is no preallocation.
        %
        % If there are no zero or negative sizes for a block, the full size is
        % preallocated (going by numlines), and any unused lines are removed
        % after the block is read.
        %
        % If the entry is a variable length vector, a sparse matrix is
        % used to store the data
        %
        % Otherwise a cell array with one index (line number) is used to store
        % the data -- the contents of each cell are an array with indices
        % matching the listing in the header



        found_vartype_RQ = zeros(1,linewords);
        entrytype = zeros(1,linewords); %0 = array, 1 = sparse array, 2 = cell array
        for nw = 1:linewords
            if (~isempty(RQ1s_list)) && (n_block == 2)
                found_vartype_RQ_holder = find(strcmp(RQ1s_list,names{n_block}{nw}),1,'first');
                if isempty(found_vartype_RQ_holder)
                    found_vartype_RQ(nw) = 0;
                else
                    found_vartype_RQ(nw) = 1;
                end
            else
                found_vartype_RQ(nw) = 1;
            end
            if found_vartype_RQ(nw)~=0 % if this is on the RQ1s_list
                if all(varsize{nw}>0)
                    if ~strcmp(vartype{nw},'char')
                        data{n_block}.(names{n_block}{nw}) = zeros([abs(length(numlines)) varsize{nw}]);
                    else
                        data{n_block}.(names{n_block}{nw}) = repmat(' ',[abs(length(numlines)) varsize{nw}]);
                    end
                elseif length(varsize{nw})==1 && varsize{nw}<=0 && ~strcmp(vartype{nw},'char')
                    data{n_block}.(names{n_block}{nw}) = sparse(abs(length(numlines)), abs(varsize{nw}));
                    entrytype(nw) = 1;
                else
                    data{n_block}.(names{n_block}{nw}) = cell(1,abs(length(numlines)));
                    entrytype(nw) = 2;
                end
            end
        end

        % read data
        if (numlines_in_file < 1) && (isempty(evt_list)) % if not using event list and number of lines is undefined
            max_nl = inf;
        else
            max_nl = max(numlines);
        end

        varsize_sum = 0;
        for ii=1:length(varsize)
            varsize_sum = varsize_sum + prod(varsize{ii})*varbytes(ii);
        end

        nl = min(numlines)-1;
        if chk==1
            begin_read = ftell(fid);
            block2_offset = length(numlines)*varsize_sum + begin_read;  % amount to skip when moving to header of next block
            fclose(fid);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        
        name = names{n_block};
        data_mex = data{n_block};
        %keyboard
        if n_block==2
        LUXLoadRQ1s_mex ;
        %keyboard
        
            data{n_block} = data_mex ;
        end
            %  if n_block == 2
        % keyboard
        %end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % end loop over blocks


        %% finish up
        if n_block == 2
            if isempty(evt_list)
                evt_list = numlines_global;
            end

            for ii=1:linewords
                if (length(varsize{ii})==2) && (found_vartype_RQ(nw)~=0)
                    data{2}.(names{2}{ii}) = reshape(data{2}.(names{2}{ii}),varsize{ii}(1),4,length(numlines));

                    %		data_temp = data{2};
                    %		data{2}.(names{2}{ii}) = zeros(varsize{ii},numlines(end));
                    %		chs = varsize{ii}(2);
                    %		pulses = varsize{ii}(1);
                    %		for iii=1:varsize{ii}(2)
                    %			data{2}.(names{2}{ii})(1:pulses,iii,:) = data_temp.(names{2}{ii})((pulses*iii):(pulses*iii+chs),:)
                    %		end
                end
            end

            if 1
                %% Output dp{1} with the correct dimensionality
                fn = fieldnames(data{2}); %Get field names from data{2}
                for ii = 1:numel(fn)
                    data{2}.(char(fn{ii})) = squeeze(data{2}.(char(fn{ii})));
                    vals = data{2}.(char(fn{ii}));
                    sizes = size(vals);
                    if length(size(data{2}.(char(fn{ii})))) == 2 %transpose if 2-dimensional
                        dp_tmp{1}.(char(fn{ii})) = vals;
                    elseif length(size(data{2}.(char(fn{ii})))) == 3 %permute dimensions [1 2 3]->[2 3 1] if 3-dimensional
                        dp_tmp{1}.(char(fn{ii})) = squeeze(permute(vals,[1 2 3]));
                    end
                    if sizes(2) == 1
                        dp_tmp{1}.(char(fn{ii})) = vals';
                    end
                    %dp{1}.(char(fn{ii})) = dp{1}.(char(fn{ii}))(:,evt_list);
                end
            else
                dp_tmp{1} = data{2} ;


            end
            if chk==1
                dp = dp_tmp;
                datah{1} = data{1} ;
            else
                dp = LUX01Concatdp(dp,dp_tmp) ; %concat dp
            end
            clear data;
        end %chunking for loop
    end
end
% Add set, header
dp{1}.set = filename_prefix;
dp{1}.evt_list = evt_list;
if exist('datah','var')
    dp{1}.RQ1_header = datah{1};
end

toc
end


function [fid switch_endian settings] = fopen_endian_safe(filename)

[fn,pm,mf]=fopen(1) ; % mf is endianness of this machine

fid = fopen(filename, 'rb', 'n'); % open in big endian (PowerPC)

% endianness check
endian_check = fread(fid,1,'uint32');
if ~strcmp(dec2hex(endian_check),'1020304')
    %fprintf('endianness is wrong...\n');
    switch_endian = 1;
    fclose(fid);
    if strcmp(mf(1:7), 'ieee-le')
      %  fprintf('reading file as big endian\n');
        fid = fopen(filename, 'rb', 'b');
        endian_check = fread(fid,1,'uint32');
    else
     %   fprintf('reading file as little endian\n');
        fid = fopen(filename, 'rb', 'l');
        endian_check = fread(fid,1,'uint32');
    end
    
else
    switch_endian = 0;
    %fprintf('endianness is consistent (%s)\n', mf);
end
    xmlstringlength = fread(fid,1,'uint32');
    xmlstring = fread(fid,xmlstringlength,'int8');
    %settings = XMLParser(char(xmlstring'));

end


