% function calibrations = LUXLoadPMTCal(filename,path)
%
%   Reads the contents of a simple binary data file for PMT gain calibrations
%
% Inputs:
%       filename - dataset name (e.g. lux01_20090101T1234)
%       path     - [optional] folder where the dataset lies. Do not include filename in it!
%
% Outputs:
%       calibration
%                  .areas_mVns      - [ch events] Contains pulse areas in mVns
%                  .peak_height_mV  - [ch events] Contains pulse height in mV
%                  .peak_time_ns    - [ch events] Contains time of peak height, in 
%
% 20091002 - CHF, created. Using JJC's binary loader code backbone

function calibration = LUXLoadPMTCal(filename,path)

%% set defaults
tic

data = cell(1,1);
names = {{}};
ioerror = 0;

if ~exist('filename','var')
    filename = [];
end;
if ~exist('path','var')
    path = [];
end;


%% check input
if iscell(filename) && length(filename)==1
    filename = filename{1};
end;
if ~ischar(filename)
    disp('requires filename (string or cell containing string) as input argument');
    return;
end;

path = [path '/' filename];

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

%% Checks where this program is running
% this determines the values of some defaults
% Get name of current computer
[unix_status,hostname] = system('hostname');
if strncmp(hostname,'gsk',3) && ~strncmp(hostname(8),'.',1) % Brown cluster
    flag_gsk10_cluster = 1;
else
    flag_gsk10_cluster = 0;
end

% If data-path is not specified, check if the computer is gsk10. If it is,
% load from the common RQ1 location
if ~exist('path', 'var') || isempty(path)
   if flag_gsk10_cluster
        path = sprintf('/Volumes/paRAID02/analysis/matlab/LUX01/RQ1s/%s/',filename);
   else
        path = sprintf('Temporary/LUX01/RQ1s/%s/',filename);
   end 
   else
%    		path = [path '/' filename '/'] ;
   		path = [path '/'] ;
end

% Needs here other path locations to be automatically assigned
% depending on the system


%% Find files
% file_list = dir([path filename '*.rq1.mrg']); % Look for merged files
file_list = dir([path filename '*.pmtcal']); % Look for merged files
if isempty(file_list)
	dis(['There are no merged files on the RQ1s path: ',sprintf('%s',[path filename])])
	return;
end
name = [path file_list(end).name]; % For now, loads the last one - assumes that's the one with the most files

%% open data file, load to memory
fprintf('Loading %s\n',name);
dual_fid = OpenBinaryFile(name,1,100*2^20); % load at most 100 MB at once

if dual_fid.fid<1
    ioerror=-1;
    disp('Failed to open file in ReadNtupleFile');
    return;
end;

%% set some useful variables
endiantest = hex2dec('01020304');
reverse_endiantest = hex2dec('04030201');

%% loop over blocks

for n_block = 1:2 % loop over only 2 blocks, assuming that we're reading a single file w/ 2 blocks (i.e. a merged file)
% read and parse header
    blockhead = dual_fread(dual_fid,1,'uint32');
    if isempty(blockhead) % reached end of file!
        break
    end
    if blockhead == reverse_endiantest
        MEM_SwitchEndianness(dual_fid.fid);
    elseif blockhead ~= endiantest
        dual_fseek(dual_fid,-4,0);
    end

    headerlength=dual_fread(dual_fid,1,'uint16');
    header=char(dual_fread(dual_fid,headerlength,'*uint8')');
    
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
                vartype{nw},filename));
            ioerror = 1;
            dual_fclose(dual_fid);
            return;
        end
        varbytes(nw) = variable_type_bytes(found_vartype);
    end

% read number of lines in this block
    numlines_in_file = dual_fread(dual_fid,1,'int32');

% If the user specified the number of events to load
% then limit the number of lines to load from the block

    numlines = 1:abs(numlines_in_file);

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
    
    entrytype = zeros(1,linewords); %0 = array, 1 = sparse array, 2 = cell array
    for nw = 1:linewords

    found_vartype_RQ = 1;
    
    if ~isempty(found_vartype_RQ) % if this is on the RQ1s_list
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
     begin_read = dual_ftell(dual_fid);
    nl = min(numlines)-1;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
   
    while nl < max_nl
        nl = nl + 1;
        line_offset = varsize_sum*(nl-1);
        %dual_fseek(dual_fid, begin_read, 'bof'); % go to beginning of reading block data
       % if nl==min(numlines)
        	dual_fseek(dual_fid, begin_read+line_offset, 'bof'); % go to data line to read
        %end
        for nw = 1:linewords

        	found_vartype_RQ = 1;

        if ~isempty(found_vartype_RQ) % if this is on the RQ1s_list
            switch entrytype(nw)
                case 0
                    if ~strcmp(vartype{nw},'char')
                        data{n_block}.(names{n_block}{nw})(nl,:) = ...
                            dual_fread(dual_fid,prod(varsize{nw}),vartype{nw});
                    else
                        data{n_block}.(names{n_block}{nw})(nl,:) = ...
                            char(dual_fread(dual_fid,prod(varsize{nw}),'*uint8'));
                    end                        
                case 1
                    thissize = dual_fread(dual_fid,1,'uint16');
                    data{n_block}.(names{n_block}{nw})(nl,:) = ...
                        dual_fread(dual_fid,thissize,vartype{nw});
                case 2
                    thissize = varsize{nw};
                    thissize(thissize<=0) = dual_fread(dual_fid,length(find(thissize<=0)),'uint16');
                    if ~strcmp(vartype{nw},'char')
                        data{n_block}.(names{n_block}{nw}){nl} = ...
                            reshape(dual_fread(dual_fid,prod(thissize),vartype{nw}),thissize);
                    else
                        data{n_block}.(names{n_block}{nw}){nl} = ...
                            char(reshape(dual_fread(dual_fid,prod(thissize),'*uint8'),thissize));
                    end
            end
		else % skip RQ
            dual_fseek(dual_fid, prod(varsize{nw})*varbytes(nw), 'cof');
        end
        end % end of line
        if numlines_in_file < 1
            % check if this is last line
            endcheck = dual_fread(dual_fid,1,'*uint8'); % is there another line after this?
            if isempty(endcheck) % reached end of file!
                if nl < abs(numlines) % need to clean up excess allocated space
                    for nw=1:linewords
                        switch entrytype(nw)
                            case 0
                                thissize = size(data{n_block}.(names{n_block}{nw}));
                                data{n_block}.(names{n_block}{nw}) = ...
                                    reshape(data{n_block}.(names{n_block}{nw})(1:nl,:), ...
                                    [nl thissize(2:end)]);
                            case 1
                                data{n_block}.(names{n_block}{nw}) = ...
                                    data{n_block}.(names{n_block}{nw})(1:nl,:);
                            case 2
                                data{n_block}.(names{n_block}{nw}) = ...
                                    data{n_block}.(names{n_block}{nw})(1:nl);
                        end
                    end
                end
                break
            else % not at end, so go back one byte from checking
                dual_fseek(dual_fid,-1,0);
            end
        end
    end
    % end while loop over lines    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% end loop over blocks
end

dual_fclose(dual_fid);
%% Output calibration structure with the correct dimensionality

calibration.areas_mVns      = data{2}.areas_mVns';
calibration.peak_height_mV  = data{2}.peak_height_mV';
calibration.peak_time_ns    = data{2}.peak_time_ns';

fprintf('Data loaded in %5.1f seconds\n',toc)
