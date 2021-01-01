% function [data,names,ioerror]=ReadNtupleFile(filename, skiplines, readlines)
%
%   Reads the contents of a simple binary data file.
%
% Inputs:
%       filename        - slow control data file name with path.
%
% Outputs:
%       names           - 2D cell array with variable names
%       data          - structure with data values.  The syntax is:
%                           data{nb}.(names{nb,nv})(nl,:)
%                           where nb indicates the block, nv the variable, 
%                           and nl the dataline
%       ioerror         - ~= 0 if problem.
%
%   05/12/08, CED

function [data,names,ioerror]=ReadNtupleFile(filename)

%% set defaults
data = cell(1,1);
names = {{}};
ioerror = 0;

if nargin<1
    filename = [];
end;

%% check input
if iscell(filename) && length(filename)==1
    filename = filename{1};
end;

if ~ischar(filename)
    disp('requires filename (string or cell containing string) as input argument');
    return;
end;

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
    'float128' ...
    };

variable_type_bytes=[1 1 2 4 8 1 2 4 8 4 8 16];

%% open data file, load to memory
dual_fid = OpenBinaryFile(filename,1,100*2^20); % load at most 100 MB at once

if dual_fid.fid<1
    ioerror=-1;
    disp('Failed to open file in ReadNtupleFile');
    return;
end;

%% set some useful variables
endiantest = hex2dec('01020304');
reverse_endiantest = hex2dec('04030201');

%% loop over blocks
n_block = 0;

while true % while loop continues until end of file is reached
    n_block = n_block + 1;
%% read and parse header
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
    
    delims=[0 regexp(header,'[\.,:; ]')];

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

%% determine variable types
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

%% read number of lines in this block
    numlines = dual_fread(dual_fid,1,'int16');
    
%% pre-allocate space for this block
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
        if all(varsize{nw}>0)
            if ~strcmp(vartype{nw},'char')
                data{n_block}.(names{n_block}{nw}) = zeros([abs(numlines) varsize{nw}]);
            else
                data{n_block}.(names{n_block}{nw}) = repmat(' ',[abs(numlines) varsize{nw}]);
            end
        elseif length(varsize{nw})==1 && varsize{nw}<=0 && ~strcmp(vartype{nw},'char')
            data{n_block}.(names{n_block}{nw}) = sparse(abs(numlines), abs(varsize{nw}));
            entrytype(nw) = 1;
        else
            data{n_block}.(names{n_block}{nw}) = cell(1,abs(numlines));
            entrytype(nw) = 2;
        end
    end

%% read data
    if numlines < 1
            max_nl = inf;
        else
            max_nl = numlines;
    end

    nl = 0;
    while nl < max_nl
        nl = nl + 1;
        for nw = 1:linewords
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
        end
        if numlines < 1
            % check if this is last line
            endcheck = dual_fread(dual_fid,1,'*uint8');
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
            else
                dual_fseek(dual_fid,-1,0);
            end
        end
    end
        

%% end loop over blocks
end

%% finish up
dual_fclose(dual_fid);


