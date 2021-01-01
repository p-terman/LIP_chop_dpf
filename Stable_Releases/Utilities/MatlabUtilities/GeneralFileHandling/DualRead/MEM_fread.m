% function [data count_read] = MEM_fread(MEM_fid, count, precision, skip, machineformat);
%
% This function is to be used with files loaded into memory via MEM_fopen, 
% and handles exactly like fread, with the following exceptions:
%
% -only big and little endian machine formats are accepted (if no
% machineformat is entered, or the machine format enteredis not understood,
% the format specified at MEM_fopen is used).
%
% -precision strings should be names of matlab classes (i.e., no int*4,
% rather int32.  non-integer byte formats are not supported.  The
% type1=>type2 and *type forms are supported, as is n*type1 and
% n*type1=>type2, where n is the count of elements to be read between skips
% (only used if skip is non-zero)
%
% -as in fread, count may be Inf, or [m,Inf], but this will be very slow
% for zipped files which are not loaded up to the EOF.
%
% 07/22/06, CED
% Maintained in SVN repository at
% svn+ssh://lxe@darkanalysis.case.edu/Users/xedmcase/Documents/MatlabCodeRepositories/ReleasedCode

function [data count_read] = MEM_fread(MEM_fid, count, precision, skip, machineformat)

%% declare gloabl variables
global gbl_MEM_FID_INFO;

%% define defaults
if nargin<5
    machineformat = '';
end;

if nargin<4
    skip = 0;
end;

if nargin<3
    precision = 'uint8=>char';
end;

if nargin<2
    disp('Warning -- as apposed to fread, in MEM_fread if you don`t specify count, the rest of the currently loaded block is read');
    count = 0;
end;

data = [];
count_read = 0;

%% check input
if ~(isnumeric(MEM_fid) && length(MEM_fid)==1 && MEM_fid==fix(MEM_fid) && MEM_fid>0)
    disp('The first input to MEM_fread should be a single MEM_fid (positive integer)');
    return;
end;

if ~(isnumeric(count) && length(count)<=2 && all(count==fix(count)) && all(count >= 0))
    disp('The second input (count) to MEM_fread should be a positive integer, or a vector of two positive integers, or inf');
    return;
end;

if ~ischar(precision)
    disp('The third argument to MEM_fread should be a precision string.');
    return;
end;

if ~(isnumeric(skip) && length(skip)==1 && skip==fix(skip) && skip >= 0)
    disp('The fourth input (skip) to MEM_fread should be a non-negative integer');
    return;
end;

if ~ischar(machineformat)
    disp('The fifth argument to MEM_fread should be a machine format string.');
    return;
end;

%% check if file is open for reading
if MEM_fid>length(gbl_MEM_FID_INFO) || isempty(gbl_MEM_FID_INFO(MEM_fid).maxblocksize)
    disp('In MEM_ftell:  The MEM_fid entered does not correspond to an open file');
    return;
end;

%% parse precision string
precision_tokens = regexp(precision,'^(\d*)(\*?)(\w+)(=>\w+)?$','tokens');
if length(precision_tokens)<1
    disp('Invalid precision string in MEM_fread -- see fread help, and use matlab precision names');
end;

if ~isempty(precision_tokens{1}{4})
    precision_tokens{1}{4}=precision_tokens{1}{4}(3:end);
end;

skip_interval = 1;
if ~isempty(precision_tokens{1}{1})
    if isempty(precision_tokens{1}{2})
        precision_tokens{1}{3} = [precision_tokens{1}{1} precision_tokens{1}{3}];
    elseif skip==0
        disp('ignoring multiplier in precision string -- only meaningful is skip > 0');
    else
        skip_interval = str2double(precision_tokens{1}{1});
    end;
elseif ~isempty(precision_tokens{1}{2})
    if ~isempty(precision_tokens{1}{4})
        disp('ignoring * at start of precision string, since output precision is explicitly given');
    else
        precision_tokens{1}{4} = precision_tokens{1}{3};
    end;
end;

if isempty(precision_tokens{1}{4})
    precision_tokens{1}{4} = 'double';
end;

if strcmp(precision_tokens{1}{3},'char')
    precision_tokens{1}{3} = 'uint8';
end;

%% interpret dataclass
dataclass = precision_tokens{1}{3};
if strcmp(dataclass,'double') || strcmp(dataclass,'int64') || strcmp(dataclass,'uint64')
    bytes_per_element = 8;
elseif strcmp(dataclass,'single') || strcmp(dataclass,'int32') || strcmp(dataclass,'uint32')
    bytes_per_element = 4;
elseif strcmp(dataclass,'int16') || strcmp(dataclass,'uint16')
    bytes_per_element = 2;
else
    bytes_per_element = 1;
end;

%% parse machine format string
swapbytes = gbl_MEM_FID_INFO(MEM_fid).swapbytes;
if ~isempty(machineformat)
    endiantest = hex2dec('01020304');
    bytewise_endiantest = double(TypeCastAndSwap(uint32(endiantest),'uint8',false,false));
    currently_big_endian = (bytewise_endiantest(1)==1);
    if strcmp(machineformat,'b') || strcmp(machineformat,'ieee-be')
        swapbytes = ~currently_big_endian;
    elseif strcmp(machineformat,'l') || strcmp(machineformat,'ieee-le')
        swapbytes = currently_big_endian;
    else
        disp('machineformat input not understood, using format specified at MEM_fopen');
    end;
end;

%% read data
data = [];
numbytestoread = (prod(count) + skip * fix(prod(count)/skip_interval)) * bytes_per_element;
numbytes_left_to_save = prod(count) * bytes_per_element;
if gbl_MEM_FID_INFO(MEM_fid).filepos==(gbl_MEM_FID_INFO(MEM_fid).offset + gbl_MEM_FID_INFO(MEM_fid).currentblocksize) && ...
        ~isempty(gbl_MEM_FID_INFO(MEM_fid).data) && ...
        strcmp(gbl_MEM_FID_INFO(MEM_fid).merror,'NOT_AT_EOF')
    MEM_freopen(MEM_fid,gbl_MEM_FID_INFO(MEM_fid).filepos);
end;

while true
    if numbytestoread~=0 && numbytestoread~=Inf
        if skip==0
            skipzero_datamask.start = ...
                max(gbl_MEM_FID_INFO(MEM_fid).filepos+1-gbl_MEM_FID_INFO(MEM_fid).offset,1);
            skipzero_datamask.end = ...
                min(gbl_MEM_FID_INFO(MEM_fid).currentblocksize, ...
                gbl_MEM_FID_INFO(MEM_fid).filepos+numbytestoread-gbl_MEM_FID_INFO(MEM_fid).offset);
            numbytes_left_to_save = numbytes_left_to_save - (skipzero_datamask.end+1-skipzero_datamask.start);
        else
            fpos = gbl_MEM_FID_INFO(MEM_fid).offset:(gbl_MEM_FID_INFO(MEM_fid).currentblocksize + gbl_MEM_FID_INFO(MEM_fid).offset - 1);
            datamask = ...
                fpos>=gbl_MEM_FID_INFO(MEM_fid).filepos & ...
                mod(fpos-gbl_MEM_FID_INFO(MEM_fid).filepos,(skip+skip_interval)*bytes_per_element) < (skip_interval*bytes_per_element) & ...
                fpos < (gbl_MEM_FID_INFO(MEM_fid).filepos + numbytestoread);
            numbytes_left_to_save = numbytes_left_to_save - length(find(datamask));
        end;
    else
        if skip==0
            skipzero_datamask.start = ...
                max(gbl_MEM_FID_INFO(MEM_fid).filepos+1-gbl_MEM_FID_INFO(MEM_fid).offset,1);
            skipzero_datamask.end = ...
                gbl_MEM_FID_INFO(MEM_fid).currentblocksize;
        else
            fpos = gbl_MEM_FID_INFO(MEM_fid).offset:(gbl_MEM_FID_INFO(MEM_fid).currentblocksize + gbl_MEM_FID_INFO(MEM_fid).offset - 1);
            datamask = ...
                fpos>=gbl_MEM_FID_INFO(MEM_fid).filepos & ...
                mod(fpos-gbl_MEM_FID_INFO(MEM_fid).filepos,(skip+skip_interval)*bytes_per_element) < (skip_interval*bytes_per_element);
        end;
    end;
    if skip==0
        data = [data ; gbl_MEM_FID_INFO(MEM_fid).data(skipzero_datamask.start:skipzero_datamask.end)];
    else
        data = [data ; gbl_MEM_FID_INFO(MEM_fid).data(datamask)];
    end;

    if numbytes_left_to_save==0 || ...
            isempty(gbl_MEM_FID_INFO(MEM_fid).data) || ...
            ~strcmp(gbl_MEM_FID_INFO(MEM_fid).merror,'NOT_AT_EOF')
        break;
    end;
    MEM_freopen(MEM_fid,gbl_MEM_FID_INFO(MEM_fid).currentblocksize + gbl_MEM_FID_INFO(MEM_fid).offset);
end;

count_read = fix(length(data) / bytes_per_element);
if length(data)~=(count_read*bytes_per_element)
    data = data(1:(count_read*bytes_per_element));
end;

bytes_read = (fix(count_read / skip_interval)*(skip_interval + skip) + mod(count_read,skip_interval)) * bytes_per_element;
bytes_read = min(bytes_read,gbl_MEM_FID_INFO(MEM_fid).offset + gbl_MEM_FID_INFO(MEM_fid).currentblocksize - gbl_MEM_FID_INFO(MEM_fid).filepos);
gbl_MEM_FID_INFO(MEM_fid).filepos = gbl_MEM_FID_INFO(MEM_fid).filepos + bytes_read;

%% typecast data and swap bytes
swapbytes = swapbytes && bytes_per_element > 1;
if ~strcmp(dataclass,'uint8')
    data = TypeCastAndSwap(data,dataclass,false,swapbytes);
end;

%% reshape data (if length(count) == 2)
if length(count)==2
    if count(2)==0 || count(2)==inf
        count(2) = ceil(count_read/count(1));
    end;
    if length(data)<prod(count)
        data = [data;zeros(prod(count)-length(data),1)];
    end;
    data = reshape(data, count);
end;

%% convert data
data=cast(data,precision_tokens{1}{4});

