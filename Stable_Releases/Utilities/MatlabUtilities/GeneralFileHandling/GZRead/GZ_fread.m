% function [data count_read] = GZ_fread(GZ_fid, count, precision, skip, machineformat);
%
% This function interfaces with the GZRead function, and handles exactly
% like fread, with the following exceptions:
%
% -only big and little endian machine formats are accepted (if no
% machineformat is entered, or the machine format enteredis not understood,
% the format specified at GZ_fopen is used).
%
% -precision strings should be names of matlab classes (i.e., no int*4,
% rather int32.  non-integer byte formats are not supported.  The
% type1=>type2 and *type forms are supported, as is n*type1 and
% n*type1=>type2, where n is the count of elements to be read between skips
% (only used if skip is non-zero)
%
% -as in fread, count may be Inf, or [m,Inf], but this will be very slow
% for zipped files.  This is because in GZRead, the total space is
% allocated before reading begins, which requires finding the filesize
% before calling GZRead.  The filesize is found be running GZ_filesize,
% which effectively unzips the entire file.
%
% 06/13/06, CED
% Release Version 1.0

function [data count_read] = GZ_fread(GZ_fid, count, precision, skip, machineformat);

%% declare gloabl variables
global gbl_GZ_FID_INFO;

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
    disp('Warning -- as apposed to fread, in GZ_fread if you don`t specify count, a single character is returned.');
    count = 1;
end;

data = [];
count_read = 0;

%% check input
if ~(isnumeric(GZ_fid) && length(GZ_fid)==1 && GZ_fid==fix(GZ_fid) && GZ_fid>0)
    disp('The first input to GZ_fread should be a single GZ_fid (positive integer)');
    return;
end;

if ~(isnumeric(count) && length(count)<=2 && all(count==fix(count)) && all(count >= 0))
    disp('The second input (count) to GZ_fread should be a positive integer, or a vector of two positive integers, or inf');
    return;
end;

if ~ischar(precision)
    disp('The third argument to GZ_fread should be a precision string.');
    return;
end;

if ~(isnumeric(skip) && length(skip)==1 && skip==fix(skip) && skip >= 0)
    disp('The fourth input (skip) to GZ_fread should be a non-negative integer');
    return;
end;

if ~ischar(machineformat)
    disp('The fifth argument to GZ_fread should be a machine format string.');
    return;
end;

%% check if file is open for reading
if GZ_fid>length(gbl_GZ_FID_INFO) || isempty(gbl_GZ_FID_INFO(GZ_fid).gzfid) || ~strcmp(gbl_GZ_FID_INFO(GZ_fid).mode,'r')
    disp('In GZ_ftell:  The GZ_fid entered does not correspond to a file open for reading');
    return;
end;

%% parse precision string
precision_tokens = regexp(precision,'^(\d*)(\*?)(\w+)(=>\w+)?$','tokens');
if length(precision_tokens)<1
    disp('Invalid precision string in GZ_fread -- see fread help, and use matlab precision names');
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
        skip_interval = str2num(precision_tokens{1}{1});
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

%% parse count
if count(end)==Inf
    filesize = GZ_filesize(GZ_fopen(GZ_fid));
    count_to_read = fix((filesize-GZ_ftell(GZ_fid))/bytes_per_element);
    count_to_read = skip_interval*fix(count_to_read/(skip+skip_interval)) + ...
        min(skip_interval,mod(count_to_read,skip+skip_interval));
    if length(count)==1
        count = count_to_read;
    else
        count(2) = ceil(count_to_read/count(1));
    end;
end;

%% parse machine format string
swapbytes = gbl_GZ_FID_INFO(GZ_fid).swapbytes;
if ~isempty(machineformat)
    endiantest = hex2dec('01020304');
    bytewise_endiantest = double(TypeCastAndSwap(uint32(endiantest),'uint8',false,false));
    currently_big_endian = (bytewise_endiantest(1)==1);
    if strcmp(machineformat,'b') || strcmp(machineformat,'ieee-be')
        swapbytes = ~currently_big_endian;
    elseif strcmp(machineformat,'l') || strcmp(machineformat,'ieee-le')
        swapbytes = currently_big_endian;
    else
        disp('machineformat input not understood, using format specified at GZ_fopen');
    end;
end;

%% read data
numtoread = prod(count) + skip * fix(prod(count)/skip_interval);
startpos = double(gbl_GZ_FID_INFO(GZ_fid).gzfid(4));
[gbl_GZ_FID_INFO(GZ_fid).gzfid data] = GZRead('read',gbl_GZ_FID_INFO(GZ_fid).gzfid,uint32(numtoread),precision_tokens{1}{3},swapbytes);
bytesread = double(gbl_GZ_FID_INFO(GZ_fid).gzfid(4))-startpos;

count_read = fix(bytesread / bytes_per_element);
count_read = skip_interval*fix(count_read/(skip+skip_interval)) + min(skip_interval,mod(count_read,skip+skip_interval));

%% parse data (if skip > 0)
if numtoread~=prod(count)
    skipcount = 0;
    datamask = [];
    while length(datamask)<prod(count)
        datamask = [datamask (1+skipcount*(skip_interval+skip)):(skip_interval+skipcount*(skip_interval+skip))];
        skipcount = skipcount + 1;
    end;
    datamask = datamask(1:prod(count));
    data = data(datamask);
end;

%% reshape data (if length(count) == 2)
if length(count)==2
    data = reshape(data, count);
end;

%% convert data
data=cast(data,precision_tokens{1}{4});

