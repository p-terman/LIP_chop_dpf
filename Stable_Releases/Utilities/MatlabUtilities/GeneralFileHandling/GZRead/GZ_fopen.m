% function [GZ_fid mode_out machineformat_out] = GZ_fopen(filename, mode,
%                                                   machineformat, maxbuffersize)
%
% This function interfaces with the GZRead and GZWrite mex files, and 
% handles precisely like fopen, with the following exceptions:
%
% -the fid returned in GZ_fid has an independent numbering system from
%    that of the normal fopen
%
% -mode may only be read 'r' or append 'a', as these are the only modes
%    supported by the zlib -- and thus far only 'r' is implemented.
%
% -machineformat may be true/false (flip endian or not), 1/0 (the same), or
%    'b'/'l' (big/little endian)
%
% -maxbuffersize is the maximum buffer used for bunzip2-ing while
%    navigating the file with GZ_fseek (can only be used in 'r' mode)
%
% The defaults for the above are read, no endian flip, and 5 MB,
% respectively
%
% Like fopen, this may also be called as:
%   [filename mode machineformat] = GZ_fopen(GZ_fid)
%    where GZ_fid may be a vector (in which case the outputs are cell
%    arrays)
%
% and as
%
% GZ_fid_list = GZ_fopen('all')
%
% 06/11/06, CED
% Release Version 1.1
% 1.0 -> 1.1  Fixed machineformat -- in 1.0, it actually opened in the
%               opposite endian-ness from the one specified

function [GZ_fid mode_out machineformat_out] = GZ_fopen(filename, mode, machineformat, maxbuffersize)

%% declare gloabl variables
global gbl_GZ_FID_INFO;
if ~isstruct(gbl_GZ_FID_INFO)
    gbl_GZ_FID_INFO = struct('filename',{},'mode',{},'swapbytes',{},'gzfid',{});
end;

%% define defaults
if nargin < 4
    maxbuffersize = 5e6;
end;

if nargin < 3
    machineformat = false;
end;

if nargin < 2
    mode = 'r';
end;

GZ_fid = [];
mode_out = '';
machineformat_out = [];

%% check inputs for supported modes
endiantest = hex2dec('01020304');
bytewise_endiantest = double(TypeCastAndSwap(uint32(endiantest),'uint8',false,false));
currently_big_endian = (bytewise_endiantest(1)==1);

if islogical(machineformat)
    swapbytes = machineformat(1);
elseif isnumeric(machineformat)
    swapbytes = logical(machineformat(1));
elseif ischar(machineformat)
    if strcmp(machineformat,'b')
        swapbytes = ~currently_big_endian;
    elseif strcmp(machineformat,'l')
        swapbytes = currently_big_endian;
    else
        disp('invalid machine format -- char machine formats must be l or b');
        return;
    end;
else
    disp('invalid machine format, format must be logical (swapbytes or not), numeric (1 or 0), or char (l or b)');
    return;
end;

if ~(strcmp(mode,'a') || strcmp(mode,'r'))
    disp('Only r(ead) and a(ppend) modes are supported for GZ_fopen');
    return;
end;

if strcmp(mode,'a')
    disp('No write support yet for GZ files -- if you need this, complain to cdahl@princeton.edu');
    return;
end;

%% if first input is numeric, user is asking for info on GZ_opened file
if isnumeric(filename) && (numel(filename)==length(filename))
    if length(filename)>1
        for f=1:length(filename)
            filename(f) = fix(filename(f));
            if filename(f)<1 || filename(f)>length(gbl_GZ_FID_INFO) || ...
                    isempty(gbl_GZ_FID_INFO(filename(f)).gzfid)
                GZ_fid{f} = '';
                mode_out{f} = '';
                machineformat_out{f} = '';
            else
                GZ_fid{f} = gbl_GZ_FID_INFO(filename(f)).filename;
                mode_out{f} = gbl_GZ_FID_INFO(filename(f)).mode;
                if xor(currently_big_endian, gbl_GZ_FID_INFO(filename(f)).swapbytes)
                    machineformat_out{f} = 'b';
                else
                    machineformat_out{f} = 'l';
                end;
            end;
        end;
    else
        filename = fix(filename);
        f = filename;
        if filename<1 || filename>length(gbl_GZ_FID_INFO) || isempty(gbl_GZ_FID_INFO(filename).gzfid)
            GZ_fid = '';
            mode_out = '';
            machineformat_out = '';
        else
            GZ_fid = gbl_GZ_FID_INFO(f).filename;
            mode_out = gbl_GZ_FID_INFO(f).mode;
            if xor(currently_big_endian, gbl_GZ_FID_INFO(f).swapbytes)
                machineformat_out = 'b';
            else
                machineformat_out = 'l';
            end;
        end;
    end;
    return;
end;

%% if first input is not a string, we don't know what to do with it now
if ~ischar(filename)
    disp('first input to GZ_fopen must either be a filename, the string "all", or a number or vector of numbers');
    return;
end;

%% if input is the string 'all', return list of valid GZ_fid's
if strcmp(filename,'all')
    for f=1:length(gbl_GZ_FID_INFO)
        if ~isempty(gbl_GZ_FID_INFO(f).gzfid)
            GZ_fid = [GZ_fid f];
        end;
    end;
    return;
end;

%% now open file for GZ reading or writing (only reading for now)
if length(filename)<3 || ~strcmp(filename(end-2:end),'.gz')
    if exist(filename,'file') || ~strcmp(mode,'r')
        disp('Warning -- filename passed to GZ_fopen does not end in .gz -- treating as zipped file anyway.');
    else
        filename = [filename '.gz'];
        disp('File not found, looking for .gz version of file');
    end;
end;

if strcmp(mode,'r')
    %% check that file exists
    if ~exist(filename,'file')
        disp(['Requested file ' filename ' for GZ reading does not exist.']);
        mode_out = 'file not found';
        return;
    end;
    
    %% assign GZ_fid to this file
    for f=1:length(gbl_GZ_FID_INFO)
        if isempty(gbl_GZ_FID_INFO(f).gzfid)
            GZ_fid = f;
            break;
        end;
    end;
    if isempty(GZ_fid)
        GZ_fid = length(gbl_GZ_FID_INFO) + 1;
    end;
    
    %% open file and set info in gbl_GZ_FID_INFO
    gbl_GZ_FID_INFO(GZ_fid).filename = filename;
    gbl_GZ_FID_INFO(GZ_fid).mode = mode;
    gbl_GZ_FID_INFO(GZ_fid).swapbytes = swapbytes;
    gbl_GZ_FID_INFO(GZ_fid).gzfid = GZRead('open',filename,uint32(maxbuffersize));
    if ischar(gbl_GZ_FID_INFO(GZ_fid).gzfid)
        gbl_GZ_FID_INFO(GZ_fid).gzfid = [];
        mode_out = 'file not found by GZRead';
        disp([filename ' not found by GZRead -- try using full path.']);
    end;
    gzerror = double(TypeCastAndSwap(gbl_GZ_FID_INFO(GZ_fid).gzfid(3),'int32',false,false));
    if gzerror < 0
        disp(sprintf('Encountered gzerror %d (%s) in the file with GZ_fid=%d', ...
            gzerror,GZ_ferror(gzerror),GZ_fid));
        mode_out = sprintf('gzerror %d', gzerror);
    end;
    return;
end;