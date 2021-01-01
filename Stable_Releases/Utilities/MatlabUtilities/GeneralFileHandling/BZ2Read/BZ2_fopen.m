% function [BZ2_fid mode_out machineformat_out] = BZ2_fopen(filename, mode,
%                                                   machineformat, maxbuffersize)
%
% This function interfaces with the BZ2Read and BZ2Write mex files, and 
% handles precisely like fopen, with the following exceptions:
%
% -the fid returned in BZ2_fid has an independent numbering system from
%    that of the normal fopen
%
% -mode may only be read 'r' or append 'a', as these are the only modes
%    supported by the libbz2 -- and thus far only 'r' is implemented.
%
% -machineformat may be true/false (flip endian or not), 1/0 (the same), or
%    'b'/'l' (big/little endian)
%
% -maxbuffersize is the maximum buffer used for bunzip2-ing while
%    navigating the file with BZ2_fseek (can only be used in 'r' mode)
%
% The defaults for the above are read, no endian flip, and 5 MB,
% respectively
%
% Like fopen, this may also be called as:
%   [filename mode machineformat] = BZ2_fopen(BZ2_fid)
%    where BZ2_fid may be a vector (in which case the outputs are cell
%    arrays)
%
% and as
%
% BZ2_fid_list = BZ2_fopen('all')
%
% 06/11/06, CED
% Release Version 1.1
% 1.0 -> 1.1  Fixed machineformat -- in 1.0, it actually opened in the
%               opposite endian-ness from the one specified

function [BZ2_fid mode_out machineformat_out] = BZ2_fopen(filename, mode, machineformat, maxbuffersize)

%% declare gloabl variables
global gbl_BZ2_FID_INFO;
if ~isstruct(gbl_BZ2_FID_INFO)
    gbl_BZ2_FID_INFO = struct('filename',{},'mode',{},'swapbytes',{},'bzfid',{});
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

BZ2_fid = [];
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
    disp('Only r(ead) and a(ppend) modes are supported for BZ2_fopen');
    return;
end;

if strcmp(mode,'a')
    disp('No write support yet for BZ2 files -- if you need this, complain to cdahl@princeton.edu');
    return;
end;

%% if first input is numeric, user is asking for info on BZ2_opened file
if isnumeric(filename) && (numel(filename)==length(filename))
    if length(filename)>1
        for f=1:length(filename)
            filename(f) = fix(filename(f));
            if filename(f)<1 || filename(f)>length(gbl_BZ2_FID_INFO) || ...
                    isempty(gbl_BZ2_FID_INFO(filename(f)).bzfid)
                BZ2_fid{f} = '';
                mode_out{f} = '';
                machineformat_out{f} = '';
            else
                BZ2_fid{f} = gbl_BZ2_FID_INFO(filename(f)).filename;
                mode_out{f} = gbl_BZ2_FID_INFO(filename(f)).mode;
                if xor(currently_big_endian, gbl_BZ2_FID_INFO(filename(f)).swapbytes)
                    machineformat_out{f} = 'b';
                else
                    machineformat_out{f} = 'l';
                end;
            end;
        end;
    else
        filename = fix(filename);
        f = filename;
        if filename<1 || filename>length(gbl_BZ2_FID_INFO) || isempty(gbl_BZ2_FID_INFO(filename).bzfid)
            BZ2_fid = '';
            mode_out = '';
            machineformat_out = '';
        else
            BZ2_fid = gbl_BZ2_FID_INFO(f).filename;
            mode_out = gbl_BZ2_FID_INFO(f).mode;
            if xor(currently_big_endian, gbl_BZ2_FID_INFO(f).swapbytes)
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
    disp('first input to BZ2_fopen must either be a filename, the string "all", or a number or vector of numbers');
    return;
end;

%% if input is the string 'all', return list of valid BZ2_fid's
if strcmp(filename,'all')
    for f=1:length(gbl_BZ2_FID_INFO)
        if ~isempty(gbl_BZ2_FID_INFO(f).bzfid)
            BZ2_fid = [BZ2_fid f];
        end;
    end;
    return;
end;

%% now open file for BZ2 reading or writing (only reading for now)
if length(filename)<4 || ~strcmp(filename(end-3:end),'.bz2')
    if exist(filename,'file') || ~strcmp(mode,'r')
        disp('Warning -- filename passed to BZ2_fopen does not end in .bz2 -- treating as zipped file anyway.');
    else
        filename = [filename '.bz2'];
        disp('File not found, looking for .bz2 version of file');
    end;
end;

if strcmp(mode,'r')
    %% check that file exists
    if ~exist(filename,'file')
        disp(['Requested file ' filename ' for BZ2 reading does not exist.']);
        mode_out = 'file not found';
        return;
    end;
    
    %% assign BZ2_fid to this file
    for f=1:length(gbl_BZ2_FID_INFO)
        if isempty(gbl_BZ2_FID_INFO(f).bzfid)
            BZ2_fid = f;
            break;
        end;
    end;
    if isempty(BZ2_fid)
        BZ2_fid = length(gbl_BZ2_FID_INFO) + 1;
    end;
    
    %% open file and set info in gbl_BZ2_FID_INFO
    gbl_BZ2_FID_INFO(BZ2_fid).filename = filename;
    gbl_BZ2_FID_INFO(BZ2_fid).mode = mode;
    gbl_BZ2_FID_INFO(BZ2_fid).swapbytes = swapbytes;
    gbl_BZ2_FID_INFO(BZ2_fid).bzfid = BZ2Read('open',filename,uint32(maxbuffersize));
    if ischar(gbl_BZ2_FID_INFO(BZ2_fid).bzfid)
        gbl_BZ2_FID_INFO(BZ2_fid).bzfid = [];
        mode_out = 'file not found by BZ2Read';
        disp([filename ' not found by BZ2Read -- try using full path.']);
    end;
    bzerror = double(TypeCastAndSwap(gbl_BZ2_FID_INFO(BZ2_fid).bzfid(3),'int32',false,false));
    if bzerror < 0
        disp(sprintf('Encountered bzerror %d (%s) in the file with BZ2_fid=%d', ...
            bzerror,BZ2_ferror(bzerror),BZ2_fid));
        mode_out = sprintf('bzerror %d', bzerror);
    end;
    return;
end;