% function [MEM_fid mode_out machineformat_out] = MEM_fopen(filename, mode, machineformat, maxfilesize, ignore_mem_error, offset)
%
% This function loads files into memory for later reading, and 
% handles precisely like fopen, with the following exceptions:
%
% -the fid returned in MEM_fid has an independent numbering system from
%    that of the normal fopen
%
% -mode may only be read 'r'
% -machineformat may be true/false (flip endian or not), 1/0 (the same), or
%    'b'/'l' (big/little endian)
%
% -maxfilesize is the maximum buffer used for storing files in memory.  If
% other files are already loaded, maxfilesize will automatically be reduced
% by the memory already taken.
%
% The defaults for the above are read, no endian flip, and 100 MB,
% respectively
%
% Like fopen, this may also be called as:
%   [filename mode machineformat] = MEM_fopen(MEM_fid)
%    where MEM_fid may be a vector (in which case the outputs are cell
%    arrays)
%
% and as
%
% MEM_fid_list = MEM_fopen('all')
%
% 07/22/06, CED
% Maintained in SVN repository at
% svn+ssh://lxe@darkanalysis.case.edu/Users/xedmcase/Documents/MatlabCodeRepositories/ReleasedCode

function [MEM_fid mode_out machineformat_out] = MEM_fopen(filename, mode, machineformat, maxfilesize, ignore_mem_error, offset)

%% declare gloabl variables
global gbl_MEM_FID_INFO;
if ~isstruct(gbl_MEM_FID_INFO)
    gbl_MEM_FID_INFO = struct('filename',{},'mode',{},'swapbytes',{},'filepos',{},'data',{},...
        'offset',{},'ferror',{},'merror',{},'maxblocksize',{},'currentblocksize',{});
end;

%% define defaults
if nargin < 6
    offset = 0;
end;

if nargin < 5
    ignore_mem_error = false;
end;

if nargin < 4
    maxfilesize = 1e8;
end;

if nargin < 3
    machineformat = false;
end;

if nargin < 2
    mode = 'r';
end;

MEM_fid = [];
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

if ~ischar(machineformat)
    if xor(currently_big_endian, swapbytes)
        machineformat = 'b';
    else
        machineformat = 'l';
    end;
end;

if ~(strcmp(mode(1),'r'))
    disp('Only r(ead) mode is supported for MEM_fopen');
    return;
end;

%% adjust maxfilesize for already opened files
if ~ignore_mem_error
    for f=1:length(gbl_MEM_FID_INFO)
        if ~isempty(gbl_MEM_FID_INFO(f).maxblocksize);
            maxfilesize = maxfilesize - gbl_MEM_FID_INFO(f).maxblocksize;
        end;
    end;
end;

%% if first input is numeric, user is asking for info on MEM_opened file
if isnumeric(filename) && (numel(filename)==length(filename))
    if length(filename)>1
        for f=1:length(filename)
            filename(f) = fix(filename(f));
            if filename(f)<1 || filename(f)>length(gbl_MEM_FID_INFO) || ...
                    isempty(gbl_MEM_FID_INFO(filename(f)).maxblocksize)
                MEM_fid{f} = '';
                mode_out{f} = '';
                machineformat_out{f} = '';
            else
                MEM_fid{f} = gbl_MEM_FID_INFO(filename(f)).filename;
                mode_out{f} = gbl_MEM_FID_INFO(filename(f)).mode;
                if xor(currently_big_endian, gbl_MEM_FID_INFO(filename(f)).swapbytes)
                    machineformat_out{f} = 'b';
                else
                    machineformat_out{f} = 'l';
                end;
            end;
        end;
    else
        filename = fix(filename);
        f = filename;
        if filename<1 || filename>length(gbl_MEM_FID_INFO) || isempty(gbl_MEM_FID_INFO(filename).maxblocksize)
            MEM_fid = '';
            mode_out = '';
            machineformat_out = '';
        else
            MEM_fid = gbl_MEM_FID_INFO(f).filename;
            mode_out = gbl_MEM_FID_INFO(f).mode;
            if xor(currently_big_endian, gbl_MEM_FID_INFO(f).swapbytes)
                machineformat_out = 'b';
            else
                machineformat_out = 'l';
            end;
        end;
    end;
    return;
end;

%% if input is the string 'all', return list of valid MEM_fid's
if ischar(filename) && strcmp(filename,'all')
    for f=1:length(gbl_MEM_FID_INFO)
        if ~isempty(gbl_MEM_FID_INFO(f).maxblocksize)
            MEM_fid = [MEM_fid f];
        end;
    end;
    return;
end;

%% if first input is a string, open the file using the appropriate opener
if ischar(filename)
    if ~exist(filename,'file')
        if exist([filename '.bz2'],'file')
            filename = [filename '.bz2'];
        elseif exist([filename '.gz'],'file')
            filename = [filename '.gz'];
        else
            disp('file does not exist');
            return;
        end;
    end;

    if length(filename)>3 && strcmp(filename(end-3:end),'.bz2')
        fid = BZ2_fopen(filename,mode,machineformat);
        filecommand_prefix = 'BZ2_';
    elseif length(filename)>2 && strcmp(filename(end-2:end),'.gz')
        fid = GZ_fopen(filename,mode,machineformat);
        filecommand_prefix = 'GZ_';
    else
        fid = fopen(filename,mode,machineformat);
        filecommand_prefix = '';
    end;
    filename = struct('fid',fid,'filecommand_prefix',filecommand_prefix);
end;

%% if filename is not a dual_fid now, we don't know what to do
if ~isstruct(filename) || ~isfield(filename,'fid') || ~isfield(filename,'filecommand_prefix')
    disp('The first argument to MEM_fopen should either be a MEM_fid array, a single dual_fid (structure with fid and filecommand_prefix fields), a filename, or the string "all"');
    return;
end;

%% now open file for MEM reading or writing (only reading for now)
% assign MEM_fid to this file
for f=1:length(gbl_MEM_FID_INFO)
    if isempty(gbl_MEM_FID_INFO(f).maxblocksize)
        MEM_fid = f;
        break;
    end;
end;
if isempty(MEM_fid) || MEM_fid<1
    MEM_fid = length(gbl_MEM_FID_INFO) + 1;
end;

% read file, save info in gbl_MEM_FID_INFO
[gbl_MEM_FID_INFO(MEM_fid).filename gbl_MEM_FID_INFO(MEM_fid).mode machineformat] = dual_fopen(filename);
if strcmp(machineformat,'b') || strcmp(machineformat,'ieee-be')
    swapbytes = ~currently_big_endian;
elseif strcmp(machineformat,'l') || strcmp(machineformat,'ieee-le')
    swapbytes = currently_big_endian;
end;

gbl_MEM_FID_INFO(MEM_fid).swapbytes = swapbytes;
gbl_MEM_FID_INFO(MEM_fid).filepos = offset;
gbl_MEM_FID_INFO(MEM_fid).offset = offset;
gbl_MEM_FID_INFO(MEM_fid).ferror = '';
gbl_MEM_FID_INFO(MEM_fid).merror = '';
gbl_MEM_FID_INFO(MEM_fid).maxblocksize = maxfilesize;
dual_fseek(filename,offset,'bof');
if maxfilesize > 0
    [gbl_MEM_FID_INFO(MEM_fid).data count] = dual_fread(filename,maxfilesize,'*uint8');
    gbl_MEM_FID_INFO(MEM_fid).ferror = dual_ferror(filename);
    gbl_MEM_FID_INFO(MEM_fid).currentblocksize = count;

%     if length(gbl_MEM_FID_INFO(MEM_fid).data)>count
%         gbl_MEM_FID_INFO(MEM_fid).data(count+1:end) = [];
%     end;

    if isempty(filename.filecommand_prefix)
        fread(filename.fid,1,'*uint8');
        eof_string = ferror(filename.fid);
    elseif strcmp(filename.filecommand_prefix,'GZ_')
        gz_fpos=GZ_ftell(filename.fid);
        GZ_fread(filename.fid,1,'*uint8');
        if gz_fpos==GZ_ftell(filename.fid)
            eof_string='At end-of-file.';
        else
            eof_string='';
        end;
    else
        eof_string = gbl_MEM_FID_INFO(MEM_fid).ferror;
    end;

    if ~(strcmp(eof_string,'BZ_STREAM_END') || strcmp(eof_string,'At end-of-file.'))
        gbl_MEM_FID_INFO(MEM_fid).merror = 'NOT_AT_EOF';
    elseif offset==0
        gbl_MEM_FID_INFO(MEM_fid).maxblocksize = length(gbl_MEM_FID_INFO(MEM_fid).data);
    end;
else
    gbl_MEM_FID_INFO(MEM_fid).data = [];
    gbl_MEM_FID_INFO(MEM_fid).currentblocksize = 0;
    gbl_MEM_FID_INFO(MEM_fid).merror = 'OUT_OF_MEMORY';
    disp('Out of memory for loading files to memory -- clear existing files from memory using MEM_fclose, or use a larger value of maxfilesize');
    gbl_MEM_FID_INFO(MEM_fid).maxblocksize = [];
end;


dual_fclose(filename);
machineformat_out = machineformat;
mode_out = 'r';
return;
