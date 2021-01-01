% function [status] = MEM_freopen(MEM_fid, offset)
%
% This function will reopen a file stored in memory, clearing the part
% currently in memory and replacing it with a block starting at offset and
% going to the eof or at most maxblocksize bytes (maxblocksize is set in
% the initial MEM_fopen command)
%
% 07/22/06, CED
% Maintained in SVN repository at
% svn+ssh://lxe@darkanalysis.case.edu/Users/xedmcase/Documents/MatlabCodeRepositories/ReleasedCode

function [status] = MEM_freopen(MEM_fid, offset)

%% declare gloabl variables
global gbl_MEM_FID_INFO;
if ~isstruct(gbl_MEM_FID_INFO)
    gbl_MEM_FID_INFO = struct('filename',{},'mode',{},'swapbytes',{},'filepos',{},'data',{},...
        'offset',{},'ferror',{},'merror',{},'maxblocksize',{},'currentblocksize',{});
end;

%% define defaults
if nargin < 2
    offset = 0;
end;

status = -1;

%% check for open MEM_fid
if MEM_fid<0 || MEM_fid>length(gbl_MEM_FID_INFO) || isempty(gbl_MEM_FID_INFO(MEM_fid).maxblocksize)
    disp('In MEM_freopen:  MEM_fid does not correspond to open file');
    return;
end;

%% check this machine's endian-ness
endiantest = hex2dec('01020304');
bytewise_endiantest = double(TypeCastAndSwap(uint32(endiantest),'uint8',false,false));
currently_big_endian = (bytewise_endiantest(1)==1);

if xor(currently_big_endian, gbl_MEM_FID_INFO(MEM_fid).swapbytes)
    machineformat = 'b';
else
    machineformat = 'l';
end;
mode = 'r';

%% open the file using the appropriate opener
filename = gbl_MEM_FID_INFO(MEM_fid).filename;
if ~exist(filename,'file')
    if exist([filename '.bz2'],'file')
        filename = [filename '.bz2'];
    elseif exist([filename '.gz'],'file')
        filename = [filename '.gz'];
    else
        disp('In MEM_freopen:  file no longer exists');
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
gbl_MEM_FID_INFO(MEM_fid).filename = filename;
dfid = struct('fid',fid,'filecommand_prefix',filecommand_prefix);

%% now open file for BZ2 reading or writing (only reading for now)
% read file, save info in gbl_MEM_FID_INFO
% gbl_MEM_FID_INFO(MEM_fid).filepos = offset;
gbl_MEM_FID_INFO(MEM_fid).offset = offset;
gbl_MEM_FID_INFO(MEM_fid).ferror = '';
gbl_MEM_FID_INFO(MEM_fid).merror = '';
dual_fseek(dfid,offset,'bof');

[gbl_MEM_FID_INFO(MEM_fid).data count] = dual_fread(dfid,gbl_MEM_FID_INFO(MEM_fid).maxblocksize,'*uint8');
gbl_MEM_FID_INFO(MEM_fid).ferror = dual_ferror(dfid);
gbl_MEM_FID_INFO(MEM_fid).currentblocksize = count;

% if length(gbl_MEM_FID_INFO(MEM_fid).data)>count
%     gbl_MEM_FID_INFO(MEM_fid).data(count+1:end) = [];
% end;

if isempty(dfid.filecommand_prefix)
    fread(dfid.fid,1,'*uint8');
    eof_string = ferror(dfid.fid);
elseif strcmp(dfid.filecommand_prefix,'GZ_')
    gz_fpos=GZ_ftell(dfid.fid);
    GZ_fread(dfid.fid,1,'*uint8');
    if gz_fpos==GZ_ftell(dfid.fid)
        eof_string='At end-of-file.';
    else
        eof_string='';
    end;
else
    eof_string = gbl_MEM_FID_INFO(MEM_fid).ferror;
end;

if ~(strcmp(eof_string,'BZ_STREAM_END') || strcmp(eof_string,'At end-of-file.'))
    gbl_MEM_FID_INFO(MEM_fid).merror = 'NOT_AT_EOF';
end;

dual_fclose(dfid);
status = 0;
