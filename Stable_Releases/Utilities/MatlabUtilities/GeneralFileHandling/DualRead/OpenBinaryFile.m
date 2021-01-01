% function [dual_fid endianflag] = OpenBinaryFile(filename, load_to_memory, memory_cap)
%
% OpenBinaryFile opens a binary file for reading and checks the
% endian-ness.  If the first four bytes of the are not 0x01020304 in either
% endian-ness, the file is opened as little endian, with file pointer at
% beginning of file.
% Otherwise, the file is opened with the correct endian-ness and the
% dual_fid is returned -- the current position of the file is on the 5th
% byte (after the endian-ness indicator).
%
%  Input:  
%          filename       - the full filename, as you would put into fopen
%          load_to_memory - if set to true, file will be fully loaded to
%                           memory, then closed.  dual_fxxx routines will
%                           work as usual, with filecommand_prefix = 'mem_'
%                           default is false.
%          memory_cap     - the maximum memory allowed, in total, for files
%                           stored in memory -- if this memory is already
%                           filled, file will not be loaded to memory even
%                           if load_to_memory is set to true.  (To check
%                           memory-loaded status, look to see whether
%                           fid.filecommand_prefix is 'MEM_')
%
%  Output:
%          dual_fid      - the dual_format fid for the open file (see below)
%          endianflag    - 'l' or 'b', depending how file is opened
%
%  06/15/06 update -- now also opens zipped binary files.  If you specify a
%  filename ending in .bz2, it will automatically use BZ2_fopen to open the
%  file.  If you don't specify .bz2, and the desired file does not exist,
%  .bz2 is appended to the filename and BZ2_fopen is used.
%
%  With this change, the fid output is now a struct containing two fields:
%  fid.fid is the fid number for the opened file
%  fid.filecommand_prefix is the string prefix to the file commands used
%  with this file type (i.e., '' or 'BZ2_')
%
%  CED
%  Release Version 1.0
%
%  05/12/08, Version 1.1 (committed to LUXcode)
%            If endianness indicator is not found, file is opened as
%            little-endian with file position at beginning of file

function [dual_fid endianflag] = OpenBinaryFile(filename, load_to_memory, memory_cap)

%% set defaults
if nargin<3
    memory_cap = 100*2^20; % 100 MB
end

if nargin<2
    load_to_memory = false;
end

dual_fid.fid = -1;
dual_fid.filecommand_prefix = '';
endianflag = '';

%% check for file format
if length(filename)>4 && strcmp(filename(end-3:end),'.bz2')
    filecommand = 'BZ2_';
elseif length(filename)>3 && strcmp(filename(end-2:end),'.gz')
    filecommand = 'GZ_';
elseif ~exist(filename,'file') && exist([filename '.bz2'],'file')
    filecommand = 'BZ2_';
    filename = [filename '.bz2'];
elseif ~exist(filename,'file') && exist([filename '.gz'],'file')
    filecommand = 'GZ_';
    filename = [filename '.gz'];
else
    filecommand = '';
end

if ~exist(filename,'file')
    disp('file not found by OpenBinaryFile');
    return
end

%% create handles for file handling (no more eval statements)
fopen_command = str2func([filecommand 'fopen']);
fclose_command = str2func([filecommand 'fclose']);
fread_command = str2func([filecommand 'fread']);

%% set some constants
endiantest = hex2dec('01020304');

endianflag = 'l';
modeflag = 'r';
endiantesttype = 'uint32';

%% open file and check endianness
fid = fopen_command(filename,modeflag,endianflag);
endianresult = fread_command(fid,1,endiantesttype);
if endiantest~=endianresult
    fclose_command(fid);
    endianflag = 'b';
    fid = fopen_command(filename,modeflag,endianflag);
    endianresult = fread_command(fid,1,endiantesttype);
    if endiantest~=endianresult
        fclose_command(fid);
        endianflag = 'l'; % default is little-endian
        fid = fopen_command(filename,modeflag,endianflag);
    end
end
dual_fid.fid = fid;
dual_fid.filecommand_prefix = filecommand;

%% load to memory if desired
if load_to_memory && memory_cap>MEM_MemoryUsed()
    dual_fid.fid = MEM_fopen(dual_fid,'r',0,memory_cap);
    dual_fid.filecommand_prefix = 'MEM_';
end
