% function status = MEM_fseek(MEM_fid, offset, origin)
% 
% This function is used with files loaded to memory via MEM_fopen, and
% handles precisely like fseek, excpet that 'eof' is not recommended as a
% file origin for zipped files, as it may require unzipping the entire file
% to determine the file size, and then unzipping again from the beginning 
% to the desired location.  (If the file is completely loaded in memory,
% this is not a problem -- but if maxblocksize is smaller than the file
% size, 'eof' is not recommended for zipped files).
%
% The first argument should be a MEM_fid (given by MEM_fopen, not the same as fid's from fopen).
%
% The function returns the 0 successful, -1 otherwise otherwise.
%
% 07/22/06 CED
% Maintained in svn repository at
% svn+ssh://lxe@darkanalysis.case.edu/Users/xedmcase/Documents/MatlabCodeRepositories/ReleasedCode

function status = MEM_fseek(MEM_fid, offset, origin)

%% declare gloabl variables
global gbl_MEM_FID_INFO;

%% define defaults
if nargin<3
    origin = 'bof';
end;

status = -1;

%% check input
if ~(isnumeric(MEM_fid) && length(MEM_fid)==1 && MEM_fid==fix(MEM_fid) && MEM_fid>0)
    disp('The first input to MEM_fseek should be a single MEM_fid (positive integer)');
    return;
end;

if ~(isnumeric(offset) && length(offset)==1 && offset==fix(offset))
    disp('The second input (offset) to MEM_fseek should be an integer');
    return;
end;

if ~((isnumeric(origin) && (origin==-1 || origin==0 || origin==1)) || ...
        (ischar(origin) && (strcmp(origin,'bof') || strcmp(origin,'cof') || strcmp(origin,'eof'))))
    disp('The third argument to MEM_fseek should be either "bof" (or -1), or "cof" (or 0).  "eof" (or 1) is also supported but should be avoided');
    return;
end;

%% check if file is open for reading
if MEM_fid>length(gbl_MEM_FID_INFO) || isempty(gbl_MEM_FID_INFO(MEM_fid).maxblocksize)
    disp('In MEM_fseek:  The MEM_fid entered does not correspond to an open file');
    return;
end;

%% determine new position
if ischar(origin)
    if strcmp(origin,'bof')
        origin = -1;
    elseif strcmp(origin,'eof')
        origin = 1;
    else
        origin = 0;
    end;
end;

if origin==1
    if strcmp(gbl_MEM_FID_INFO(MEM_fid).merror,'NOT_AT_EOF')
        if length(gbl_MEM_FID_INFO(MEM_fid).filename)>3 && strcmp(gbl_MEM_FID_INFO(MEM_fid).filename(end-3:end),'.bz2')
            filesize = BZ2_filesize(gbl_MEM_FID_INFO(MEM_fid).filename);
        elseif length(gbl_MEM_FID_INFO(MEM_fid).filename)>2 && strcmp(gbl_MEM_FID_INFO(MEM_fid).filename(end-3:end),'.gz')
            filesize = GZ_filesize(gbl_MEM_FID_INFO(MEM_fid).filename);
        else
            fileinfo = dir(gbl_MEM_FID_INFO(MEM_fid).filename);
            filesize = fileinfo.bytes;
        end;
        new_position = filesize + offset;
    else
        new_position = gbl_MEM_FID_INFO(MEM_fid).offset + gbl_MEM_FID_INFO(MEM_fid).currentblocksize + offset;
    end;
else
    new_position = (origin+1)*gbl_MEM_FID_INFO(MEM_fid).filepos + offset;
end;
    
%% set position
if new_position >= gbl_MEM_FID_INFO(MEM_fid).offset && ...
        (new_position < (gbl_MEM_FID_INFO(MEM_fid).offset + gbl_MEM_FID_INFO(MEM_fid).currentblocksize) || ...
        isempty(gbl_MEM_FID_INFO(MEM_fid).data) || ~strcmp(gbl_MEM_FID_INFO(MEM_fid).merror,'NOT_AT_EOF'))
    gbl_MEM_FID_INFO(MEM_fid).filepos = new_position;
    status = 0;
else
    status = MEM_freopen(MEM_fid,new_position);
    gbl_MEM_FID_INFO(MEM_fid).filepos = gbl_MEM_FID_INFO(MEM_fid).offset;
end;

