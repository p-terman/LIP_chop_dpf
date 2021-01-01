% function status = GZ_fseek(GZ_fid, offset, origin)
% 
% This function interfaces with the GZRead mex file, and
% handles precisely like fseek, excpet that 'eof' is not recommended as a
% file origin, as it requires unzipping the entire file to determine the
% file size, and then unzipping again from the beginning to the desired
% location.
%
% The first argument should be a GZ_fid (given by GZ_fopen, not the same as fid's from fopen).
%
% The function returns the 0 successful, -1 otherwise otherwise.
%
% 06/12/06 CED
% Release Version 1.0

function status = GZ_fseek(GZ_fid, offset, origin)

%% declare gloabl variables
global gbl_GZ_FID_INFO;

%% define defaults
if nargin<3
    origin = 'bof';
end;

status = -1;

%% check input
if ~(isnumeric(GZ_fid) && length(GZ_fid)==1 && GZ_fid==fix(GZ_fid) && GZ_fid>0)
    disp('The first input to GZ_fseek should be a single GZ_fid (positive integer)');
    return;
end;

if ~(isnumeric(offset) && length(offset)==1 && offset==fix(offset))
    disp('The second input (offset) to GZ_fseek should be an integer');
    return;
end;

if ~((isnumeric(origin) && (origin==-1 || origin==0 || origin==1)) || ...
        (ischar(origin) && (strcmp(origin,'bof') || strcmp(origin,'cof') || strcmp(origin,'eof'))))
    disp('The third argument to GZ_fseek should be either "bof" (or -1), or "cof" (or 0).  "eof" (or 1) is also supported but should be avoided');
    return;
end;

%% check if file is open for reading
if GZ_fid>length(gbl_GZ_FID_INFO) || isempty(gbl_GZ_FID_INFO(GZ_fid).gzfid) || ~strcmp(gbl_GZ_FID_INFO(GZ_fid).mode,'r')
    disp('In GZ_fseek:  The GZ_fid entered does not correspond to a file open for reading');
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
    new_position = GZ_filesize(GZ_fid) + offset;
else
    new_position = (origin+1)*double(gbl_GZ_FID_INFO(GZ_fid).gzfid(4)) + offset;
end;
    
%% set position
gbl_GZ_FID_INFO(GZ_fid).gzfid = GZRead('setpos',gbl_GZ_FID_INFO(GZ_fid).gzfid,uint32(new_position));

%% test new position
if new_position==double(gbl_GZ_FID_INFO(GZ_fid).gzfid(4))% && double(TypeCastAndSwap(gbl_GZ_FID_INFO(GZ_fid).gzfid(3),'int32',false,false))>=0
    status = 0;
end;

