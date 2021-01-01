% function status = BZ2_fseek(BZ2_fid, offset, origin)
% 
% This function interfaces with the BZ2Read mex file, and
% handles precisely like fseek, excpet that 'eof' is not recommended as a
% file origin, as it requires unzipping the entire file to determine the
% file size, and then unzipping again from the beginning to the desired
% location.
%
% The first argument should be a BZ2_fid (given by BZ2_fopen, not the same as fid's from fopen).
%
% The function returns the 0 successful, -1 otherwise otherwise.
%
% 06/12/06 CED
% Release Version 1.0

function status = BZ2_fseek(BZ2_fid, offset, origin)

%% declare gloabl variables
global gbl_BZ2_FID_INFO;

%% define defaults
if nargin<3
    origin = 'bof';
end;

status = -1;

%% check input
if ~(isnumeric(BZ2_fid) && length(BZ2_fid)==1 && BZ2_fid==fix(BZ2_fid) && BZ2_fid>0)
    disp('The first input to BZ2_fseek should be a single BZ2_fid (positive integer)');
    return;
end;

if ~(isnumeric(offset) && length(offset)==1 && offset==fix(offset))
    disp('The second input (offset) to BZ2_fseek should be an integer');
    return;
end;

if ~((isnumeric(origin) && (origin==-1 || origin==0 || origin==1)) || ...
        (ischar(origin) && (strcmp(origin,'bof') || strcmp(origin,'cof') || strcmp(origin,'eof'))))
    disp('The third argument to BZ2_fseek should be either "bof" (or -1), or "cof" (or 0).  "eof" (or 1) is also supported but should be avoided');
    return;
end;

%% check if file is open for reading
if BZ2_fid>length(gbl_BZ2_FID_INFO) || isempty(gbl_BZ2_FID_INFO(BZ2_fid).bzfid) || ~strcmp(gbl_BZ2_FID_INFO(BZ2_fid).mode,'r')
    disp('In BZ2_fseek:  The BZ2_fid entered does not correspond to a file open for reading');
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
    new_position = BZ2_filesize(BZ2_fid) + offset;
else
    new_position = (origin+1)*double(gbl_BZ2_FID_INFO(BZ2_fid).bzfid(4)) + offset;
end;
    
%% set position
gbl_BZ2_FID_INFO(BZ2_fid).bzfid = BZ2Read('setpos',gbl_BZ2_FID_INFO(BZ2_fid).bzfid,uint32(new_position));

%% test new position
if new_position==double(gbl_BZ2_FID_INFO(BZ2_fid).bzfid(4)) && double(TypeCastAndSwap(gbl_BZ2_FID_INFO(BZ2_fid).bzfid(3),'int32',false,false))>=0
    status = 0;
end;

