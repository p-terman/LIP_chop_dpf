% function position = GZ_ftell(GZ_fid)
%
% This function interfaces with the GZRead and GZWrite mex files, and
% handles precisely like ftell.  The argument should be a GZ_fid
% (given by GZ_fopen, not the same as fid's from fopen).
%
% The function returns the poisition in the file if successful, -1
% otherwise.
%
% 06/12/06 CED
% Release Version 1.0

function position = GZ_ftell(GZ_fid)

%% declare gloabl variables
global gbl_GZ_FID_INFO;

%% define defaults
position = -1;

%% check input
if ~(isnumeric(GZ_fid) && length(GZ_fid)==1 && GZ_fid==fix(GZ_fid) && GZ_fid>0)
    disp('The input to GZ_ftell should  be a single GZ_fid (positive integer)');
    return;
end;

%% check if file is open
if GZ_fid>length(gbl_GZ_FID_INFO) || isempty(gbl_GZ_FID_INFO(GZ_fid).gzfid)
    disp('In GZ_ftell:  The GZ_fid entered does not correspond to an open file');
    return;
end;

position = double(gbl_GZ_FID_INFO(GZ_fid).gzfid(4));

