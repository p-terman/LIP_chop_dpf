% function position = BZ2_ftell(BZ2_fid)
%
% This function interfaces with the BZ2Read and BZ2Write mex files, and
% handles precisely like ftell.  The argument should be a BZ2_fid
% (given by BZ2_fopen, not the same as fid's from fopen).
%
% The function returns the poisition in the file if successful, -1
% otherwise.
%
% 06/12/06 CED
% Release Version 1.0

function position = BZ2_ftell(BZ2_fid)

%% declare gloabl variables
global gbl_BZ2_FID_INFO;

%% define defaults
position = -1;

%% check input
if ~(isnumeric(BZ2_fid) && length(BZ2_fid)==1 && BZ2_fid==fix(BZ2_fid) && BZ2_fid>0)
    disp('The input to BZ2_ftell should  be a single BZ2_fid (positive integer)');
    return;
end;

%% check if file is open
if BZ2_fid>length(gbl_BZ2_FID_INFO) || isempty(gbl_BZ2_FID_INFO(BZ2_fid).bzfid)
    disp('In BZ2_ftell:  The BZ2_fid entered does not correspond to an open file');
    return;
end;

position = double(gbl_BZ2_FID_INFO(BZ2_fid).bzfid(4));

