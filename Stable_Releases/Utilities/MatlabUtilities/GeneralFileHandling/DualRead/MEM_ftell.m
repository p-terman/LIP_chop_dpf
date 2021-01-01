% function position = MEM_ftell(MEM_fid)
%
% This function is used with files stored in memory via MEM_fopen, and
% handles precisely like ftell.  The argument should be a MEM_fid
% (given by MEM_fopen, not the same as fid's from fopen).
%
% The function returns the poisition in the file if successful, -1
% otherwise.
%
% 07/22/06 CED
% Maintained in svn repository at
% svn+ssh://lxe@darkanalysis.case.edu/Users/xedmcase/Documents/MatlabCodeRepositories/ReleasedCode

function position = MEM_ftell(MEM_fid)

%% declare gloabl variables
global gbl_MEM_FID_INFO;

%% define defaults
position = -1;

%% check input
if ~(isnumeric(MEM_fid) && length(MEM_fid)==1 && MEM_fid==fix(MEM_fid) && MEM_fid>0)
    disp('The input to MEM_ftell should  be a single MEM_fid (positive integer)');
    return;
end;

%% check if file is open
if MEM_fid>length(gbl_MEM_FID_INFO) || isempty(gbl_MEM_FID_INFO(MEM_fid).maxblocksize)
    disp('In MEM_ftell:  The MEM_fid entered does not correspond to an open file');
    return;
end;

position = gbl_MEM_FID_INFO(MEM_fid).filepos;

