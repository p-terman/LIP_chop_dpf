% function message = MEM_ferror(MEM_fid)
%
% This function reports MEM_error codes -- the input is an MEM_fid
%
% 07/22/06 CED
% Maintained in svn repository at
% svn+ssh://lxe@darkanalysis.case.edu/Users/xedmcase/Documents/MatlabCodeRepositories/ReleasedCode

function message = MEM_ferror(MEM_fid)

%% declare gloabl variables
global gbl_MEM_FID_INFO;

%% define defaults
message = 'failed to read error';

%% check input
if ~(isnumeric(MEM_fid) && length(MEM_fid)==1 && MEM_fid==fix(MEM_fid))
    disp('The input to MEM_ferror should either be a single MEM_fid or a bzerror code (negative integer)');
    return;
end;

%% check that MEM_fid is valid
if MEM_fid<1
    message = 'Invalid MEM_fid';
elseif MEM_fid>length(gbl_MEM_FID_INFO)
    message = 'The MEM_fid entered does not correspond to an open file';
else
    message = ['File error:  ' gbl_MEM_FID_INFO(MEM_fid).ferror ...
        '     Memory error:  ' gbl_MEM_FID_INFO(MEM_fid).merror];
    if isempty(gbl_MEM_FID_INFO(MEM_fid).maxblocksize)
        message = [message '(File no longer open)'];
    end;
end;

