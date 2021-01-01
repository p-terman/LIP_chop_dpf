% function status = MEM_fclose(MEM_fid)
% 
% This function is for use with files loaded into memory via MEM_fopen, and
% handles precisely like fclose.  The argument may be any scalar or vector
% of MEM_fid's (given by MEM_fopen, not the same as fid's from fopen), or
% the string 'all', in which case all files opened with MEM_fopen are
% closed.
%
% The function returns 0 if successful, -1 otherwise.  If you pass it a
% vector of MEM_fid's, the output is a corresponding vector of 0's and -1's
%
% 07/22/06 CED
% Maintained in svn repository at
% svn+ssh://lxe@darkanalysis.case.edu/Users/xedmcase/Documents/MatlabCodeRepositories/ReleasedCode

function status = MEM_fclose(MEM_fid)

%% declare gloabl variables
global gbl_MEM_FID_INFO;

%% define defaults
status = -1;

%% check input
if ~((ischar(MEM_fid) && strcmp(MEM_fid,'all')) || (isnumeric(MEM_fid) && ...
        numel(MEM_fid)==length(MEM_fid)))
    disp('The input to MEM_fclose should either be the string "all" or a vector of MEM_fids');
    return;
end;

%% close requested files
if isnumeric(MEM_fid)
    MEM_fid = fix(MEM_fid);
    status = -(ones(1,length(MEM_fid)));
    for f=1:length(MEM_fid)
        if MEM_fid(f)>0 && MEM_fid(f)<=length(gbl_MEM_FID_INFO)
            gbl_MEM_FID_INFO(MEM_fid(f)).data = [];
            gbl_MEM_FID_INFO(MEM_fid(f)).maxblocksize = [];
            status(f) = 0;
        end;
    end;
else
    status = 0;
    for f=1:length(gbl_MEM_FID_INFO)
        gbl_MEM_FID_INFO(f).data = [];
        gbl_MEM_FID_INFO(f).maxblocksize = [];
    end;
end;

%% clean up gbl_MEM_FID_INFO
while length(gbl_MEM_FID_INFO)>0 && isempty(gbl_MEM_FID_INFO(end).maxblocksize)
    gbl_MEM_FID_INFO = gbl_MEM_FID_INFO(1:end-1);
end;
