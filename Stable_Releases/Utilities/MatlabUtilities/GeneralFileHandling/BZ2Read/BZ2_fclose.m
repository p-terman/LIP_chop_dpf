% function status = BZ2_fclose(BZ2_fid)
% 
% This function interfaces with the BZ2Read and BZ2Write mex files, and
% handles precisely like fclose.  The argument may be any scalar or vector
% of BZ2_fid's (given by BZ2_fopen, not the same as fid's from fopen), or
% the string 'all', in which case all files opened with BZ2_fopen are
% closed.
%
% The function returns 0 if successful, -1 otherwise.  If you pass it a
% vector of BZ2_fid's, the output is a corresponding vector of 0's and -1's
%
% write support dependent on existance of BZ2Write, which at the time of
% this comment has not been written.
% 06/11/06, CED
% Release Version 1.0

function status = BZ2_fclose(BZ2_fid)

%% declare gloabl variables
global gbl_BZ2_FID_INFO;

%% define defaults
status = -1;

%% check input
if ~((ischar(BZ2_fid) && strcmp(BZ2_fid,'all')) || (isnumeric(BZ2_fid) && ...
        numel(BZ2_fid)==length(BZ2_fid)))
    disp('The input to BZ2_fclose should either be the string "all" or a vector of BZ2_fids');
    return;
end;

%% close requested files
if isnumeric(BZ2_fid)
    BZ2_fid = fix(BZ2_fid);
    status = -(ones(1,length(BZ2_fid)));
    for f=1:length(BZ2_fid)
        if BZ2_fid(f)>0 && BZ2_fid(f)<=length(gbl_BZ2_FID_INFO) && ...
                length(gbl_BZ2_FID_INFO(BZ2_fid(f)).bzfid)>0
            if strcmp(gbl_BZ2_FID_INFO(BZ2_fid(f)).mode,'r')
                gbl_BZ2_FID_INFO(BZ2_fid(f)).bzfid = BZ2Read('close',gbl_BZ2_FID_INFO(BZ2_fid(f)).bzfid);
            else
                gbl_BZ2_FID_INFO(BZ2_fid(f)).bzfid = BZ2Write('close',gbl_BZ2_FID_INFO(BZ2_fid(f)).bzfid);
            end;
            bzerror = double(TypeCastAndSwap(gbl_BZ2_FID_INFO(BZ2_fid(f)).bzfid(3),'int32',false,false));
            if bzerror < 0
                disp(sprintf('Encountered bzerror %d (%s) in the file with BZ2_fid=%d', ...
                    bzerror,BZ2_ferror(bzerror),BZ2_fid(f)));
            end;
            gbl_BZ2_FID_INFO(BZ2_fid(f)).bzfid = [];
            status(f) = 0;
        end;
    end;
else
    status = 0;
    for f=1:length(gbl_BZ2_FID_INFO)
        if length(gbl_BZ2_FID_INFO(f).bzfid)>0
            if strcmp(gbl_BZ2_FID_INFO(f).mode,'r')
                gbl_BZ2_FID_INFO(f).bzfid = BZ2Read('close',gbl_BZ2_FID_INFO(f).bzfid);
            else
                gbl_BZ2_FID_INFO(f).bzfid = BZ2Write('close',gbl_BZ2_FID_INFO(f).bzfid);
            end;
            bzerror = double(TypeCastAndSwap(gbl_BZ2_FID_INFO(f).bzfid(3),'int32',false,false));
            if bzerror < 0
                disp(sprintf('Encountered bzerror %d (%s) in the file with BZ2_fid=%d', ...
                    bzerror,BZ2_ferror(bzerror),f));
            end;
            gbl_BZ2_FID_INFO(f).bzfid = [];
        end;
    end;
end;

%% clean up gbl_BZ2_FID_INFO
while length(gbl_BZ2_FID_INFO)>0 && isempty(gbl_BZ2_FID_INFO(end).bzfid)
    gbl_BZ2_FID_INFO = gbl_BZ2_FID_INFO(1:end-1);
end;