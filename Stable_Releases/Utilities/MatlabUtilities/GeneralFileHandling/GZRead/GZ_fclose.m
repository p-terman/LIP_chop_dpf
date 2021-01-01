% function status = GZ_fclose(GZ_fid)
% 
% This function interfaces with the GZRead and GZWrite mex files, and
% handles precisely like fclose.  The argument may be any scalar or vector
% of GZ_fid's (given by GZ_fopen, not the same as fid's from fopen), or
% the string 'all', in which case all files opened with GZ_fopen are
% closed.
%
% The function returns 0 if successful, -1 otherwise.  If you pass it a
% vector of GZ_fid's, the output is a corresponding vector of 0's and -1's
%
% write support dependent on existance of GZWrite, which at the time of
% this comment has not been written.
% 06/11/06, CED
% Release Version 1.0

function status = GZ_fclose(GZ_fid)

%% declare gloabl variables
global gbl_GZ_FID_INFO;

%% define defaults
status = -1;

%% check input
if ~((ischar(GZ_fid) && strcmp(GZ_fid,'all')) || (isnumeric(GZ_fid) && ...
        numel(GZ_fid)==length(GZ_fid)))
    disp('The input to GZ_fclose should either be the string "all" or a vector of GZ_fids');
    return;
end;

%% close requested files
if isnumeric(GZ_fid)
    GZ_fid = fix(GZ_fid);
    status = -(ones(1,length(GZ_fid)));
    for f=1:length(GZ_fid)
        if GZ_fid(f)>0 && GZ_fid(f)<=length(gbl_GZ_FID_INFO) && ...
                length(gbl_GZ_FID_INFO(GZ_fid(f)).gzfid)>0
            if strcmp(gbl_GZ_FID_INFO(GZ_fid(f)).mode,'r')
                gbl_GZ_FID_INFO(GZ_fid(f)).gzfid = GZRead('close',gbl_GZ_FID_INFO(GZ_fid(f)).gzfid);
            else
                gbl_GZ_FID_INFO(GZ_fid(f)).gzfid = GZWrite('close',gbl_GZ_FID_INFO(GZ_fid(f)).gzfid);
            end;
            gzerror = double(TypeCastAndSwap(gbl_GZ_FID_INFO(GZ_fid(f)).gzfid(3),'int32',false,false));
            if gzerror < 0
                disp(sprintf('Encountered gzerror %d (%s) in the file with GZ_fid=%d', ...
                    gzerror,GZ_ferror(gzerror),GZ_fid(f)));
            end;
            gbl_GZ_FID_INFO(GZ_fid(f)).gzfid = [];
            status(f) = 0;
        end;
    end;
else
    status = 0;
    for f=1:length(gbl_GZ_FID_INFO)
        if length(gbl_GZ_FID_INFO(f).gzfid)>0
            if strcmp(gbl_GZ_FID_INFO(f).mode,'r')
                gbl_GZ_FID_INFO(f).gzfid = GZRead('close',gbl_GZ_FID_INFO(f).gzfid);
            else
                gbl_GZ_FID_INFO(f).gzfid = GZWrite('close',gbl_GZ_FID_INFO(f).gzfid);
            end;
            gzerror = double(TypeCastAndSwap(gbl_GZ_FID_INFO(f).gzfid(3),'int32',false,false));
            if gzerror < 0
                disp(sprintf('Encountered gzerror %d (%s) in the file with GZ_fid=%d', ...
                    gzerror,GZ_ferror(gzerror),f));
            end;
            gbl_GZ_FID_INFO(f).gzfid = [];
        end;
    end;
end;

%% clean up gbl_GZ_FID_INFO
while length(gbl_GZ_FID_INFO)>0 && isempty(gbl_GZ_FID_INFO(end).gzfid)
    gbl_GZ_FID_INFO = gbl_GZ_FID_INFO(1:end-1);
end;