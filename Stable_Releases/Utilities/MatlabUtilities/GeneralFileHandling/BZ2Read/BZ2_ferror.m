% function message = BZ2_ferror(BZ2_fid)
%
% This function interprets BZ2_error codes -- the input is either a BZ2_fid
% or error code number (values <=0 are interpreted as error codes, > 0 are
% interpreted as BZ2_fid's), and the output is a string.
%
% 06/11/06 CED
% Release Version 1.0

function message = BZ2_ferror(BZ2_fid)

%% declare gloabl variables
global gbl_BZ2_FID_INFO;

%% define defaults
message = 'failed to read error';

%% check input
if ~(isnumeric(BZ2_fid) && length(BZ2_fid)==1 && BZ2_fid==fix(BZ2_fid))
    disp('The input to BZ2_ferror should either be a single BZ2_fid or a bzerror code (negative integer)');
    return;
end;

%% check that BZ2_fid is valid
if BZ2_fid<1
    bzerror = BZ2_fid;
elseif BZ2_fid>length(gbl_BZ2_FID_INFO) || isempty(gbl_BZ2_FID_INFO(BZ2_fid).bzfid)
    disp('In BZ2_ferror:  The BZ2_fid entered does not correspond to an open file');
    return;
else
    bzerror = double(TypeCastAndSwap(gbl_BZ2_FID_INFO(BZ2_fid).bzfid(3),'int32',false,false));
end;

%% interpret error
if bzerror==0
    message = 'BZ_OK';
elseif bzerror==1
    message = 'BZ_RUN_OK';
elseif bzerror==2
    message = 'BZ_FLUSH_OK';
elseif bzerror==3
    message = 'BZ_FINISH_OK';
elseif bzerror==4
    message = 'BZ_STREAM_END';
elseif bzerror==-1
    message = 'BZ_SEQUENCE_ERROR';
elseif bzerror==-2
    message = 'BZ_PARAM_ERROR';
elseif bzerror==-3
    message = 'BZ_MEM_ERROR';
elseif bzerror==-4
    message = 'BZ_DATA_ERROR';
elseif bzerror==-5
    message = 'BZ_DATA_ERROR_MAGIC';
elseif bzerror==-6
    message = 'BZ_IO_ERROR';
elseif bzerror==-7
    message = 'BZ_UNEXPECTED_EOF';
elseif bzerror==-8
    message = 'BZ_OUTBUFF_FULL';
elseif bzerror==-9
    message = 'BZ_CONFIG_ERROR';
else
    message = 'invalid bzerror';
end;

