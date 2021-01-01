% function message = GZ_ferror(GZ_fid)
%
% This function interprets GZ_error codes -- the input is either a GZ_fid
% or error code number (values <=0 are interpreted as error codes, > 0 are
% interpreted as GZ_fid's), and the output is a string.
%
% 06/11/06 CED
% Release Version 1.0

function message = GZ_ferror(GZ_fid)

%% declare gloabl variables
global gbl_GZ_FID_INFO;

%% define defaults
message = 'failed to read error';

%% check input
if ~(isnumeric(GZ_fid) && length(GZ_fid)==1 && GZ_fid==fix(GZ_fid))
    disp('The input to GZ_ferror should either be a single GZ_fid or a gzerror code (negative integer)');
    return;
end;

%% check that GZ_fid is valid
if GZ_fid<1
    gzerror = GZ_fid;
elseif GZ_fid>length(gbl_GZ_FID_INFO) || isempty(gbl_GZ_FID_INFO(GZ_fid).gzfid)
    disp('In GZ_ferror:  The GZ_fid entered does not correspond to an open file');
    return;
else
    gzerror = double(TypeCastAndSwap(gbl_GZ_FID_INFO(GZ_fid).gzfid(3),'int32',false,false));
end;

%% interpret error
if gzerror==0
    message = 'BZ_OK';
elseif gzerror==1
    message = 'BZ_RUN_OK';
elseif gzerror==2
    message = 'BZ_FLUSH_OK';
elseif gzerror==3
    message = 'BZ_FINISH_OK';
elseif gzerror==4
    message = 'BZ_STREAM_END';
elseif gzerror==-1
    message = 'BZ_SEQUENCE_ERROR';
elseif gzerror==-2
    message = 'BZ_PARAM_ERROR';
elseif gzerror==-3
    message = 'BZ_MEM_ERROR';
elseif gzerror==-4
    message = 'BZ_DATA_ERROR';
elseif gzerror==-5
    message = 'BZ_DATA_ERROR_MAGIC';
elseif gzerror==-6
    message = 'BZ_IO_ERROR';
elseif gzerror==-7
    message = 'BZ_UNEXPECTED_EOF';
elseif gzerror==-8
    message = 'BZ_OUTBUFF_FULL';
elseif gzerror==-9
    message = 'BZ_CONFIG_ERROR';
else
    message = 'invalid gzerror';
end;

