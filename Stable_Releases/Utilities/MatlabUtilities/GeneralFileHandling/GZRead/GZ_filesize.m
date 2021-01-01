% function filesize = GZ_filesize(GZ_fid)
%
% This function interfaces with GZRead to return the size of a GZ_file
% opened for reading.  It does this by unzipping the file to the end of the
% file, so should be avoided when processing many large files.  Afterwards,
% the file position is reset to the start of the file.  GZ_fid may also be
% a filename, in which case the file is GZ_fopen'ed, read through, and
% then GZ_fclose'd.
%
% 06/15/06, CED
% Release Version 1.0

function filesize = GZ_filesize(GZ_fid)

%% declare gloabl variables
global gbl_GZ_FID_INFO;

%% set defaults
filesize = -1;

%% check input
if ~((isnumeric(GZ_fid) && length(GZ_fid)==1 && GZ_fid==fix(GZ_fid) && GZ_fid>0) || ischar(GZ_fid))
    disp('The first input to GZ_filsize should be a single GZ_fid (positive integer), or a filename');
    return;
end;

%% check if file is open for reading
if ischar(GZ_fid)
    closefile = true;
    GZ_fid = GZ_fopen(GZ_fid);
else
    closefile = false;
end;

if GZ_fid>length(gbl_GZ_FID_INFO) || isempty(gbl_GZ_FID_INFO(GZ_fid).gzfid) || ~strcmp(gbl_GZ_FID_INFO(GZ_fid).mode,'r')
    disp('In GZ_filesize:  The GZ_fid entered does not correspond to a open file for reading');
    return;
end;

%% go to end of file
gbl_GZ_FID_INFO(GZ_fid).gzfid = GZRead('setpos',gbl_GZ_FID_INFO(GZ_fid).gzfid,uint32(2^31-1),true);

% %% verify at end of file
% if ~strcmp(GZ_ferror(GZ_fid),'BZ_STREAM_END');
%     disp('GZ error in GZ_filesize -- did not reach end of file');
%     return;
% end;

%% read filesize
filesize = GZ_ftell(GZ_fid);

%% reset file position to start of file, or close the file
if closefile
    GZ_fclose(GZ_fid);
else
    GZ_fseek(GZ_fid, 0, 'bof');
end;

