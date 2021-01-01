% function filesize = BZ2_filesize(BZ2_fid)
%
% This function interfaces with BZ2Read to return the size of a BZ2_file
% opened for reading.  It does this by unzipping the file to the end of the
% file, so should be avoided when processing many large files.  Afterwards,
% the file position is reset to the start of the file.  BZ2_fid may also be
% a filename, in which case the file is BZ2_fopen'ed, read through, and
% then BZ2_fclose'd.
%
% 06/15/06, CED
% Release Version 1.0

function filesize = BZ2_filesize(BZ2_fid)

%% declare gloabl variables
global gbl_BZ2_FID_INFO;

%% set defaults
filesize = -1;

%% check input
if ~((isnumeric(BZ2_fid) && length(BZ2_fid)==1 && BZ2_fid==fix(BZ2_fid) && BZ2_fid>0) || ischar(BZ2_fid))
    disp('The first input to BZ2_filsize should be a single BZ2_fid (positive integer), or a filename');
    return;
end;

%% check if file is open for reading
if ischar(BZ2_fid)
    closefile = true;
    BZ2_fid = BZ2_fopen(BZ2_fid);
else
    closefile = false;
end;

if BZ2_fid>length(gbl_BZ2_FID_INFO) || isempty(gbl_BZ2_FID_INFO(BZ2_fid).bzfid) || ~strcmp(gbl_BZ2_FID_INFO(BZ2_fid).mode,'r')
    disp('In BZ2_filesize:  The BZ2_fid entered does not correspond to a open file for reading');
    return;
end;

%% go to end of file
gbl_BZ2_FID_INFO(BZ2_fid).bzfid = BZ2Read('setpos',gbl_BZ2_FID_INFO(BZ2_fid).bzfid,uint32(1e99),true);

%% verify at end of file
if ~strcmp(BZ2_ferror(BZ2_fid),'BZ_STREAM_END');
    disp('BZ2 error in BZ2_filesize -- did not reach end of file');
    return;
end;

%% read filesize
filesize = BZ2_ftell(BZ2_fid);

%% reset file position to start of file, or close the file
if closefile
    BZ2_fclose(BZ2_fid);
else
    BZ2_fseek(BZ2_fid, 0, 'bof');
end;

