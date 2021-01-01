% function [fid_info_out mode_out machineformat_out] = dual_fopen(fid_info)
%
% If passed an fid, calls the appropriate fopen function to get file info,
% if passed the string 'all', calls all fopen functions and returns list of
% fid's with filecommand_prefix'es
%
% This code cannot be used to actually open a file -- use OpenBinaryFile instead.
%
% CED
% Release Version 1.0

function [fid_info_out mode_out machineformat_out] = dual_fopen(fid_info)

%% set defaults
mode_out = '';
machineformat_out = '';
fid_info_out = [];

%% check input
if isstruct(fid_info)
%% get fid info
    fid = fid_info.fid;
    prefix = fid_info.filecommand_prefix;

%% make call
    switch prefix
        case ''
            [fid_info_out mode_out machineformat_out] = fopen(fid);
        case 'BZ2_'
            [fid_info_out mode_out machineformat_out] = BZ2_fopen(fid);
        case 'GZ_'
            [fid_info_out mode_out machineformat_out] = GZ_fopen(fid);
        case 'MEM_'
            [fid_info_out mode_out machineformat_out] = MEM_fopen(fid);
        otherwise
            disp('unknown filecommand_prefix in dual_fopen');
    end;
    
%% otherwise get all fid lists
elseif ischar(fid_info) && strcmp(fid_info,'all')
    mode_out = '';
    machineformat_out = '';
    fids = fopen('all');
    BZ2_fids = BZ2_fopen('all');
    GZ_fids = GZ_fopen('all');
    MEM_fids = MEM_fopen('all');
    for f=1:length(fids)
        fid_info_out(f).fid = fids(f);
        fid_info_out(f).filecommand_prefix = '';
    end;
    for f=1:length(BZ2_fids)
        fid_info_out(f+length(fids)).fid = BZ2_fids(f);
        fid_info_out(f+length(fids)).filecommand_prefix = 'BZ2_';
    end;
    for f=1:length(GZ_fids)
        fid_info_out(f+length(fids)).fid = GZ_fids(f);
        fid_info_out(f+length(fids)).filecommand_prefix = 'GZ_';
    end;
    for f=1:length(MEM_fids)
        fid_info_out(f+length(fids)).fid = MEM_fids(f);
        fid_info_out(f+length(fids)).filecommand_prefix = 'MEM_';
    end;
else
    disp('Invalid input to dual_fopen -- to actually open a file, use OpenBinaryFile');
end;