% function [data count_read] = dual_fread(fid_info, count, precision, skip,
%                                           machineformat)
%
% Calls the appropriate fread, given by fid_info.filecommand_prefix
%
% CED
% Release Version 1.0

function [data count_read] = dual_fread(fid_info, count, precision, skip, machineformat)

%% declare global variables
global gbl_MEM_FID_INFO

%% make call
switch fid_info.filecommand_prefix
    case 'MEM_'
        switch nargin
            case 3
                [data count_read bytes_read] = MEM_fread_fast(fid_info.fid, count, precision);
                if isempty(bytes_read)
                    [data count_read] = MEM_fread(fid_info.fid, count, precision);
                else
                    gbl_MEM_FID_INFO(fid_info.fid).filepos = gbl_MEM_FID_INFO(fid_info.fid).filepos + bytes_read;
                end;
            case 2
                [data count_read bytes_read] = MEM_fread_fast(fid_info.fid, count);
                if isempty(bytes_read)
                    [data count_read] = MEM_fread(fid_info.fid, count);
                else
                    gbl_MEM_FID_INFO(fid_info.fid).filepos = gbl_MEM_FID_INFO(fid_info.fid).filepos + bytes_read;
                end;
            case 1
                [data count_read] = MEM_fread(fid_info.fid);
            case 4
                [data count_read] = MEM_fread(fid_info.fid, count, precision, skip);
            otherwise
                [data count_read] = MEM_fread(fid_info.fid, count, precision, skip, machineformat);
        end;
    case ''
        switch nargin
            case 3
                [data count_read] = fread(fid_info.fid, count, precision);
            case 2
                [data count_read] = fread(fid_info.fid, count);
            case 1
                [data count_read] = fread(fid_info.fid);
            case 4
                [data count_read] = fread(fid_info.fid, count, precision, skip);
            otherwise
                [data count_read] = fread(fid_info.fid, count, precision, skip, machineformat);
        end;
    case 'BZ2_'
        switch nargin
            case 3
                [data count_read] = BZ2_fread(fid_info.fid, count, precision);
            case 2
                [data count_read] = BZ2_fread(fid_info.fid, count);
            case 1
                [data count_read] = BZ2_fread(fid_info.fid);
            case 4
                [data count_read] = BZ2_fread(fid_info.fid, count, precision, skip);
            otherwise
                [data count_read] = BZ2_fread(fid_info.fid, count, precision, skip, machineformat);
        end;
    case 'GZ_'
        switch nargin
            case 3
                [data count_read] = GZ_fread(fid_info.fid, count, precision);
            case 2
                [data count_read] = GZ_fread(fid_info.fid, count);
            case 1
                [data count_read] = GZ_fread(fid_info.fid);
            case 4
                [data count_read] = GZ_fread(fid_info.fid, count, precision, skip);
            otherwise
                [data count_read] = GZ_fread(fid_info.fid, count, precision, skip, machineformat);
        end;
    otherwise
        data = [];
        count_read = 0;
end;