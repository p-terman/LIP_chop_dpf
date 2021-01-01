% function status = dual_fseek(fid_info, offset, origin)
%
% Calls the appropriate fseek, given by fid_info.filecommand_prefix
% CED
% Release Version 1.0

function status = dual_fseek(fid_info, offset, origin)

%% get fid info
fid = fid_info.fid;
prefix = fid_info.filecommand_prefix;

%% make call
switch prefix
    case ''
        switch nargin
            case 1
                status = fseek(fid);
            case 2
                status = fseek(fid, offset);
            otherwise
                status = fseek(fid, offset, origin);
        end;
    case 'BZ2_'
        switch nargin
            case 1
                status = BZ2_fseek(fid);
            case 2
                status = BZ2_fseek(fid, offset);
            otherwise
                status = BZ2_fseek(fid, offset, origin);
        end;
    case 'GZ_'
        switch nargin
            case 1
                status = GZ_fseek(fid);
            case 2
                status = GZ_fseek(fid, offset);
            otherwise
                status = GZ_fseek(fid, offset, origin);
        end;
    case 'MEM_'
        switch nargin
            case 1
                status = MEM_fseek(fid);
            case 2
                status = MEM_fseek(fid, offset);
            otherwise
                status = MEM_fseek(fid, offset, origin);
        end;
    otherwise
        status = -1;
end;

