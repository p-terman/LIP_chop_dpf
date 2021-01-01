% function position = dual_ftell(fid_info)
%
% calls the appropriate ftell function,
%
% CED
% Release Version 1.0

function position = dual_ftell(fid_info)

%% get fid info
fid = fid_info.fid;
prefix = fid_info.filecommand_prefix;

%% make call
switch prefix
    case ''
        position = ftell(fid);
    case 'BZ2_'
        position = BZ2_ftell(fid);
    case 'GZ_'
        position = GZ_ftell(fid);
    case 'MEM_'
        position = MEM_ftell(fid);
    otherwise
        position = -1;
end;