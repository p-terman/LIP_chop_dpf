% function message = dual_ferror(fid_info)
%
% calls the appropriate ferror function
%
% CED
% Release Version 1.0

function message = dual_ferror(fid_info)

%% get fid info
fid = fid_info.fid;
prefix = fid_info.filecommand_prefix;

%% make call
switch prefix
    case ''
        message = ferror(fid);
    case 'BZ2_'
        message = BZ2_ferror(fid);
    case 'GZ_'
        message = GZ_ferror(fid);
    case 'MEM_'
        message = MEM_ferror(fid);
    otherwise
        message = 'Unknown filecommand_prefix';
end;
