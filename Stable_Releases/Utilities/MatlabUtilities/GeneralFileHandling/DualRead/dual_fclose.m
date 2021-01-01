% function status = dual_fclose(fid_info)
%
% If passed an fid, calls the appropriate fclose function,
% if passed the string 'all', calls all fclose functions with argument
% 'all'.
%
% CED
% Release Version 1.0

function status = dual_fclose(fid_info)

%% check input
if isstruct(fid_info)
%% get fid info
    fid = fid_info.fid;
    prefix = fid_info.filecommand_prefix;

%% make call
    switch prefix
        case ''
            status = fclose(fid);
        case 'BZ2_'
            status = BZ2_fclose(fid);
        case 'GZ_'
            status = GZ_fclose(fid);
        case 'MEM_'
            status = MEM_fclose(fid);
        otherwise
            status = -1;
    end;
%% otherwise get all fid lists
elseif ischar(fid_info) && strcmp(fid_info,'all')
    status = 0;
    status = status + fclose('all');
    status = status + BZ2_fclose('all');
    status = status + GZ_fclose('all');
    status = status + MEM_fclose('all');
else
    disp('Invalid input to dual_fclose');
end;
