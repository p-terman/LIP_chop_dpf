function status = load_saved_matlab_path(path_to_save_file)
% This function loads a matlab path from a .mat file
% function status = load_saved_matlab_path(path_to_save_file)
%
%  Inputs:
%   path_to_save_file - path to saved .mat file (including filename)
%
% Outputs:
%              status - status (either current path or empty
%
% 2012-11-20 - JRV - Created

%%

status = [];

try
    a = load(path_to_save_file);
    path(a.str);
    status = a.str;
catch
    fprintf('ERROR: Could not load path\n');
end

