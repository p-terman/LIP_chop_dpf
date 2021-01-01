function status = create_save_matlab_path(path_to_save_file, dp_path)
% This function creates the DataProcessing path and saves it to a file
% function status = create_save_matlab_path(path_to_save_file, dp_path)
%
%  Inputs:
%   path_to_save_file - path to save the .mat file (including filename)
%             dp_path - path to the DataProcessing directory
%
% Outputs:
%              status - status (either path to saved file or empty)
%
% 2012-11-20 - JRV - Created

%%

status = [];

try
    disp(dp_path)
    disp(path_to_save_file)
    addpath(genpath([dp_path]));
    str = path;
    save(path_to_save_file,'str');
    status = path_to_save_file;
catch
    fprintf('ERROR: Could not save path\n');
end

