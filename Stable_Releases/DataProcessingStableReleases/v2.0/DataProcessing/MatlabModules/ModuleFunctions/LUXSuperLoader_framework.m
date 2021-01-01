function [settings, fid_out offset] = LUXSuperLoader_framework(file_name, file_path, options)
% [settings, fid] = LUXSuperLoader_framework(file_name, file_path, options)
% 
% This function creates a structure "settings" containing information about the
% dataset to be analyzed;  you normally need to run this before looking at a dataset
%
% Inputs:
%  file_name - name of file to load
%  file_path - path (relative or absolute) to file. Leave blank to search
%   current directory
%  options - structure with the following recognized fields:
%   .load_settings - load the settings from the file prefix and return in
%     the settings structure. Set to 0 to speed up execution, if dataset 
%     settings have already been loaded. [Default: 1]
%   .query_lug - do a LUG query and append settings to settings structure.
%     NOTE: LUG will not be queried if .load_settings = 0.
%     [Default: 0]
%   .use_xml_file - use xml file for loading data
% 
% Outputs:
%  settings - XML settings structure loaded from file
%  fid - file pointer of type used by f* functions (i.e. fread)
%   pointing to beginning of data in file. File will be closed if fid is
%   not assigned at output
% 
% 2010-04-15 DCM v1.0
% 

%% Defaults

load_settings = true;
query_lug = false;

settings = [];
fid = [];


%% Input checking

if nargin < 1
    disp('Please specify the file to load');
    return
end

if nargin < 2 || isempty(file_path)
    file_path = pwd;
end

% Evaluate options structure
if nargin >= 3
    field_names = fieldnames(options);
    for ii=1:length(field_names)
        eval([field_names{ii} ' = options.' field_names{ii} ';']);
    end
end

%% Open file

file_full_name = [file_path filesep file_name];

if ~exist(file_full_name,'file')
    % Didn't find this exact name -- do a search
    files_list = dir([file_full_name '*']);
    if ~isempty(files_list)
        file_full_name = [file_path filesep files_list(1).name];
    else
        error('Could not locate file name matching "%s" at %s', ...
            file_name, file_path);
        return
    end
end


%file_full_name
fid = fopen(file_full_name,'rb','l');
fread(fid,1,'uint32');


%% Read XML settings and do LUG query

size_of_xml = fread(fid, 1, 'uint32');
xml_string = char(fread(fid, size_of_xml, 'char'))';

if load_settings
    settings = XMLParser_framework(char(xml_string));
end

%% Assign output and close file
offset = ftell(fid);
if nargout > 1
    fid_out = fid;
else
    fclose(fid);
end


