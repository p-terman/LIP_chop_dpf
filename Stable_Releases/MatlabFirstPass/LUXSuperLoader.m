function [settings, fid_out] = LUXSuperLoader(file_name, file_path, options)
% [settings, fid] = LUXSuperLoader(file_name, file_path, options)
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



fid = fopen(file_full_name,'rb','l');
fread(fid,1,'uint32');


%% Read XML settings and do LUG query

size_of_xml = fread(fid, 1, 'uint32');
xml_string = char(fread(fid, size_of_xml, 'char'))';

if load_settings
    
    settings = XMLParser(char(xml_string));
    
    % Do LUG query if requested and file is the proper type
    if query_lug
        
        % Artificial insemination for testing purposes
        settings.daq_settings.filename = file_name;
        
        filename_prefix = strtok(settings.daq_settings.filename, 'f');
        filename_prefix = filename_prefix(1:end-1);
        
        lug_out = LUX01LoadState(filename_prefix);
        
        % Need to fix this on the LUX01LoadState end!
        if ~isfinite(lug_out.amplification.preamp)
            lug_out.amplification.preamp = 1;
        end
        if ~isfinite(lug_out.amplification.postamp.out1_height)
            lug_out.amplification.postamp.out1_height = 1;
        end
        if ~isfinite(lug_out.amplification.postamp.out1_area)
            lug_out.amplification.postamp.out1_area = 1;
        end
        for ii=1:length(lug_out.ch)
            if ~isfinite(lug_out.ch(ii).sphe_area_mVns)
                error('Dataset %s - LUG calibration for PMT %d has a bad value', ...
                    filename_prefix, ii);
            end
        end
        
        settings.lug_settings = lug_out;
        
    end
    
end

%% Assign output and close file

if nargout > 1
    fid_out = fid;
else
    fclose(fid);
end


end

