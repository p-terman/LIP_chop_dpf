function LUXWriteRQMFile(output_filename,output_path,dp,headerinfo)
%
% LUXWriteRQMFile(output_filename,output_path,dp,headerinfo)
%
% This function saves an RQM file (no dp structure, so you can pick
% individual variables to load)
%
% INPUTS:
%       output_filename - filename to write (same as input to LUXBinaryRQWriter.m)
%        output_path - location to save
%                 dp - structure with first pass variables
%         headerinfo - admin settings
%
% Versioning:
%   20111111 CHF - Created
%   20111115 PS  - takes input 'output_filename' rather than 'filename_rqm', for consistency
%                  with LUXBinaryRQWriter.m
%   20130306 JJC - now writes rqm from .rq (instead of .rq1). Small change.
%   20130515 JRV - Now writes .rq.mat instead of .rqm.mat
%%

clear varname

if isfield(dp,'admin')
    dp = rmfield(dp,'admin');
end

temp_var_names = nstruct2cell(dp);

for nn = 1:length(temp_var_names)
    dotinds = strfind(temp_var_names{nn},'.');
    varname{nn} = temp_var_names{nn}(dotinds(end)+1:end);
    eval(sprintf('%s = %s;',varname{nn},temp_var_names{nn}));
end

admin = headerinfo;
varname{end+1} = 'admin';

% Save .rqm file
%filename_rqm = strrep(output_filename,'.rq1','.rqm');
%filename_rqm = strrep(output_filename,'.rq','.rqm'); % now does nothing - for consistency
filename_rqm = output_filename;
save([output_path filesep filename_rqm '.mat'],varname{:});
fprintf('\nwrote file:\n %s\n',[output_path filesep filename_rqm '.mat']); drawnow;
