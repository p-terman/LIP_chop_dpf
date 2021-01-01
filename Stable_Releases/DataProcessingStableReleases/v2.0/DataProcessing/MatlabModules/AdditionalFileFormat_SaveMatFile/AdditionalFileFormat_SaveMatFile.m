function status = AdditionalFileFormat_SaveMatFile(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)
%
% Simple module to read .rq file and save it as a rqm.mat file which is
% easily loaded into matlab. Run at end of DP framework when all RQs are
% calculated.
%
% status = AdditionalFileFormat_SaveMatFile(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)
%
%
%
% Required RQs: some.
%
%
% Versioning:
%       2013-04-11 JJC - copied from similar code in Run02 style
%       processing.
%       2013-04-11 JRV - Renamed to conform with module naming scheme
%       2013-05-14 JRV - Now calles them .rq.mat files for consistency
%
% RQ versions:
%
%
%% Load .rq file

status = [];

myname = 'SaveRQFile2MATFile';
fprintf('\n\n *** Starting module %s\n',myname);

dp = LUXLoadRQ1s_framework(filename_rq, data_path_rq);
dp.livetime_latch_samples = dp.admin.livetime.livetime_latch_samples;
dp.livetime_end_samples = dp.admin.livetime.livetime_end_samples;
%% Check for matfiles directory in rq directory
matfiles_dir = dir([data_path_rq '/matfiles']);
if isempty(matfiles_dir)
    mkdir(data_path_rq, '/matfiles/');
end

%% Write .mat file
LUXWriteRQMFile_framework(filename_rq, [data_path_rq '/matfiles/'], dp, dp.admin);
fprintf('saving mat file %s.mat\n',filename_rq);


fprintf('Done!\n')

