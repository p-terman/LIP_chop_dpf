% function [status] = CompileMexFiles(pathbase);
%
% Compiles all the LUX Analysis Mex Functions.
%
% Input: 
% pathbase: the top level of the LUXcode working copy
%
% Output:
% status:  =0 for successful run, <0 for unsuccesful run, >0 for warnings
%
% 02/12/08, CED

function [status] = CompileMexFiles(pathbase)

%% set defaults
if nargin<1 || ~ischar(pathbase) || ~exist(pathbase,'dir')
    status = -1;
    return
end

%% set mex path list
mexdirs = { ...
    'FirstPassAnalysis/MatlabFirstPass/mex_files', ...
    'Utilities/MatlabUtilities/tcp_udp_ip', ...
    'Utilities/MatlabUtilities/GeneralFileHandling/DualRead', ...
    'Utilities/MatlabUtilities/GeneralFileHandling/GZRead', ...
    'Utilities/MatlabUtilities/GeneralFileHandling/BZ2Read'};

%% save current directory
starting_directory = pwd;

%% loop over mex paths, and compile any uncompiled or changed mex functions
for d=1:length(mexdirs)
    cd([pathbase filesep mexdirs{d}])
    clist = dir(pwd);
    for f=1:length(clist)
        % all files that can be mexed should have a .c extension
        if length(clist(f).name)>2 && strcmp(clist(f).name(end-1:end),'.c')
            % determine which file matlab assosciates with this name-  this
            % will be either a .m file or a compiled mex file with the
            % local architecture extension
            %
            % NOTE -- if no .m file or compiled mex file exists, this file
            % will not be compiled.  Every mex (.c) file MUST have an
            % assosciated .m file -- this .m file contains the help comment
            % for the mex file, and should issue a comment when run along
            % the lines of 'You need to compile filename.c'
            foundfile = which(clist(f).name(1:end-2));
            if length(foundfile)>2
                finfo = dir(foundfile);
                % if the found file is a .m file, then the mex file has not
                % been compiled -- compile it.  If the found file is older
                % than the mex file, recompile it.
                if strcmp(foundfile(end-1:end),'.m') || (datenum(finfo.date)<datenum(clist(f).date))
                    try
                        % some mex files require special flags when
                        % compiling
                        switch clist(f).name
                            case 'BZ2Read.c'
                                if ~ispc
                                    mex -l bz2 BZ2Read.c;
                                    disp('mexed BZ2Read.c');
                                else
                                    disp('Cannot mex BZ2Read.c on this platform -- hopefully the BZ2Read.mexw32 file in the repository will work for you.');
                                end
                            case 'GZRead.c'
                                if ~ispc
                                    mex -l z GZRead.c;
                                    disp('mexed GZRead.c');
                                else
                                    disp('Cannot mex GZRead.c on this platform -- hopefully the GZRead.mexw32 file in the repository will work for you.');
                                end
                            case 'pnet.c'
                                if ~ispc
                                    mex -O pnet.c
                                else
                                    matlab_install_dir = ...
                                        fileparts(fileparts(fileparts(fileparts(which('mex')))));
                                    mex('-O', 'pnet.c', ...
                                        [matlab_install_dir '\sys\lcc\lib\wsock32.lib'], '-DWIN32');
                                end
                            otherwise
                                mex(clist(f).name);
                                disp(['mexed ' clist(f).name]);
                        end
                    catch
                        disp(sprintf('********************\nFailed to mex %s\n****************',clist(f).name));
                    end
                end
            end
        end
    end
end

%% finish up
cd(starting_directory);
status=0;
