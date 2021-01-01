% This script should be called by your startup.m file to load the LUX
% matlab libraries.
%
% This will add folders containing LUX code to your PATH, and will also
% compile any MEX functions that have either not been complied or have
% changed since their last compile.  (Note:  compiled code is not stored in
% the SVN, only in your working copy, so if you update the SVN to a version
% with new or changed mex files have been changed, this script will
% re-compile those mex files).
%
% 02/12/08, CED
% 2010-11-08 JJC now just add Trunk by default

thisfolder=fileparts(which('LoadLuxMatlabLibraries'));

%pathbase = [thisfolder filesep '..' filesep '..' filesep];
pathbase = '~/matlab/library/LUXcode/Trunk/';

warning off;

%% Compile all MEX functions
CompileMexFiles(pathbase);
sprintf('Mex files compiled\n');

%% Add libraries to path
addpath(genpath(pathbase));

disp('LUX Libraries added to Path');

%% finish up
clear pathbase
warning on;

