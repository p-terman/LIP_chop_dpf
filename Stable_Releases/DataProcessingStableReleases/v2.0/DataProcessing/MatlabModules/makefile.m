%%

% assume in home for now...
[flag,hdir] = unix(['echo $HOME']);
hdir = hdir(1:end-1); % Remove LF at end
dp_dir = [hdir '/LUXcode/Trunk/DataProcessing/']; 
addpath(genpath(dp_dir));
matlab_modules_dir = [dp_dir '/MatlabModules/'];


matlab_bin_dir = [matlab_modules_dir '/bin/'];
if ~exist(matlab_bin_dir,'dir')
    mkdir(matlab_bin_dir)
end

modules_to_compile = dir([matlab_modules_dir '/*_*']);

%%

for ii = 1:length(modules_to_compile)
    sdir = [matlab_modules_dir '/' modules_to_compile(ii).name '/'];
    sfile = [modules_to_compile(ii).name '.m'];
    cd(sdir);
    
    bdir = [matlab_bin_dir '/' modules_to_compile(ii).name '/'];
    try
        rmdir(bdir,'s');
    catch
    end
    
    mkdir(bdir);
    
    fprintf('\n\nCreating directory: %s\n\n',bdir);
    try
        eval(sprintf('mcc -R -nojvm -R -nodesktop -R -nosplash -m ./%s -d %s', sfile, bdir));
    catch
        fprintf('\n\nERROR: Compilation of %s failed; deleting %s\n\n', sfile, bdir);
        rmdir(bdir,'s');
    end
end
