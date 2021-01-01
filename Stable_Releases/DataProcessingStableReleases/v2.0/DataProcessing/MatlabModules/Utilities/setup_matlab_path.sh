#!/bin/bash
#
# This script writes the file matlab_path.mat containing your Matlab path
# for the DP framework. It should be run every time you change your
# DataProcessing SVN.
#
# Inputs:
#   $1 - The path to your DataProcessing installation
#   $2 - The path to your matlab executable

dp_path=$1
matlab_path=$2

path_to_save_file="$dp_path/MatlabModules/Utilities/matlab_path.mat"

cd $dp_path/MatlabModules/Utilities/

run_command="$matlab_path -nojvm -nodesktop -nosplash -r create_save_matlab_path('$path_to_save_file','$dp_path');quit"

$run_command
