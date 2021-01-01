#!/bin/bash


primaryFolder="/users/pterman/scratch/"${1}""

mkdir "$primaryFolder"
mkdir ""${primaryFolder}"/evt"
mkdir ""${primaryFolder}"/matfiles"

#ssh to gsk 60 to look for files and send. 
ssh lux_mirror_user@gsk-60.het.brown.edu "./find_and_send_files.sh "${1}""

gunzip -f  "${primaryFolder}"/evt/*.evt.gz

sbatch BASH_RUN_SRP "${1}"
