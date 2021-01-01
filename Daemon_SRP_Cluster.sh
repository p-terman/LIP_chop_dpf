#!/bin/bash


primaryFolder="/users/pterman/scratch/"${1}""

mkdir "$primaryFolder"
mkdir ""${primaryFolder}"/evt"
mkdir ""${primaryFolder}"/matfiles"

#ssh to gsk 60 to look for files and send. 
ssh lux_mirror_user@gsk-60.het.brown.edu "./find_and_send_files.sh "${1}""

gunzip -f "${primaryFolder}"/evt/*.evt.gz

# how many files need to be reprocessed
numFiles=$(ls -1q /users/pterman/scratch/"${1}"/matfiles/*.rq.mat | wc -l)

# have 100 files processed in each array, +1 is for remainder
numArray=$(expr "$numFiles" / 100)
numArray=$(expr "$numArray" + 1)

#make BASH_RUN_SRP

dollar='$'
quote="'"
doublequote='"'
chip='`'
rm BASH_SRP_CLUSTER

echo -e "#!/bin/bash">> "BASH_SRP_CLUSTER"
echo -e "#">> "BASH_SRP_CLUSTER"
echo -e "# New job submission scheme for CCV queue LIP select reprocessing">> "BASH_SRP_CLUSTER"
echo -e "#">> "BASH_SRP_CLUSTER"
echo -e "# 20130114 - JRV - Modified from submit_batch_job to run on new Oscar cluster using SLURM">> "BASH_SRP_CLUSTER"
echo -e "# 20190402 PAT rewrote old submission file for LIP Reprocessing">> "BASH_SRP_CLUSTER"
echo -e "">> "BASH_SRP_CLUSTER"
echo -e "################################################################">> "BASH_SRP_CLUSTER"
echo -e "# YOU SHOULD NOT EDIT THIS FILE AT ALL!!!                      #">> "BASH_SRP_CLUSTER"
echo -e "# There is no more need to modify or duplicate run_analysis.m  #">> "BASH_SRP_CLUSTER"
echo -e "################################################################">> "BASH_SRP_CLUSTER"
echo -e "">> "BASH_SRP_CLUSTER"
echo -e "# name your job to help you and others identify it in the queue">> "BASH_SRP_CLUSTER"
echo -e "#SBATCH --job-name=SRP${1##*_cp}">> "BASH_SRP_CLUSTER"
echo -e "">> "BASH_SRP_CLUSTER"
echo -e "# set mail address">> "BASH_SRP_CLUSTER"
echo -e "#SBATCH --mail-user=terman.paul@yahoo.com">> "BASH_SRP_CLUSTER"
echo -e "#SBATCH --mail-type=END">> "BASH_SRP_CLUSTER"
echo -e "">> "BASH_SRP_CLUSTER"
echo -e "# estimate the time needed by your job, as HH:MM:SS">> "BASH_SRP_CLUSTER"
echo -e "#SBATCH --time=1:30:00">> "BASH_SRP_CLUSTER"
echo -e "">> "BASH_SRP_CLUSTER"
echo -e "# set which partition">> "BASH_SRP_CLUSTER"
echo -e "#SBATCH --partition=batch">> "BASH_SRP_CLUSTER"
echo -e "">> "BASH_SRP_CLUSTER"
echo -e "## set priority">> "BASH_SRP_CLUSTER"
echo -e "">> "BASH_SRP_CLUSTER"
echo -e "# set cores per job">> "BASH_SRP_CLUSTER"
echo -e "#SBATCH --ntasks=1">> "BASH_SRP_CLUSTER"
echo -e "#SBATCH --ntasks-per-node=1">> "BASH_SRP_CLUSTER"
echo -e "">> "BASH_SRP_CLUSTER"
echo -e "# set memory per cpu">> "BASH_SRP_CLUSTER"
echo -e "#SBATCH --mem-per-cpu=5G">> "BASH_SRP_CLUSTER"
echo -e "">> "BASH_SRP_CLUSTER"
echo -e "# constrain to intel nodes only">> "BASH_SRP_CLUSTER"
echo -e "#SBATCH --constraint=intel">> "BASH_SRP_CLUSTER"
echo -e "">> "BASH_SRP_CLUSTER"
echo -e "# set number of jobs in array">> "BASH_SRP_CLUSTER"
echo -e "#SBATCH --array=1-${numArray}">> "BASH_SRP_CLUSTER"
echo -e "">> "BASH_SRP_CLUSTER"
echo -e "# [optional] combine stdout and stderr streams into one output file">> "BASH_SRP_CLUSTER"
echo -e "">> "BASH_SRP_CLUSTER"
echo -e "# [optional] specify a different output file than the default JobName.out">> "BASH_SRP_CLUSTER"
echo -e "#SBATCH -o MySerialJob-%j.out">> "BASH_SRP_CLUSTER"
echo -e "#SBATCH -e MySerialJob-%j.out">> "BASH_SRP_CLUSTER"
echo -e "# set path to auto run file">> "BASH_SRP_CLUSTER"
echo -e "PATH_TO_AUTO_RUN=/users/pterman/run_select_dpf_brown_Cluster.py">> "BASH_SRP_CLUSTER"
echo -e "">> "BASH_SRP_CLUSTER"
echo -e "# set CP">> "BASH_SRP_CLUSTER"
echo -e "CP=${dollar}{1##*_cp}">> "BASH_SRP_CLUSTER"
echo -e "">> "BASH_SRP_CLUSTER"
echo -e "# set the filename prefix">> "BASH_SRP_CLUSTER"
echo -e "Filename_prefix=${dollar}{1%%_cp*}">> "BASH_SRP_CLUSTER"
echo -e "EVT_DIRECTORY=/users/pterman/scratch/${dollar}{1}/evt/">> "BASH_SRP_CLUSTER"
echo -e "RQ_DIRECTORY=/users/pterman/scratch/${dollar}{1}/matfiles/">> "BASH_SRP_CLUSTER"
echo -e "">> "BASH_SRP_CLUSTER"
echo -e "# set other paths">> "BASH_SRP_CLUSTER"
echo -e "LUXcode_path=/users/pterman/LUXCode/Stable_Releases/DataProcessingStableReleases/v2.0/DataProcessing/">> "BASH_SRP_CLUSTER"
echo -e "dp_settings=/users/pterman/data_processing_settings_run03_oscar.xml">> "BASH_SRP_CLUSTER"
echo -e "">> "BASH_SRP_CLUSTER"
echo -e "# set the number of cores (NODES * ppn): e.g. NODES = 0-3 & ppn=4 =>">> "BASH_SRP_CLUSTER"
echo -e "NUM_CORES=1">> "BASH_SRP_CLUSTER"
echo -e " ">> "BASH_SRP_CLUSTER"
echo -e "#######################################################">> "BASH_SRP_CLUSTER"
echo -e "# YOU DO NOT NEED TO CHANGE ANYTHING BELOW THIS POINT #">> "BASH_SRP_CLUSTER"
echo -e "#######################################################">> "BASH_SRP_CLUSTER"
echo -e " ">> "BASH_SRP_CLUSTER"
echo -e "srun bash -l -c ${quote}${quote}${quote}">> "BASH_SRP_CLUSTER"
echo -e "    cd ${dollar}PWD">> "BASH_SRP_CLUSTER"
echo -e "    export ID=${dollar}((SLURM_TASKS_PER_NODE*(SLURM_ARRAY_TASK_ID-1)+SLURM_PROCID));">> "BASH_SRP_CLUSTER"
echo -e "    echo ${doublequote}Starting job ${dollar}ID on ${dollar}HOSTNAME (CPU ${dollar}SLURM_PROCID)${doublequote};">> "BASH_SRP_CLUSTER"
echo -e " ">> "BASH_SRP_CLUSTER"
echo -e "    export PATH_TO_AUTO_RUN=${chip}echo -e ${doublequote}${dollar}0${doublequote}${chip}">> "BASH_SRP_CLUSTER"
echo -e "    export Filename_prefix=${chip}echo -e ${doublequote}${dollar}1${doublequote}${chip}">> "BASH_SRP_CLUSTER"
echo -e "    export EVT_DIRECTORY=${chip}echo -e ${doublequote}${dollar}2${doublequote}${chip}">> "BASH_SRP_CLUSTER"
echo -e "    export RQ_DIRECTORY=${chip}echo -e ${doublequote}${dollar}3${doublequote}${chip}">> "BASH_SRP_CLUSTER"
echo -e "    export LUXcode_path=${chip}echo -e ${doublequote}${dollar}4${doublequote}${chip}">> "BASH_SRP_CLUSTER"
echo -e "    export dp_settings=${chip}echo -e ${doublequote}${dollar}5${doublequote}${chip}">> "BASH_SRP_CLUSTER"
echo -e "    export CP=${chip}echo -e ${doublequote}${dollar}6${doublequote}${chip}">> "BASH_SRP_CLUSTER"
echo -e "    export LD_LIBRARY_PATH=${dollar}LD_LIBRARY_PATH:/users/jverbus/mysql/">> "BASH_SRP_CLUSTER"
echo -e " ">> "BASH_SRP_CLUSTER"
echo -e "    python ${dollar}PATH_TO_AUTO_RUN ${dollar}Filename_prefix ${dollar}EVT_DIRECTORY ${dollar}RQ_DIRECTORY ${dollar}LUXcode_path ${dollar}dp_settings ${dollar}CP ${dollar}SLURM_ARRAY_TASK_ID">> "BASH_SRP_CLUSTER"
echo -e "${quote}${quote}${quote} ${doublequote}${dollar}PATH_TO_AUTO_RUN${doublequote} ${doublequote}${dollar}Filename_prefix${doublequote} ${doublequote}${dollar}EVT_DIRECTORY${doublequote} ${doublequote}${dollar}RQ_DIRECTORY${doublequote} ${doublequote}${dollar}LUXcode_path${doublequote} ${doublequote}${dollar}dp_settings${doublequote} ${doublequote}${dollar}CP${doublequote} ${doublequote}${dollar}SLURM_ARRAY_TASK_ID${doublequote} ">> "BASH_SRP_CLUSTER"

sbatch BASH_SRP_CLUSTER "${1}"
