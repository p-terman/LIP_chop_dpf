#!/usr/bin/python
'''
auto_run_dp_framework_oscar.py

This script is intended for production DP running on the Oscar cluster.
This script is intended to be submitted to each core in parallel running on Oscar.

Arguements:
	evt_directory - name of evt directory (e.x.: lux10_20130205T1154_eb00006)
	cp - cp analysis ID
    job_num - this thread's number in the job array
	number_of_cores - total number of cores for the entire job array

2012-12-10 - JRV - Created
2013-01-24 - JRV - Finds your home directory and uses that by default
'''

# call script to set python paths
import os
home_dir = os.path.expanduser("~")
execfile(home_dir + '/LUXcode/Stable_Releases/DataProcessingStableReleases/v2.0/DataProcessing/DataProcessingFramework/Utilities/LUXSetPythonPaths.py')

# import required modules
import sys
import re
import math
import LUXGenerateRQsPyMod
from subprocess import Popen
from LUXModuleWrapperPyMod import GetGlobalSettings
from ReportLUXcodePathPyMod import ReportLUXcodePath
from glob import glob

# get input arguements
evt_directory = sys.argv[1]
cp = int(sys.argv[2])
job_num = int(sys.argv[3])
number_of_cores = int(sys.argv[4])
slurm_procid = int(sys.argv[5])

filename_prefix = re.match('^(lux.._[0-9]+T[0-9]+)', evt_directory).group()

# set directories
evt_dir = home_dir + '/scratch/evt/' + evt_directory + '/'
rq_dir = home_dir + '/scratch/rq/' + filename_prefix + '_cp%0.5d/' % cp

# find all .evt files in evt_dir and throw them into a list
evt_file_list = list()
evt_files = glob(evt_dir + '/*evt')
evtgz_files = glob(evt_dir + '/*evt.gz')

for f in evt_files:
    evt_file_list.append(os.path.basename(f))
for f in evtgz_files:
    if f.rstrip('.gz') not in evt_file_list:
        evt_file_list.append(os.path.basename(f))
evt_file_list.sort()

# figure out which files this particular thread is responsible for processing
# this was ported from auto_run_analysis.m - this probably can be done more
# pythonically

tot_num_files = len(evt_file_list)
files_per_job = math.floor(float(tot_num_files) / number_of_cores)
remaining_files = int(tot_num_files -  number_of_cores * files_per_job)
remaining_file_nums = range(tot_num_files - remaining_files + 1, tot_num_files + 1)

start_idx = int(job_num * files_per_job)
data_set_list = range(start_idx, int(start_idx + files_per_job))

if (job_num < remaining_files):
	data_set_list.append(remaining_file_nums[job_num] - 1)

evt_file_list_subset = [evt_file_list[ii] for ii in data_set_list]

# run the data processing framework
output_logs = LUXGenerateRQsPyMod.LUXGenerateRQs(evt_dir, rq_dir, cp, evt_file_list_subset)
