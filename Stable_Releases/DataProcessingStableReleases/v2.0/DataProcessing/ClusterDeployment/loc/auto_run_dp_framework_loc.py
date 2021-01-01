#!/usr/bin/python
'''
auto_run_dp_framework_loc.py

This script is intended for production DP running on the Oscar cluster.
This script is intended to be submitted to each core in parallel running on Oscar.

Arguements:
  evt_directory - name of evt directory (e.x.: lux10_20130205T1154_eb00006)
	cp - cp analysis ID
    job_num - this thread's number in the job array
	number_of_cores - total number of cores for the entire job array

2013-02-05 - MW - Initilial Submission (Stolen)
2013-02-11 - MW - Update to use eb files.
2014-03-31 - AL - Moved code from LUXGenerateRQs into here so that the XML settings 
files for LUG IQs and DP settings are now created only once (in the head node), when 
submitting the jobs. The location of these temporary files (which are created in
home_dir + '/LOC_DataProcessing/SavedDatasets/xml/' followng the path convention for
evts and rqs) is passed on to LUXGenerateRQs. This change should allow clusters without 
node internet connection to be able to run the DP framework.
2014-07-15 - AL - Reverting to the old style of getting XML settings: each node does
its own communication with the LUG

'''

# call script to set python paths
import os
home_dir = os.path.expanduser("~")
execfile('/home/lux_admin/LOC_DataProcessing/SpecificRevisions/4109/DataProcessingFramework/Utilities/LUXSetPythonPaths.py')

# import required modules
import sys
import math
import re
import LUXGenerateRQsPyMod
from subprocess import Popen
from LUXModuleWrapperPyMod import GetGlobalSettings
from ReportLUXcodePathPyMod import ReportLUXcodePath

# get input arguements
#filename_prefix = sys.argv[1]
evt_directory = sys.argv[1]
cp = int(sys.argv[2])
job_num = int(sys.argv[3])
number_of_cores = int(sys.argv[4])
global_job_id = int(sys.argv[5])

# Determine the lux10_blahblah_eb_blahblah filename prefix.
filename_prefix = re.match('^(lux10_[0-9]+T[0-9]+)', evt_directory).group()

# set directories
# Not yet...
evt_dir = '/home/lux_admin/LOC_DataProcessing/SavedDatasets/evt/' + evt_directory + '/'
rq_dir = '/home/lux_admin/LOC_DataProcessing/SavedDatasets/rq/' + filename_prefix + '_cp%0.5d/' % cp
xml_dir = home_dir + '/LOC_DataProcessing/SavedDatasets/xml/'
#evt_dir = '/home/lux_admin/working_space/develop_job_submission/evt/' + filename_prefix + '/'
#rq_dir = '/home/lux_admin/working_space/develop_job_submission/rq/' + filename_prefix + '_cp%0.5d/' % cp


# find all .evt files in evt_dir and throw them into a list
evt_file_list = list()
files = os.listdir(evt_dir)
for f in files:
    if (f.startswith('lux') and f.endswith('.evt')):
        evt_file_list.append(f)
    if (f.startswith('lux') and f.endswith('.evt.gz')):
        evt_file_list.append(f)
evt_file_list.sort()

# figure out which files this particular thread is responsible for processing
# this was ported from auto_run_analysis.m - this probably can be done more
# pythonically

tot_num_files = len(evt_file_list)
files_per_job = math.floor(float(tot_num_files) / number_of_cores)
remaining_files = int(tot_num_files -  number_of_cores * files_per_job)
remaining_file_nums = range(tot_num_files - remaining_files + 1, tot_num_files + 1)

start_idx = int((job_num-1) * files_per_job)
data_set_list = range(start_idx, int(start_idx + files_per_job))
#print "number_of_cores:     ", number_of_cores
#print "tot_num_files:       ", tot_num_files
#print "files_per_job:       ", files_per_job
#print "remaining_files:     ", remaining_files
#print "remaining_file_nums: ", remaining_file_nums
#print "job_num:             ", job_num
#print "start_idx:           ", start_idx
#print "data_set_list:       ", data_set_list

# Print out the files we're using.
mfirst_file = (job_num - 1)*files_per_job
mlast_file = mfirst_file + files_per_job
print "I should at least do files", mfirst_file, "through", mlast_file
if (job_num-1 < remaining_files):
  data_set_list.append(remaining_file_nums[job_num-1] - 1)
  print "\tin addition I'm taking file", remaining_file_nums[job_num-1] - 1

evt_file_list_subset = [evt_file_list[ii] for ii in data_set_list]

print ""
sys.stdout.flush()

# run the data processing framework
import time
start_time = time.time()
output_logs = LUXGenerateRQsPyMod.LUXGenerateRQs(evt_dir, rq_dir, cp, evt_file_list_subset)
end_time = time.time()

print "RUN TIME FOR", len(evt_file_list_subset), "FILES IS", end_time - start_time
