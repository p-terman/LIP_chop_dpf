#!/usr/bin/python
"""
DataProcessingDaemonLoc.py

This script runs in a screen on the loc head node.
In v1.0 of the DPF it looks for new datasets, creates
the submit_batch_job file, and calls qsub to submit the
job.

In v2.0 of the DPF it will be expanded to do a number
of LUG queries for processing requests/settings.

2012-12-05 - MW  - Created (stolen)
"""

# Set Python Paths
import os
import glob

execfile('/home/lux_admin/LOC_DataProcessing/SpecificRevisions/4109/DataProcessingFramework/Utilities/LUXSetPythonPaths.py')

# import modules
import datetime
import time
import re
from subprocess import Popen, PIPE
from LUXLUGQueriesPyMod import LUXPrepareDPSettings
from ReportLUXcodePathPyMod import ReportLUXcodePath, ReportDataProcessingPath

# Get the LUXcode directory
luxcode_path = ReportLUXcodePath()
dp_path = ReportDataProcessingPath()

# define things

EVT_DIR = '/home/lux_admin/LOC_DataProcessing/SavedDatasets/evt'
LOC_DP_STARTED_FLAG = 'loc_dp_started'
LOC_DP_ERROR_FLAG = 'loc_dp_error'
LOC_SYNC_DONE_FLAG = 'loc_sync_done'
BASE_LOG_DIR = "/share/data/logs/dplogs/"
DEFAULT_DP_SETTINGS_FILE = '/home/lux_admin/LOC_DataProcessing/SpecificRevisions/4109/PACMAN_settings.xml'
PRIORITY = 'priority'

def create_sarray_settings_file(evt_directory, cp, evt_path, priority):
    """
    This function creates the sarray settings file. It is the equivalent
    to create_submit_batch_job.sh in the Run02 analysis framework.
    """

    #return
    # import required modules

    import sys
    import os
    import datetime
    home_dir = os.path.expanduser("~")

    # make sure evt_path is a directory and has loc_sync_done flag in it
    if not os.path.isdir(evt_path + '/' + evt_directory):
        sys.exit("ERROR: Filename prefix doesn't correspond to a directory in the evt path")

    if not evt_directory.startswith('lux'):
        sys.exit("ERROR: Input must be a valid filename prefix")

    if not os.path.isfile( (evt_path + '/' + evt_directory + '/' + 'loc_sync_done') ):
    	sys.exit("ERROR: LOC sync has not completed! Exiting...")

    # Determine the lux10_blahblah from lux10_blahblah_ebblahblah
    filename_prefix = re.match('^(lux10_[0-9]+T[0-9]+)', evt_directory).group() 
    print filename_prefix, cp, evt_path, priority

    # get list of evt files
    #zipped_evt_files = glob.glob((evt_path + '/' + evt_directory + '/' + '*.evt.gz'))
    #num_zipped_evt_files = len(zipped_evt_files)
    #evt_files = glob.glob((evt_path + '/' + evt_directory + '/' + '*.evt'))
    #num_evt_files = len(evt_files)
    evt_files = glob.glob((evt_path + '/' + evt_directory + '/' + '*.evt*'))
    num_evt_files = len(evt_files)
    if num_evt_files == 0:
      print "WARNING: Found zero evt files in", evt_directory
      ## write LOC_DP_STARTED_FLAG
      open(EVT_DIR + '/' + evt_directory + '/' + LOC_DP_STARTED_FLAG, 'w').close()
      open(EVT_DIR + '/' + evt_directory + '/' + LOC_DP_ERROR_FLAG, 'w').close()
      #return

    # calculate number of cores needed and walltime needed for this dataset

    # How many nodes available on Onsite Cluster.
    total_number_of_cores_on_cluster = 0
    total_number_of_cores_on_cluster += 0  #compute-0-0    **
    total_number_of_cores_on_cluster += 0  #compute-0-1    **
    total_number_of_cores_on_cluster += 4  #compute-0-2
    total_number_of_cores_on_cluster += 0  #compute-0-3    **
    total_number_of_cores_on_cluster += 4  #compute-0-4
    total_number_of_cores_on_cluster += 4  #compute-0-5
    total_number_of_cores_on_cluster += 0  #compute-0-6    **
    total_number_of_cores_on_cluster += 4  #compute-0-7
    total_number_of_cores_on_cluster += 0  #compute-0-8    **
    total_number_of_cores_on_cluster += 0  #compute-0-9    **
    total_number_of_cores_on_cluster += 4  #compute-0-10
    total_number_of_cores_on_cluster += 4  #compute-0-11
    total_number_of_cores_on_cluster += 0  #compute-0-12   **
    total_number_of_cores_on_cluster += 4  #compute-0-13
    total_number_of_cores_on_cluster += 4  #compute-0-14
    total_number_of_cores_on_cluster += 4  #compute-0-15
    total_number_of_cores_on_cluster += 0  #compute-0-16   **
    total_number_of_cores_on_cluster += 4  #compute-0-17
    all_available_nodes = 17

    if num_evt_files > total_number_of_cores_on_cluster:
        num_cores = total_number_of_cores_on_cluster
        
    else:
        num_cores = num_evt_files


    print "I have decided to use", num_cores, "cores to process these", num_evt_files, "evt files."

    proc_time_minutes = num_evt_files * 15 / num_cores
    walltime_str = str(datetime.timedelta(minutes=proc_time_minutes))

    job_name = filename_prefix.split("T")[-1] + "_" + evt_directory + "_cp" + str(cp)
    
    #print "num_cores: ", num_cores, type(num_cores)
    #print "proc_time: ", proc_time_minutes
    #print "wall_time: ", walltime_str


    #
    # The submit_batch_job script is written below this point
    #

    #path_to_auto_run = home_dir + '/LUXcode/Trunk/DataProcessing/ClusterDeployment/Oscar/auto_run_dp_framework_oscar.py'

    #path_to_auto_run = luxcode_path + '/Trunk/DataProcessing/ClusterDeployment/'
    path_to_auto_run = dp_path + '/ClusterDeployment/'
    path_to_auto_run += '/loc/auto_run_dp_framework_loc.py'

    # Make the directory where the logs will go.
    os.mkdir(BASE_LOG_DIR + job_name)

    # Build name of output_file (qsub submission file)
    output_file = BASE_LOG_DIR + job_name + "/job_file.sh"

    # rm old script and open new script
    if os.path.isfile(output_file):
        os.remove(output_file)

    of = open(output_file, 'w')
    """
    """

    #of.write("TEST")
    # write contents to file
    of.write("""#!/bin/bash

# New job submission scheme for LOC queue.
# 
#
# 20130128 - MW  - Modified from submit_batch_job to run on loc

################################################################
# YOU SHOULD NOT EDIT THIS FILE AT ALL!!!                      #
################################################################

#!/bin/bash
#$ -cwd
#$ -j yes
#$ -N T%s
#$ -S /bin/bash
#$ -o /share/data/logs/dplogs/%s
#$ -t 1-%d

# This file handles files:
#     $file_index_stripped - $end_file

source /share/apps/environment.sh

echo "SGE_TASK_ID:    " $SGE_TASK_ID
echo "SGE_TASK_FIRST: " $SGE_TASK_FIRST
echo "SGE_TASK_LAST:  " $SGE_TASK_LAST
echo "SGE_STEP_SIZE:  " $SGE_STEP_SIZE
echo "JOB_ID:         " $JOB_ID
echo "ID:             " $ID

export ID=$((SLURM_TASKS_PER_NODE*(SLURM_ARRAYID-1)+SLURM_PROCID));
echo "Starting job $ID on $HOSTNAME (CPU $SLURM_PROCID)";
echo "Starting job $SGE_TASK_ID on $HOSTNAME (JOB $JOB_ID)";


export EVT_DIRECTORY=%s
export NUM_CORES=%d
export AUTO_RUN_PATH=%s
export CP=%d

echo `which python`
echo "/share/apps/python/pylux/bin/python "

python -u $AUTO_RUN_PATH $EVT_DIRECTORY $CP $SGE_TASK_ID $NUM_CORES $JOB_ID

  """ % (job_name, job_name, num_cores, evt_directory, num_cores, path_to_auto_run, cp))
    of.close()
    return output_file 

  # Args for the above are:
    # 1) The job name which is based on the evt directory name
    # 2) Total number of cores/jobs you'll be making in this job array
    # 3) evt directory name
    # 4) number of cores again.
    # 5) path to the auto_run_dp_framework_loc.py file.
    # 6) The cp number.



def are_we_synced_with_primary_mirror(evt_dir_path, evt_directory):
  # Get number of uncompressed evt files localy.
  local_evt_files = glob.glob(evt_dir_path + "/" + evt_directory + '/' + '*.evt')
  num_local_evt_files = len(local_evt_files)
  local_evtgz_files = glob.glob(evt_dir_path + "/" + evt_directory + '/' + '*.evt.gz')
  num_local_evtgz_files = len(local_evtgz_files)

  # Get number of files on primary mirror
  p = Popen("./get_remote_num_evts.sh " +  evt_directory, shell=True, stdout=PIPE, stderr=PIPE)
  p.wait()
  stdout, stderr = p.communicate()
  #print "Raw:", stdout, stderr
  if stderr:
    print "There was an error in the subprocess:"
    print stderr
  try:
    evt_and_evtgz = stdout.strip().split()
    num_remote_evt = int(evt_and_evtgz[0])
    num_remote_evtgz = int(evt_and_evtgz[1])
    #num_remote_local_evt_files = int(stdout.strip())
  except:
    print "This is where the weird error happened!"
    import traceback
    print traceback.format_exc()
    print "[evt_and_evtgz]:"
    print evt_and_evtgz
    print "[stdout]:"
    print stdout
    print "[stderr]:"
    print stderr
    print "I see", num_local_evt_files, "when this happened."
    print "Setting number of remote files to -2 so this loop continues."
    print ""
    print ""
    #num_remote_local_evt_files = -2
    num_remote_evt = -2
    num_remote_evtgz = -2

  if num_remote_evt < 0 or num_remote_evtgz < 0:
    #print "The primary mirror hasn't finished syncing %s from Lead." % evt_directory
    print "PM still transferring from Lead."
    return False
  #print "Local vs Remote:\t", num_local_evt_files, "vs.", num_remote_local_evt_files
  print "\nLocal vs. Remote:\t%d evt (%d evt.gz) vs. %d evt (%d evt.gz)" % (
      num_local_evt_files, num_local_evtgz_files,
      num_remote_evt, num_remote_evtgz)
  evt_match, evtgz_match = False, False
  evt_match = num_local_evt_files == num_remote_evt
  evtgz_match = num_local_evtgz_files == num_remote_evtgz
  if num_local_evt_files == 0:
    evt_match = False
  if num_local_evtgz_files == 0:
    evtgz_match = False
  if evt_match or evtgz_match:
    return True
  else:
    return False

# run through all directories and if the sync is done
# and the analysis not yet started than start it
SLEEP_SECS_1 = 5
SLEEP_SECS_2 = 5
while True:
    evt_dirs = [f for f in os.listdir(EVT_DIR) \
            if os.path.isdir(os.path.join(EVT_DIR, f))]

    for evt_directory in evt_dirs:
        evt_directory = evt_directory.replace("LUX", "lux")
        filename_prefix = re.match('^(lux10_[0-9]+T[0-9]+)', evt_directory).group()
        # While the loc is running on data synced from the primary mirror, we
        # don't have the option of having a "LOC_SYNC_DONE_FLAG" synced over,
        # so we have to make our own.
        if not os.path.isfile(EVT_DIR + '/' + evt_directory + '/' + LOC_SYNC_DONE_FLAG):
            print "Sync not known to be finished for", evt_directory, "... Checking against primary mirror..."
            sync_status = are_we_synced_with_primary_mirror(EVT_DIR, evt_directory)
            if sync_status == True:
                print "The sync is now done, placing", LOC_SYNC_DONE_FLAG, "flag."
                open(EVT_DIR + '/' + evt_directory + '/' + LOC_SYNC_DONE_FLAG, 'w').close()
        if os.path.isfile(EVT_DIR + '/' + evt_directory + '/' + LOC_SYNC_DONE_FLAG) and not os.path.isfile(EVT_DIR + '/' + evt_directory + '/' + LOC_DP_STARTED_FLAG):

            # print out time and dataset name
            now_time = str(datetime.datetime.now())
            print now_time
            print 'Starting analysis on', evt_directory, 'at', now_time + '...'

            # Load default dp settings file as a string
            data_processing_xml_string = ''.join(open(DEFAULT_DP_SETTINGS_FILE).readlines())
            #print data_processing_xml_string

            # Get eb from file
            if evt_directory == filename_prefix:
                eb = 0;
            else:
                eb = int(re.search('[0-9]+$', evt_directory).group())

            # LUXPrepareDPSettings here
            """
            try:
                cp = LUXPrepareDPSettings(filename_prefix, eb, data_processing_xml_string)
                # create sarray settings file
                create_sarray_settings_file(evt_directory, cp, EVT_DIR, SARRAY_SETTINGS_FILE, PRIORITY) 
            except:
                print 'LUXPrepareDPSettings for ' + evt_directory + ' failed!'
                """
            # NAH.. Let's do this.
            cp = 0
            print "="*80
            print evt_directory
            print "="*80
            try:
              cp = LUXPrepareDPSettings(filename_prefix, eb, data_processing_xml_string)
              print "cp: ", type(cp), cp
              qsub_file = create_sarray_settings_file(evt_directory, cp, EVT_DIR, PRIORITY) 
              #write cp number to event file
              cpFile = open('/share/data/logs/cplogs/'+evt_directory,'a+')
              cpFile.write(str(cp) + '\n')
              cpFile.close()
              # write LOC_DP_STARTED_FLAG
              open(EVT_DIR + '/' + evt_directory + '/' + LOC_DP_STARTED_FLAG, 'w').close()
              # run qsub on job file
              p = Popen(['qsub', qsub_file])
              p.wait()
              print "Job submitted:", evt_directory, "-->", qsub_file
            except Exception as inst:
                print "!"*75
                print 'LUXPrepareDPSettings for ' + evt_directory + ' failed!'
                print 'Cause of error:'
                import traceback
                print traceback.format_exc()
                open(EVT_DIR + '/' + evt_directory + '/' + LOC_DP_STARTED_FLAG, 'w').close()
                open(EVT_DIR + '/' + evt_directory + '/' + LOC_DP_ERROR_FLAG, 'w').close()
                print "Flags printed. This dataset will be skipped."
                print "!"*75

            # This section was to slow down stuff and study the cause of flow breakdown.
            """
            resp = raw_input("Continue? ")
            if resp == "n" or resp == "":
              import sys
              sys.exit()
            if resp == "f":
              open(EVT_DIR + '/' + evt_directory + '/' + LOC_DP_STARTED_FLAG, 'w').close()
            """

            # sleep before submitting next job
            time.sleep(SLEEP_SECS_2)
            print '\n'

    # sleep before looping through directories again
    time.sleep(SLEEP_SECS_1)
