#!/usr/bin/python
"""
DataProcessingDaemonOscar.py

This script runs in a screen on an Oscar login node.
In v1.0 of the DPF it looks for new datasets, creates
the submit_batch_job file, and calls qsub to submit the
job. It is the equivalent of ccv_daemon.sh in the Run02
framework.

In v2.0 of the DPF it will be expanded to do a number
of LUG queries for processing requests/settings.

2012-12-05 - JRV - Created
2012-12-10 - CHF - Added functionality for preparing the DP settings
2013-01-17 - JRV - Added function to produce submit_batch_job_SLURM files
2013-04-11 - JRV - Now uses humanize_time() to fix DP bug #35
"""

# Set Python Paths
import os
home_dir = os.path.expanduser("~")
execfile(home_dir + '/LUXcode/Trunk/DataProcessing/DataProcessingFramework/Utilities/LUXSetPythonPaths.py')

# import modules
import datetime
import time
from subprocess import Popen
from LUXLUGQueriesPyMod import LUXPrepareDPSettings
import re
from ReportLUXcodePathPyMod import ReportDataProcessingPath
from xml2dict import xml2dict
from LUXDatabaseSQLPyMod import LUXDatabaseSQL
from str2numPyMod import str2num

# define things

EVT_DIR = home_dir + '/scratch/evt/'
CCV_DP_STARTED_FLAG = 'ccv_dp_started'
CCV_SYNC_DONE_FLAG = 'ccv_sync_done'
MANUAL_SUBMIT_FLAG = 'ccv_dp_manual'
SARRAY_SETTINGS_FILE = home_dir + '/submit_batch_job_SLURM'
DEFAULT_DP_SETTINGS_FILE = home_dir + '/LUXcode/Trunk/DataProcessing/default_data_processing_settings.xml'
#PRIORITY = 'pri-cah3'
PRIORITY = 'pri-jverbus'

def humanize_time(secs):
    """
    Converts seconds to a string formatted as 'HH:MM:SS'
    """
    mins, secs = divmod(secs, 60)
    hours, mins = divmod(mins, 60)
    return '%02d:%02d:%02d' % (hours, mins, secs)


def create_sarray_settings_file(evt_directory, cp, evt_path, output_file, priority):
    """
    This function creates the sarray settings file. It is the equivalent
    to create_submit_batch_job.sh in the Run02 analysis framework.
    """

    # import required modules

    import sys
    import os
    import glob
    import datetime
    home_dir = os.path.expanduser("~")

    # make sure evt_path is a directory and has ccv_sync_done flag in it
    if not os.path.isdir(evt_path + '/' + evt_directory):
        sys.exit("ERROR: Filename prefix doesn't correspond to a directory in the evt path")

    if not evt_directory.startswith('lux'):
        sys.exit("ERROR: Input must be a valid filename prefix")

    if not os.path.isfile( (evt_path + '/' + evt_directory + '/' + 'ccv_sync_done') ):
    	sys.exit("ERROR: CCV sync has not completed! Exiting...")

    filename_prefix = re.match('^(lux10_[0-9]+T[0-9]+)', evt_directory).group()

    # get list of evt files
    evt_files = glob.glob((evt_path + '/' + evt_directory + '/' + '*.evt'))
    num_evt_files = len(evt_files)

    # calculate number of cores needed and walltime needed for this dataset
    if num_evt_files > 800:
        num_nodes = 32
        cores_per_node = 8
        node_string = ('1-' + str(num_nodes))
    elif num_evt_files > 400:
        num_nodes = 16
        cores_per_node = 8
        node_string = ('1-' + str(num_nodes))
    elif num_evt_files > 200:
        num_nodes = 8
        cores_per_node = 8
        node_string = ('1-' + str(num_nodes))
    elif num_evt_files > 32:
        num_nodes = 4
        cores_per_node = 8
        node_string = ('1-' + str(num_nodes))
    elif num_evt_files > 16:
        num_nodes = 2
        cores_per_node = 8
        node_string = ('1-' + str(num_nodes))
    elif num_evt_files > 8:
        num_nodes = 1
        cores_per_node = 8
        node_string = '1'
    elif num_evt_files > 2:
        num_nodes = 1
        cores_per_node = 2
        node_string = '1'
    else:
        num_nodes = 1
        cores_per_node = 1
        node_string = '1'

    num_cores = num_nodes * cores_per_node
    proc_time_seconds = 60 * num_evt_files * 10 / num_cores  # magic!
    walltime_str = humanize_time(proc_time_seconds)

    #
    # The submit_batch_job script is written below this point
    #

    path_to_auto_run = home_dir + '/LUXcode/Trunk/DataProcessing/ClusterDeployment/Oscar/auto_run_dp_framework_oscar.py'

    # rm old script and open new script
    if os.path.isfile(output_file):
        os.remove(output_file)

    of = open(output_file, 'w')

    # write contents to file
    of.write("""#!/bin/bash
#
# New job submission scheme for CCV queue.
# Only modify this file, there is no need to modify the .m file.
#
# 20130114 - JRV - Modified from submit_batch_job to run on new Oscar cluster using SLURM

################################################################
# YOU SHOULD NOT EDIT THIS FILE AT ALL!!!                      #
# There is no more need to modify or duplicate run_analysis.m  #
################################################################

# name your job to help you and others identify it in the queue
#SBATCH --job-name=%s

# set mail address
#SBATCH --mail-user=james_verbus@brown.edu

# estimate the time needed by your job, as HH:MM:SS
#SBATCH --time=%s

# set which partition
##SBATCH -p sandy-batch

## set priority
#SBATCH --qos=%s
##SBATCH --qos=rgaitske-sb-condo

# set cores per job
#SBATCH --ntasks=%d
#SBATCH --ntasks-per-node=%d

# set memory per cpu
#SBATCH --mem-per-cpu=4G

# constrain to intel nodes only
#SBATCH --constraint="intel"

# set number of jobs in array
#SARRAY --range=%s

# [optional] combine stdout and stderr streams into one output file

# [optional] specify a different output file than the default JobName.out

# set path to auto run file
PATH_TO_AUTO_RUN='%s'

# set the filename prefix
EVT_DIRECTORY='%s'

# set the number of cores (NODES * ppn): e.g. NODES = 0-3 & ppn=4 =>
NUM_CORES=%d

# cp id
CP=%d

#######################################################
# YOU DO NOT NEED TO CHANGE ANYTHING BELOW THIS POINT #
#######################################################

srun bash -l -c ''' 
    cd $PWD
    export ID=$((SLURM_TASKS_PER_NODE*(SLURM_ARRAYID-1)+SLURM_PROCID));
    echo "Starting job $ID on $HOSTNAME (CPU $SLURM_PROCID)";
                          
    export PATH_TO_AUTO_RUN=`echo -e "$0"`                      
    export EVT_DIRECTORY=`echo -e "$1"`
    export CP=`echo -e "$2"`
    export NUM_CORES=`echo -e "$3"`

    python $PATH_TO_AUTO_RUN $EVT_DIRECTORY $CP $ID $NUM_CORES
''' "$PATH_TO_AUTO_RUN" "$EVT_DIRECTORY" "$CP" "$NUM_CORES"
""" % (filename_prefix, walltime_str, priority, cores_per_node, cores_per_node, node_string, path_to_auto_run, evt_directory, num_cores, cp))

    # order is: filename_prefix, walltime_str, priority, cores_per_node, cores_per_node, node_string, filename_prefix, num_cores, cp

    of.close()


def check_run_dp(filename_prefix):

    xml_cred = xml2dict(ReportDataProcessingPath() + 'DataProcessingFramework/Utilities/lug_query_credentials_readonly.xml')
    db = LUXDatabaseSQL(xml_cred['credentials']['host'], xml_cred['credentials']['user'], xml_cred['credentials']['password'], xml_cred['credentials']['database'], str2num(xml_cred['credentials']['port']))
    db.connect()

    try:
        query1 = "SELECT entry_id FROM lug_acquisitions WHERE filename_prefix = '%s' LIMIT 1" % filename_prefix
        result1 =  db.query_db(query1)
        lug_entry_id = int(result1[0][0])

        query2 = "SELECT run_dp_framework FROM daq_control WHERE lug_entry_id = %d LIMIT 1" % lug_entry_id
        result2 = db.query_db(query2)
        run_dp_flag = int(result2[0][0])
    except:
        print "ERROR: check_run_dp LUG queries failed (was the DAQ GUI used?). Not submitting..."
        run_dp_flag = 0

    return run_dp_flag


# run through all directories and if the sync is done
# and the analysis not yet started than start it
SLEEP_SECS_1 = 5
SLEEP_SECS_2 = 5
while True:
    evt_dirs = [f for f in os.listdir(EVT_DIR) \
            if os.path.isdir(os.path.join(EVT_DIR, f))]

    for evt_directory in evt_dirs:
        filename_prefix = re.match('^(lux10_[0-9]+T[0-9]+)', evt_directory).group()
        if os.path.isfile(EVT_DIR + '/' + evt_directory + '/' + CCV_SYNC_DONE_FLAG) and not os.path.isfile(EVT_DIR + '/' + evt_directory + '/' + CCV_DP_STARTED_FLAG):

            # check to see if we should process this dataset
            run_dp_flag = check_run_dp(filename_prefix)
            if run_dp_flag == 0 and not os.path.isfile(EVT_DIR + '/' + evt_directory + '/' + MANUAL_SUBMIT_FLAG):
                print 'Skipping ' + evt_directory + ' as requested in DAQ GUI'
                # write CCV_DP_STARTED_FLAG -- we probably should have a special flag for this case
                open(EVT_DIR + '/' + evt_directory + '/' + CCV_DP_STARTED_FLAG, 'w').close()
                continue

            # print out time and dataset name
            now_time = str(datetime.datetime.now())
            print now_time
            print 'Starting analysis on ' + filename_prefix + '...'

            # Load default dp settings file as a string
            data_processing_xml_string = ''.join(open(DEFAULT_DP_SETTINGS_FILE).readlines())

            # Get eb from file
            if evt_directory == filename_prefix:
                eb = 0;
            else:
                eb = int(re.search('[0-9]+$', evt_directory).group())

            # LUXPrepareDPSettings here
            try:
                cp = LUXPrepareDPSettings(filename_prefix, eb, data_processing_xml_string)
                create_sarray_settings_file(evt_directory, cp, EVT_DIR, SARRAY_SETTINGS_FILE, PRIORITY) 
            except:
                print 'LUXPrepareDPSettings for ' + filename_prefix + ' failed!'
                if os.path.isfile(SARRAY_SETTINGS_FILE):
                    os.unlink(SARRAY_SETTINGS_FILE)

            # write CCV_DP_STARTED_FLAG
            open(EVT_DIR + '/' + evt_directory + '/' + CCV_DP_STARTED_FLAG, 'w').close()

            # run sarray
            p = Popen(['sarray', SARRAY_SETTINGS_FILE])
            p.wait()

            # sleep before submitting next job
            time.sleep(SLEEP_SECS_2)
            print '\n'

    # sleep before looping through directories again
    time.sleep(SLEEP_SECS_1)
