#!/usr/bin/python
'''
This script is a very simple example outlining how to run
an early version of the dp framwork (v0.9). All lines that
need to be changed for your environment are labeled.

2012-11-28 - JRV - Created
2012-12-20 - JRV - Added LUG queries (v0.9)
'''

# call script to set python paths
execfile('/private/tmp/eli/DataProcessing/DataProcessingFramework/Utilities/LUXSetPythonPaths.py')

# import required modules
import os
import re
import LUXGenerateRQsPyMod
from subprocess import Popen
from LUXModuleWrapperPyMod import GetGlobalSettings
from ReportLUXcodePathPyMod import ReportLUXcodePath
import LUXLUGQueriesPyMod
from glob import glob

# Create syntax function
def help():
  print "This is a simple program for running a single dataset at a time"
  print "through the data processing framework."
  print "SYNTAX:"
  print "\t%s [-h/--help] dataset path_to_evt_files path_to_put_rqs" % os.sys.argv[0]
  print "EXAMPLE:"
  print "\t%s lux10_20130213T0852_eb00009 ~/myevts/ ~/myrqs/" % os.sys.argv[0]
  os.sys.exit()

# Parse input values
argv = os.sys.argv
if len(argv) is not 4:
  help()
if "-h" in argv or "--help" in argv:
  help()
dataset_name = argv[1]
evt_root_path = argv[2]
rq_root_path = argv[3]

# you need to modify these paths for your enviornment
filename_prefix = re.match('^(lux.._[0-9]+T[0-9]+)', dataset_name).group()
#filename_prefix = 'lux10_20121214T1327'  # YOU NEED TO CHANGE THIS
dp_settings_xml_path = '/tmp/eli/DataProcessing/PACMAN_settings.xml'

# do all required LUXPrepare LUG queries
eb = 0  # event builder record number is 0
data_processing_xml_string = ''.join(open(dp_settings_xml_path).readlines())
cp = LUXLUGQueriesPyMod.LUXPrepareDPSettings(filename_prefix, eb, data_processing_xml_string)

evt_file_path = evt_root_path + '/' + dataset_name + '/'
rq_file_path = rq_root_path + '/' + filename_prefix + '_cp%0.5d/' % cp

# find all .evt files in evt_file_path and throw them into a list
evt_file_list = list()
evt_files = glob(evt_file_path + '/*evt')
evtgz_files = glob(evt_file_path + '/*evt.gz')

for f in evt_files:
    evt_file_list.append(os.path.basename(f))
for f in evtgz_files:
    if f.rstrip('.gz') not in evt_file_list:
        evt_file_list.append(os.path.basename(f))
evt_file_list.sort()

# run the data processing framework
output_logs = LUXGenerateRQsPyMod.LUXGenerateRQs(evt_file_path, rq_file_path, cp, evt_file_list)
