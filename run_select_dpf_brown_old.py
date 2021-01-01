#!/usr/bin/python
'''
This script is a very simple example outlining how to run
an early version of the dp framwork (v0.9). All lines that
need to be changed for your environment are labeled.

2012-11-28 - JRV - Created
2012-12-20 - JRV - Added LUG queries (v0.9)
2015-02-17 - AL  - Modified to maintain compatibility with LUXPrepareDPSettings (last argument changed from being a string
                    containing the full xml settings to being just the path to the xml settings file)
20190217 Edited for LIP chopping PAT
'''

# call script to set python paths
execfile('~/DataProcessing/DataProcessingFramework/Utilities/LUXSetPythonPaths.py')

# import required modules
import os
import re
#import LUXGenerateRQsPyMod PAT removed
from subprocess import Popen
from LUXModuleWrapperPyMod import GetGlobalSettings
from ReportLUXcodePathPyMod import ReportLUXcodePath
import LUXLUGQueriesPyMod
from glob import glob
import matlab.engine
import gzip
import shutil
#PAT added matlab import and gzip shutil

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
filename_prefix = argv[1]
evt_file_path= argv[2]
data_path_rq = argv[3] # this will have to be the full path to the rq files
LUXcode_path = argv[4]
data_processing_xml_path = argv[5]

# you need to modify these paths for your enviornment
#filename_prefix = re.match('^(lux.._[0-9]+T[0-9]+)', dataset_name).group() PAT removed
#filename_prefix = 'lux10_20121214T1327'  # YOU NEED TO CHANGE THIS
# '/tmp/eli/DataProcessing/PACMAN_settings.xml'

# do all required LUXPrepare LUG queries
eb = 0  # event builder record number is 0
cp = LUXLUGQueriesPyMod.LUXPrepareDPSettings(filename_prefix, eb, data_processing_xml_path)

#data_path_evt = evt_root_path + '/' + dataset_name + '/'
#rq_file_path = rq_root_path + '/' + filename_prefix + '_cp%0.5d/' % cp PAT removed-


# find all .evt files in evt_file_path and throw them into a list 
#evt_file_list = list()
#evt_files = glob(evt_file_path + '/*evt')
evtgz_files = glob(evt_file_path + '/*evt.gz')

#for f in evt_files:
#    evt_file_list.append(os.path.basename(f))
for f in evtgz_files: #PAT added / adapted to unzip files if there are any
	g = gzip.open(f, 'rb')
	file_content = g.read()
	file_new_name = f.replace('.gz', '')
	file_new = open(file_new_name, "w")
	file_new.write(file_content)
	file_new.close()
	g.close()

os.mkdir(os.getcwd() + '/' + filename_prefix)

evt_files = glob(evt_file_path + '/*evt')

for q in evt_files: #PAT added to copy unzipped files to the working directory for evt rq and new folders
	shutil.copy2(q, os.getcwd() + '/' + filename_prefix)

rq_files = glob(data_path_rq + '/matfiles/*.rq.mat' )
os.mkdir(os.getcwd() + '/' + filename_prefix + '/matfiles/')

for r in rq_files:
	shutil.copy2(q, os.getcwd() + '/' + filename_prefix + '/matfiles/')

#then all of the files should be in the same folder in home dir 
#evt files in the folder with the pathbase name and within that a folder called matfiles

	
#    if f.rstrip('.gz') not in evt_file_list:
#        evt_file_list.append(os.path.basename(f))
#evt_file_list.sort()

# run the data processing framework
#output_logs = LUXGenerateRQsPyMod.LUXGenerateRQs(evt_file_path, rq_file_path, cp, evt_file_list)

#does anything need to be moved/copied? Will need to change file paths then

#PAT added these two lines
eng = matlab.engine.start_matlab()
eng.select_reprocessing_dpf_brown(LUXcode_path, data_path_evt, data_path_rq, data_processing_xml_path, iq_xml_path)

#add something to remove unzipped and copied files created. PAT

for ff in evtgz_files:
	file_to_remove = ff.replace('.gz', '')
	os.remove(file_to_remove)

remove_evt_list = glob(os.getcwd() + '/' + filename_prefix + '/*evt')
for qq in remove_evt_list:
	os.remove(qq)

remove_rq_list = glob(q, os.getcwd() + '/' + filename_prefix + '/matfiles/*.rq.mat')
for rr in remove_rq_list:
	os.remove(rr)

#end result should be an otherwise empty folder with a 'chopped_files' folder with the chop.rq files 
#and within that there is a 'matfiles' folder with the chop.rq.mat files







