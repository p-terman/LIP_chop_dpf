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
20190402 PAT has overhauled this 
'''


# import required modules
#from subprocess import Popen
import os
import re
import sys
#import LUXGenerateRQsPyMod PAT removed
#from LUXModuleWrapperPyMod import GetGlobalSettings
#from ReportLUXcodePathPyMod import ReportLUXcodePath
#import LUXLUGQueriesPyMod
from glob import glob
#import matlab.engine
import gzip
import shutil
#PAT added matlab import and gzip shutil

# call script to set python paths
home_dir = '/gpfs_home/pterman/'
matlab_path = '/gpfs/runtime/opt/matlab/R2017b/bin/matlab'
execfile(home_dir + '/LUXCode/Trunk/DataProcessing/DataProcessingFramework/Utilities/LUXSetPythonPaths.py')

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
if len(argv) is not 8:
  help()
if "-h" in argv or "--help" in argv:
  help()
filename_prefix = argv[1]
data_path_evt = argv[2]
data_path_rq = argv[3] # this will have to be the full path to the rq files
LUXcode_path = argv[4]
data_processing_xml_path = argv[5]
cp = argv[6]
ClusterID = argv[7]

#data_path_evt = evt_root_path + '/' + dataset_name + '/'
#rq_file_path = rq_root_path + '/' + filename_prefix + '_cp%0.5d/' % cp PAT removed-


# find all .evt files in data_path_evt and throw them into a list 
#evt_file_list = list()
evt_files = glob(data_path_evt + '/*evt')
#evtgz_files = glob(data_path_evt + '/*evt.gz')
rq_files = glob(data_path_rq + '/*.rq.mat')

#for f in evt_files:
#    evt_file_list.append(os.path.basename(f))
#for f in evtgz_files: #PAT added / adapted to unzip files if there are any
#	g = gzip.open(f, 'rb')
#	file_content = g.read()
#	file_new_name = f.replace('.gz', '')
#	file_new = open(file_new_name, "w")
#	file_new.write(file_content)
#	file_new.close()
#	g.close()

iq_string = '/gpfs_home/pterman//LUXCode/Stable_Releases/DataProcessingStableReleases/v2.0/DataProcessing/DataProcessingFramework/LUGqueryIQs.py %s ' % cp
#I know this is crude but I am tired
os.system('python ' + iq_string)

print iq_string


#for q in evt_files: #PAT added to copy unzipped files to the working directory for evt rq and new folders
#	shutil.copy2(q, os.getcwd() + '/' + filename_prefix)


#then all of the files should be in the same folder in home dir 
#evt files in the folder with the pathbase name and within that a folder called matfiles

	
#    if f.rstrip('.gz') not in evt_file_list:
#        evt_file_list.append(os.path.basename(f))
#evt_file_list.sort()

iq_xml_path ='/users/pterman/lug_iqs_cp' + cp + '.xml' 
shutil.copy2('/users/pterman/lug_iqs.xml', iq_xml_path)

#eng = matlab.engine.start_matlab()
#eng.select_reprocessing_dpf_brown(LUXcode_path, data_path_evt, data_path_rq, data_processing_xml_path, iq_xml_path)
comm = list()
quote = "'"
double_quote = '"'
matlab_line_1 = "{}select_reprocessing_dpf_brown_Cluster({}{}{}, {}{}{}, {}{}{}, {}{}{}, {}{}{}, {}{}{}) {}".format(double_quote, quote,LUXcode_path,quote, quote, data_path_evt, quote, quote, data_path_rq, quote, quote, data_processing_xml_path, quote, quote, iq_xml_path, quote, quote, ClusterID ,quote, double_quote) #this line is the main program
print matlab_line_1
comm.append(matlab_line_1)
#comm.append("select_reprocessing_dpf_brown('%s','%s','%s','%s','%s'); ".format(LUXcode_path, data_path_evt, data_path_rq, data_processing_xml_path, iq_xml_path)) #this line is the main program
#comm.append(" exit")
#matlab_comm = ''.join(comm)
#Popen([matlab_path,  '-nojvm', '-nodesktop', '-nosplash', '-r', matlab_comm]) 
#I don't get why Popen isn't working. Tired of this, so I'm brute forcing
os.system('matlab -nojvm -nodesktop -nosplash -r ' + matlab_line_1)

#add something to remove unzipped and copied files created. PAT

#for ff in evtgz_files:
#	file_to_remove = ff.replace('.gz', '')
#	os.remove(file_to_remove)

#for qq in evt_files:
#	os.remove(qq)

#for rr in rq_files:
#	os.remove(rr)

dst = '/users/pterman/data/pterman/SRP_out/' + filename_prefix +'_cp' + cp 

try:
	os.mkdir(dst)
	os.mkdir(dst + '/chopped_events/')	
	os.mkdir(dst + '/chopped_events/matfiles')
except:
	print('Folder {} already exists'.format(dst))
src = '/users/pterman/scratch/' + filename_prefix +'_cp' + cp + '/chopped_events/'
chopped_rq = glob(src + '/*.rq')
chopped_mat = glob(src + '/matfiles/*.rq.mat')
print(dst + '/chopped_events/')

for q in chopped_rq:
	print(q)
	try:
	        shutil.copy2(q, dst + '/chopped_events/')
	except:
		print('RQ file {} cannot be copied'.format(q))
for r in chopped_mat:
	print(r)
	try:
	       	shutil.copy2(r, dst + '/chopped_events/matfiles/')
	except:
		print('MAT file {} cannot be copied'.format(r))

#end result should be an otherwise empty folder with a 'chopped_files' folder with the chop.rq files 
#and within that there is a 'matfiles' folder with the chop.rq.mat files



