"""
CleanUp_DefaultPyMod.py

This is a simple Python module that doesn't require any binary file I/O. It is not a
module template.

This module does basic clean up operations after the DP framework is done.
Currently:
    1) Deletes cvt file
    2) Gzip rq file
    3) Add line to progress status file for each created rq

2013-02-19 - JRV - Created
2013-11-22 - AL  - Added section to write rq file name to progress status file
"""

import os
import gzip
import sys
import re


def gz(file_path_in, del_old=1):
    file_path_out = None
    if os.path.isfile(file_path_in):
        f_in = open(file_path_in, 'rb')
        file_path_out = file_path_in + '.gz'
        f_out = gzip.open(file_path_out, 'wb')
        f_out.writelines(f_in)
        f_out.close()
        f_in.close()
        if del_old:
            os.remove(file_path_in)
    return file_path_out


def CleanUp_Default(evt_file_name, evt_file_path, rq_file_name, rq_file_path, run_module_number, dp_settings_xml_path, lug_iqs_xml_path):

    fexists = os.path.isfile(rq_file_path + '/' + rq_file_name)

    # 1) delete cvt file
    print "Deleting cvt file"
    cvt_file_name = evt_file_name.replace('evt', 'cvt')

    if os.path.isfile(evt_file_path + '/' + cvt_file_name):
        try:
            os.unlink(evt_file_path + '/' + cvt_file_name)
        except:
            print "ERROR: Couldn't delete file %s" % cvt_file_name

    # 2) gzip rq file
    print "gzipping rq"
    gz(rq_file_path + '/' + rq_file_name)

    # 3) write to the "completed jobs" list
    # First of all, we need to check if the rq file exists!
    if fexists:
        print "Write complete flag to status file"
        # Let's check if the directory exists
        #status_dir = evt_file_path + '../../JobsInProgressList/'
        # for running on Oscar, use the following line
        status_dir = '/users/jverbus/scratch/JobsInProgressList/'
        print "Status file directory: " + status_dir
        if (os.path.isdir(status_dir)):
            # let's get the dataset name from evt_file_name
            dset, dummy = re.split('_f', evt_file_name)
            # and now the cp number from rq_file_path
            dummy, cp1 = re.split('_cp', rq_file_name)
            cp, dummy = re.split('.rq', cp1)
            # now we know how to name the progress status file
            status_file_name = status_dir + dset + '_cp' + str(cp)
            status_file = open(status_file_name,'a')
            new_rq_file_name = rq_file_name[1:]
            status_file.write(new_rq_file_name + '\n')
            status_file.close()
        else:
            print "CleanUp module WARNING: Couldn't write to status file"
    else:
        print "CleanUp module WARNING: couldn't find rq file!"

