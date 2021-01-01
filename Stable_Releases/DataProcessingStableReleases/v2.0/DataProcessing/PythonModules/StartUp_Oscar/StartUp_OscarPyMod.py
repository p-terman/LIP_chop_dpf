"""
StartUp_OscarPyMod.py

This is a simple Python module that doesn't require any binary file I/O. It is not a
module template.

This module does basic per file preperation
Currently:
    1) Checks if evt file is available in evt directory
    2) If .evt file does not exist, then it syncs it from gsk-41
    3) Gunzip .evt.gz file if nesassary

2013-04-23 - JRV - Created
"""

import os
import gzip


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


def get_evt_file(remote_host, remote_evt_paths, evt_file_name, evt_file_path):
    # loop through remote_evt_paths
    for p in remote_evt_paths:
        # determine evt directory and check if it exists
        if os.path.isdir(evt_file_path):
            # rsync .evt or .evt.gz file over
            pass


def StartUp_Oscar(evt_file_name, evt_file_path, rq_file_name, rq_file_path, run_module_number, dp_settings_xml_path, lug_iqs_xml_path):

    status = []

    # 1) check for evt file
    if os.path.isfile(evt_file_path + '/' + evt_file_name):
        return status
    elif os.path.isfile(evt_file_path + '/' + evt_file_name + '.gz'):
        # gunzip file

        return status
    else:
        # sync file from
        pass

    # 2) gunzip evt file if nessasary
