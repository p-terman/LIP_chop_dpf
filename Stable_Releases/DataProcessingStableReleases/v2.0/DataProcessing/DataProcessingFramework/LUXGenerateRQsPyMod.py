"""
LUXGenerateRQs.py

This file contains the LUXGenerateRQs function.

2012-11-08 - JRV - Created
"""

# import needed python modules
import glob
import re
import LUXModuleWrapperPyMod
import LUXLUGQueriesPyMod
from LUXMachineLookupPyMod import LUXMachineLookup
import string
import random
import os
from errno import EEXIST
import gzip
import time
from shutil import move, copy, rmtree

def ugz(file_path_in, del_old=1):
    file_path_out = None
    if os.path.isfile(file_path_in) and file_path_in.endswith('.gz'):
        f_in = gzip.open(file_path_in, 'rb')
        file_path_out = file_path_in.rstrip('.gz')
        f_out = open(file_path_out, 'wb')
        f_out.writelines(f_in)
        f_out.close()
        f_in.close()
        if del_old:
            try:
                os.unlink(file_path_in)
            except:
                print "ERROR: function 'ugz()' - Cannot delete " + file_path_in
    else:
        file_path_out = file_path_in
    
    return file_path_out


class SmartRQPath:
    """ This class will transparently handle ramdisk use on
    Oscar. 
    
    20130425 - JRV - Created
    """

    def __init__(self, rq_dir, evt_dir):
        self.rq_dir = rq_dir
        self.evt_dir = evt_dir
        self.use_ramdisk = False
        self.rq_dir_to_use = self.rq_dir
        self.evt_dir_to_use = self.evt_dir
        self.create_directory(self.get_rq_path())
    
    def setup_ramdisk(self, ramdisk_path):

        self.ramdisk_path = ramdisk_path
        
        if os.access(self.ramdisk_path, os.W_OK):  # if RAM disk is writable 
            self.use_ramdisk = True
            self.rq_dir_to_use = self.ramdisk_path + os.path.basename(os.path.normpath(self.rq_dir)) + '/'
            self.create_directory(self.get_rq_path())

            self.evt_dir_to_use = self.ramdisk_path + os.path.basename(os.path.normpath(self.evt_dir)) + '/'
            self.create_directory(self.get_evt_path())

        else:
            self.use_ramdisk = False
            self.rq_dir_to_use = self.rq_dir
            self.evt_dir_to_use = self.evt_dir
            self.create_directory(self.get_rq_path())

    def using_ramdisk(self):
        return self.use_ramdisk
   
    def delete_old_ramdisk_dirs(self, age_seconds):
        limit_time_seconds = time.time() - age_seconds
        for d in os.listdir(self.ramdisk_path):
                try:
                   print d
                   print os.path.getmtime(self.ramdisk_path + '/' + d) < limit_time_seconds
                   if d.startswith('lux') and os.path.getmtime(self.ramdisk_path + '/' + d) < limit_time_seconds:
                        print "Deleteing ramdisk dir: " + d
                        rmtree(self.ramdisk_path + '/' + d)
                except:
                    print "WARNING couldn't delete old ramdisk dir " + d

    def create_directory(self, d):
        # make rq directory if it doesn't exist
        try:
            os.makedirs(d)
            os.chmod(d, 0o775)
        except OSError, e:
            if e.errno != EEXIST:
                raise
    
    def prepare_evt_path(self, evt_file):
        if self.use_ramdisk:
            try:
                copy(self.evt_dir + '/' + evt_file, self.evt_dir_to_use + '/' + evt_file) 
                os.chmod(self.evt_dir_to_use + '/' + evt_file, 0770)
            except:
                print "ERROR: Couldn't copy .evt to ramdisk"

        if self.use_ramdisk:
            unzipped_evt_path = ugz(self.evt_dir_to_use + '/' + evt_file, 1)
        else:
            unzipped_evt_path = ugz(self.evt_dir_to_use + '/' + evt_file, 0)

        return os.path.basename(unzipped_evt_path.rstrip('/'))

    def get_evt_path(self):
        return self.evt_dir_to_use

    def get_final_evt_path(self):
        return self.evt_dir
    
    def get_rq_path(self):
        return self.rq_dir_to_use

    def get_final_rq_path(self):
        return self.rq_dir

    def clean_up(self, rq_file, evt_file):
        if self.use_ramdisk:
            # cleanup .rq files
            temp_rq_files = glob.glob(self.rq_dir_to_use + '/' + os.path.splitext(rq_file)[0] + '*.rq.gz')
            temp_mat_files = glob.glob(self.rq_dir_to_use + '/matfiles/' + os.path.splitext(rq_file)[0] + '*.mat')
            temp_root_files = glob.glob(self.rq_dir_to_use + '/rootfiles/' + os.path.splitext(rq_file)[0] + '*.root')
            temp_hdf5_files = glob.glob(self.rq_dir_to_use + '/hdf5/' + os.path.splitext(rq_file)[0] + '*.hdf5')

            self.create_directory(self.rq_dir)
            self.create_directory(self.rq_dir + '/matfiles/')
            self.create_directory(self.rq_dir + '/rootfiles/')
            self.create_directory(self.rq_dir + '/hdf5/')

            for f in temp_rq_files:
                move(f, self.rq_dir + '/' + os.path.basename(f))
            for f in temp_mat_files:
                move(f, self.rq_dir + '/matfiles/' + os.path.basename(f))
            for f in temp_root_files:
                move(f, self.rq_dir + '/rootfiles/' + os.path.basename(f))
            for f in temp_hdf5_files:
                move(f, self.rq_dir + '/hdf5/' + os.path.basename(f))

            # cleanup evt files
            try:
                os.unlink(self.evt_dir_to_use + '/' + evt_file)
            except:
                print "ERROR: Couldn't delete temp evt file"


def LUXGenerateRQs(evt_dir, rq_dir, cp, evt_file_list):
    """
    LUXGenerateRQs

    Inputs:
                        evt_dir - path to the evt directory (string)
                         rq_dir - path to the rq directory (string)
                             cp - complete process record id (integer)
                  evt_file_list - evt_files for which this job is responsible (list of file names)
    
    2013-03-08 - JRV - Gzipping removed; now done in CleanUp module
    2013-04-24 - JRV - Replaced rq_dir instances with SmartRQPath class
    """

    # this prevents up to ~200 simultaenous SQL queries
    time.sleep(random.randrange(30)) 
    
    # create default rq directory
    smart_rq_path = SmartRQPath(rq_dir, evt_dir)

    # write xml files for this specific job
    data_processing_xml_string, iq_xml_string = LUXLUGQueriesPyMod.LUXGetDPSettings(cp)
    data_processing_xml_path = smart_rq_path.get_final_rq_path() + '.data_processing_settings_' + rand_str_generator() + '.xml'
    iq_xml_path = smart_rq_path.get_final_rq_path() + '.lug_iqs_' + rand_str_generator() + '.xml'

    f = open(data_processing_xml_path, 'w')
    f.write(data_processing_xml_string)
    f.close()
    os.chmod(data_processing_xml_path, 0o775)

    f = open(iq_xml_path, 'w')
    f.write(iq_xml_string)
    f.close()
    os.chmod(iq_xml_path, 0o775)

    # load data_processing_settings.xml file
    module_settings, mod_num_map = LUXModuleWrapperPyMod.GetModuleSettings(data_processing_xml_path)
    global_settings = LUXModuleWrapperPyMod.GetGlobalSettings(data_processing_xml_path)

    # setup ramdisk if possible
    if global_settings.has_key('ramdisk_path'):
        smart_rq_path.setup_ramdisk(global_settings['ramdisk_path'])
        smart_rq_path.delete_old_ramdisk_dirs(3600 * 6)  # delete if older than 6 hours (argument in seconds)
        
    # loop over all evt files for which this job is responsible
    output_logs = dict()
    for evt_file_in in evt_file_list:

        start_file_time = time.time()
        
        # this avoids crashes if two processes touch at same time
        try:
            os.utime(evt_dir, (start_file_time, start_file_time))
        except:
            pass
        try:
            os.utime(smart_rq_path.get_final_rq_path(), (start_file_time, start_file_time))
        except:
            pass

        # copy evt file to ramdisk if possible
        evt_file = smart_rq_path.prepare_evt_path(evt_file_in)

        # create rq_file string from evt_file
        if len(evt_file) == 34:
            rq_file = re.sub('\.evt$', '_cp%.5d.rq' % cp, evt_file)
        else:
            rq_file = re.sub('_eb[0-9]+\.evt$', '_cp%.5d.rq' % cp, evt_file)

        output_logs[evt_file] = dict()
        # loop over all modules
        print '\n'
        for ii in range(1, len(module_settings) + 1):
            start_module_time = time.time()
            output_logs[evt_file][module_settings[mod_num_map.index(ii)]['module_name']] = \
                    LUXModuleWrapperPyMod.LUXModuleWrapper(evt_file, smart_rq_path.get_evt_path(),  '.' + rq_file, smart_rq_path.get_rq_path(), ii, data_processing_xml_path, iq_xml_path)
            end_module_time = time.time()
            print 'MODULE TIMING: Module %s completed on file %s in %.2f seconds at unix time %d' \
            % (module_settings[mod_num_map.index(ii)]['module_name'], rq_file, end_module_time - start_module_time, end_module_time)
        
        # copy files back out of RAM disk if required
        smart_rq_path.clean_up('.' + rq_file, evt_file)
        
        # unhide files
        hidden_gz = smart_rq_path.get_final_rq_path() + '/.' + rq_file + '.gz'
        final_gz = smart_rq_path.get_final_rq_path() + '/' + rq_file + '.gz'
        if os.path.isfile(hidden_gz):
            os.rename(hidden_gz, final_gz)
            os.chmod(final_gz, 0o775)

        hidden_mat = smart_rq_path.get_final_rq_path() + '/matfiles/.' + rq_file + '.mat'
        final_mat = smart_rq_path.get_final_rq_path() + '/matfiles/' + rq_file + '.mat'
        if  os.path.isfile(hidden_mat):
            os.rename(hidden_mat, final_mat)
            os.chmod(final_mat, 0o775)
        
        hidden_root = smart_rq_path.get_final_rq_path() + '/rootfiles/.' + rq_file + '.root'
        final_root = smart_rq_path.get_final_rq_path() + '/rootfiles/' + rq_file + '.root'
        if os.path.isfile(hidden_root):
            os.rename(hidden_root, final_root)
            os.chmod(final_root, 0o775)
        
        hidden_hdf5 = smart_rq_path.get_final_rq_path() + '/hdf5/.' + rq_file + '.hdf5'
        final_hdf5 = smart_rq_path.get_final_rq_path() + '/hdf5/' + rq_file + '.hdf5'
        if os.path.isfile(hidden_hdf5):
            os.rename(hidden_hdf5, final_hdf5)
            os.chmod(final_hdf5, 0o775)
        
        # delete source evt if specified
        if global_settings.has_key('delete_source_evts'):
            if int(global_settings['delete_source_evts']):
                try:
                    os.unlink(evt_dir + evt_file_in)
                except:
                    print "WARNING: Couldn't delete source evt file " + evt_dir + '/' + evt_file_in
        
        end_file_time = time.time()
        print 'FILE TIMING: %s fully processed in %.2f seconds at unix time %d\n' % (rq_file, end_file_time - start_file_time, end_file_time)

    # delete temporary xml files
    os.remove(data_processing_xml_path)
    os.remove(iq_xml_path)

    return output_logs


def rand_str_generator(size=8, chars=string.ascii_uppercase + string.digits + string.ascii_lowercase):
    return ''.join(random.choice(chars) for x in range(size))
