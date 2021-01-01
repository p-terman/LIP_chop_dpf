"""
UseMCRPyMod.py

2013-03-30 - JRV - Created
"""

# import required modules
import os
import platform
from subprocess import Popen

class UseMCR:
    def __init__(self, global_settings):
        self.global_settings = global_settings
        self.mcr_path_to_use = []

        if global_settings.has_key('matlab_mcr_path'):
            self.mcr_root = global_settings['matlab_mcr_path']
            self.mcr_root_stripped = self.mcr_root.rstrip('/')
            self.mcr_path_to_use = mcr_root
            
            if global_settings.has_key('ramdisk_path'):
                self.ramdisk_path = global_settings['ramdisk_path']
                set_ramdisk_path()
                copy_mcr_to_ramdisk()
            set_mcr_cache()
            set_mcr_paths()


    def check_ramdisk(self, ramdisk_path):
        if os.access(ramdisk_path, os.W_OK):
            return True
        else:
            return False

    def set_ramdisk_path(self):
        if check_ramdisk(self.ramdisk_path):
            self.mcr_path_to_use = self.ramdisk_path + '/' + os.path.basename(self.mcr_root_stripped) + '/'
            copy_mcr_to_ramdisk()
        else:
            self.mcr_path_to_use = mcr_root


    def copy_mcr_to_ramdisk(self):
        if not os.path.isdir(self.mcr_path_to_use):
            p = Popen(['rsync', '-a', mcr_root_stripped, self.mcr_path_to_use + '/'])
            p.wait()
   

    def get_mcr_path(self):
        return self.mcr_path_to_use 
   

    def set_mcr_cache(self):
        temp_dir = get_mcr_path()
        smart_rq_path.create_directory(temp_dir + '/mcr_cache_root_lux/')
        os.environ['MCR_CACHE_ROOT'] = temp_dir + '/mcr_cache_root_lux/'

    
    def set_mcr_paths(mcr_path_to_use):
        # for 64 bit Linux
        if platform.machine() == 'x86_64' and platform.system() == 'Linux':
            if os.environ.has_key('LD_LIBRARY_PATH'):
                ld_library_path = os.environ['LD_LIBRARY_PATH']
            else:
                ld_library_path = ''
    
            mcr_jre = mcr_path_to_use + '/sys/java/jre/glnxa64/jre/lib/amd64'
            
            mcr_paths = [mcr_path_to_use + '/runtime/glnxa64', \
            mcr_path_to_use + '/bin/glnxa64', \
            mcr_path_to_use + '/sys/os/glnxa64', \
            mcr_jre + '/native_threads', \
            mcr_jre + '/server', \
            mcr_jre + '/client', \
            mcr_jre]
            
            xapplresdir = mcr_path_to_use + '/X11/app-defaults'
    
            for p in mcr_paths:
                if p not in ld_library_path:
                    ld_library_path += (':' + p)
            
            os.environ['XAPPLRESDIR'] = xapplresdir
            os.environ['LD_LIBRARY_PATH'] = ld_library_path
