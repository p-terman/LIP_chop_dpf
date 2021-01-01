"""
LUXModuleWrapperPyMod.py


The LUXModuleWrapper function should be used around every
module, no exceptions. This function takes in the
data_processing_settings.xml and the desired module name,
makes all required lug queries, runs the module, and reports
the status.

2012-11-12 - JRV - Created
2014-07-23 - CHF - Modified RunPythonModule to now use python module name from 'relative_path'
                    instead of hard-coded 'PyMod'
"""

# import required python modules here

from subprocess import Popen, PIPE, STDOUT
from str2numPyMod import str2num
import ReportLUXcodePathPyMod
from xml2dict import xml2dict
import os
import platform
import re

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

def LUXModuleWrapper(evt_file_name, evt_file_path, rq_file_name, rq_file_path, run_module_number, dp_settings_xml_path, lug_iqs_xml_path):
    """
    Inputs:
        evt_file_name - name of the evt file
        evt_file_path - path to the directory containing the evt file
        rq_file_name - name of the rq file
        rq_file_path - path to the directory containing the rq file
        run_module_number - the desired module number as defined in the xml
        dp_settings_xml_path - the full path to the data_processing_settings.xml file
        lug_iqs_xml_path - the full path to the lug_iqs.xml file
    """

    module_output = None
    #global_settings = GetGlobalSettings(dp_settings_xml_path)
    mod, mod_num_map = GetModuleSettings(dp_settings_xml_path)

    lang = mod[mod_num_map.index(run_module_number)]['language']
    if (lang == 'matlab'):
        # run matlab module
        module_output = RunMatlabModule(evt_file_name, evt_file_path, rq_file_name, rq_file_path, run_module_number, dp_settings_xml_path, lug_iqs_xml_path)
    elif (lang == 'matlab-compiled'):
        # run compiled matlab module
        module_output = RunMatlabCompiledModule(evt_file_name, evt_file_path, rq_file_name, rq_file_path, run_module_number, dp_settings_xml_path, lug_iqs_xml_path)
    elif (lang == 'root'):
        # run ROOT module
        module_output = RunRootModule(evt_file_name, evt_file_path, rq_file_name, rq_file_path, run_module_number, dp_settings_xml_path, lug_iqs_xml_path)
    elif (lang == 'python'):
        # run python module
        module_output = RunPythonModule(evt_file_name, evt_file_path, rq_file_name, rq_file_path, run_module_number, dp_settings_xml_path, lug_iqs_xml_path)
    elif (lang == 'julia'):
        # run julia module
        module_output = RunJuliaModule(evt_file_name, evt_file_path, rq_file_name, rq_file_path, run_module_number, dp_settings_xml_path, lug_iqs_xml_path)
    elif (lang == 'cpp'):
        # run cpp module
        module_output = RunCppModule(evt_file_name, evt_file_path, rq_file_name, rq_file_path, run_module_number, dp_settings_xml_path, lug_iqs_xml_path)

    return module_output


def GetGlobalSettings(dp_settings_xml_path):
    """ This function is used to load the global DP settings """

    global_settings = None

    # load the xml file
    xml = xml2dict(dp_settings_xml_path)

    # write global settings into a dict
    gs_keys = list()
    gs_values = list()
    for key in xml['data_processing_settings']['global'].keys():
        gs_keys.append(key)
        gs_values.append(xml['data_processing_settings']['global'][key])
    global_settings = dict(zip(gs_keys, gs_values))

    return global_settings


def GetModuleSettings(dp_settings_xml_path):
    """
    This function is used load the module settings.

    The variable nod_num_map maps betweeen the run_order number and
    the position in the xml file. This is used as follows:
        mod_num_map.index(run_module_number)

    """

    mod = None
    mod_num_map = None

    # load the xml file
    xml = xml2dict(dp_settings_xml_path)

    # If there is only one module specified, xml2dict doesn't return a list of
    # module settings. Let's toss it in to a list of just itself.
    if(type(xml['data_processing_settings']['module']) != list):
        xml['data_processing_settings']['module'] = [xml['data_processing_settings']['module']]
    # write module settings into list of dicts
    modules = list()
    for module in xml['data_processing_settings']['module']:
        mod_keys = list()
        mod_values = list()
        for key in module.keys():
            mod_keys.append(key)
            mod_values.append(module[key])
        modules.append(dict(zip(mod_keys, mod_values)))
    mod = modules

    # create mapping by module number
    mod_num_map = list()
    for m in mod:
        mod_num_map.append(str2num(m['run_order']))

    return mod, mod_num_map


def RunMatlabModule(evt_file_name, evt_file_path, rq_file_name, rq_file_path, run_module_number, dp_settings_xml_path, lug_iqs_xml_path):
    """ 
    This subfunction is used to run an individual Matlab module
    
    2013-03-05 - JRV - Changed all paths relative to DataProcessing/; fixed DP bug #19
    """

    module_output = None

    # get settings from xml file
    global_settings = GetGlobalSettings(dp_settings_xml_path)
    mod, mod_num_map = GetModuleSettings(dp_settings_xml_path)

    # determine debug mode, which disables the try/catch statement.
    try:
        debug_mode = int(global_settings['debug_mode'])
        print "Debug mode is", str(debug_mode)
    except:
        print "Debug modes disabled."
        debug_mode = 0
    # setup paths
    matlab_path = global_settings['matlab_path']
    #module_path = matlab_path + mod[mod_num_map.index(run_module_number)].relative_path
    module_name = mod[mod_num_map.index(run_module_number)]['module_name']
    dp_path = ReportLUXcodePathPyMod.ReportDataProcessingPath()

    # setup matlab command
    saved_matlab_path_file = dp_path + '/MatlabModules/Utilities/matlab_path.mat'  # this will have to change to be external
    command = list()
    if debug_mode == 0:
        command.append('try, ')
    command.append("cd('" + dp_path + "/MatlabModules/Utilities/');")
    command.append("load_saved_matlab_path('" + saved_matlab_path_file + "');")  # this goes at the beginning to set path
    command.append(module_name + "('%s','%s','%s','%s','%s','%s');" \
            % (evt_file_name, evt_file_path, rq_file_name, rq_file_path, dp_settings_xml_path, lug_iqs_xml_path))  # this line is the module call
    command.append('quit; ')  # this goes at the end to exit matlab
    if debug_mode == 0:
        command.append('catch ME, ')
        command.append("delete('%s/%s');" % (rq_file_path, rq_file_name))
        command.append('quit,')
        command.append('end')
    matlab_command = ''.join(command)  # join all commands in single string
    
    # run module in matlab
    #p = Popen([matlab_path, '-nojvm', '-nodesktop', '-nosplash', '-r', matlab_command], stdout=PIPE, stderr=STDOUT)
    p = Popen([matlab_path, '-nojvm', '-nodesktop', '-nosplash', '-r', matlab_command])

    # wait while matlab runs, then return module output
    p.wait()

    module_output = p.communicate()[0]

    return module_output


def RunMatlabCompiledModule(evt_file_name, evt_file_path, rq_file_name, rq_file_path, run_module_number, dp_settings_xml_path, lug_iqs_xml_path):
    """
    This subfunction is used to run an individual compiled Matlab module
    
    2013-04-29 - JRV - Created
    """

    module_output = None

    # get settings from xml file
    global_settings = GetGlobalSettings(dp_settings_xml_path)
    mod, mod_num_map = GetModuleSettings(dp_settings_xml_path)
    module_settings = mod[mod_num_map.index(run_module_number)]
    # setup paths
    matlab_mcr_path = global_settings['matlab_mcr_path']
    
    module_name = mod[mod_num_map.index(run_module_number)]['module_name']
    dp_path = ReportLUXcodePathPyMod.ReportDataProcessingPath()
    bin_path = module_settings['relative_path']

    command = list()
    command.append(dp_path + bin_path)
    command.append(evt_file_name)
    command.append(evt_file_path)
    command.append(rq_file_name)
    command.append(rq_file_path)
    command.append(dp_settings_xml_path)
    command.append(lug_iqs_xml_path)

    # TODO: This should use a temp directory if no ramdisk
    if global_settings.has_key('ramdisk_path'): 
        os.environ['MCR_CACHE_ROOT'] = global_settings['ramdisk_path'] + '/mcr_cache_root_lux/'
    
    set_mcr_paths(matlab_mcr_path)
    my_env = os.environ.copy()

    p = Popen(command, env=my_env)
    p.wait()
    module_output = p.communicate()[0]

    return module_output


def RunRootModule(evt_file_name, evt_file_path, rq_file_name, rq_file_path, run_module_number, dp_settings_xml_path, lug_iqs_xml_path):
    """ This subfunction is used to run an individual Root module """
    module_output = None

    # get settings from xml file
    #global_settings = GetGlobalSettings(dp_settings_xml_path)
    mod, mod_num_map = GetModuleSettings(dp_settings_xml_path)

    #print global_settings
    module_settings = mod[mod_num_map.index(run_module_number)]

    #root_path = global_settings['LUXcode_path']  # root as in fundamental.
    root_path = ReportLUXcodePathPyMod.ReportDataProcessingPath()
    bin_path = module_settings['relative_path']

    commands = [root_path + bin_path]
    commands.append(evt_file_name)
    commands.append(evt_file_path)
    commands.append(rq_file_name)
    commands.append(rq_file_path)
    commands.append(str(run_module_number))
    commands.append(dp_settings_xml_path)
    commands.append(lug_iqs_xml_path)
    #p = Popen(commands, stdout=PIPE, stderr=STDOUT)
    print commands
    p = Popen(commands)
    p.wait()
    module_output = p.communicate()[0]

    return module_output


def RunPythonModule(evt_file_name, evt_file_path, rq_file_name, rq_file_path, run_module_number, dp_settings_xml_path, lug_iqs_xml_path):
    """
    This subfunction is used to run an individual Python module
    
    2013-03-06 - JRV - Created
    2014-07-23 - CHF - Removed hard-coded 'PyMod' suffix. Now using 'relative_path' from XML module settings
                        It does a regular expression check for modulename_***.py in relative_path, 
                        then returns the module name from there.
                        - Fixed regex code. Module needs to end with PyMod
    """

    module_output = None
    
    # get settings from xml file
    global_settings = GetGlobalSettings(dp_settings_xml_path)
    mod, mod_num_map = GetModuleSettings(dp_settings_xml_path)
    module_name = mod[mod_num_map.index(run_module_number)]['module_name']
    module_relative_path = mod[mod_num_map.index(run_module_number)]['relative_path']
    python_module_name = re.findall('%s_?PyMod\.py$' % module_name,module_relative_path)[0].rstrip('.py')

    try:
        # import DP module as needed
        exec('from %s import %s as dp_module' % (python_module_name, module_name)) 
        # run specified module
        module_output = dp_module(evt_file_name, evt_file_path, rq_file_name, rq_file_path, run_module_number, dp_settings_xml_path, lug_iqs_xml_path)
    except:
        print "ERROR: Module %s execution failed on %s" % (module_name, evt_file_name)

    return module_output


def RunCppModule(evt_file_name, evt_file_path, rq_file_name, rq_file_path, run_module_number, dp_settings_xml_path, lug_iqs_xml_path):
    """ This subfunction is used to run an individual Cpp module """

    module_output = None

    # get settings from xml file
    global_settings = GetGlobalSettings(dp_settings_xml_path)
    mod, mod_num_map = GetModuleSettings(dp_settings_xml_path)

    #print global_settings
    module_settings = mod[mod_num_map.index(run_module_number)]

    root_path = ReportLUXcodePathPyMod.ReportDataProcessingPath()
    bin_path = module_settings['relative_path']

    commands = [root_path + bin_path]
    commands.append(evt_file_name)
    commands.append(evt_file_path)
    commands.append(rq_file_name)
    commands.append(rq_file_path)
    commands.append(str(run_module_number))
    commands.append(dp_settings_xml_path)
    commands.append(lug_iqs_xml_path)
    p = Popen(commands, stdout=PIPE)

    p.wait()
    module_output = p.communicate()[0]

    return module_output


def RunJuliaModule(evt_file_name, evt_file_path, rq_file_name, rq_file_path, run_module_number, dp_settings_xml_path, lug_iqs_xml_path):
    """ This subfunction is used to run an individual Julia module """

    module_output = None

    return module_output
