"""
Install script

The basic code one needs to add another setup module is the following:

@install_action
class empty:
  def __init__(self):
    pass
  def pre_install_action(self):
    pass
  def install_action(self):
    pass

Change "empty" to a descriptive name and fill the in the three section:
  __init__: is good for getting/setting variables from the global dict
  pre_install_action: can be used for confirmation or extra info.
  install_action: is where things should actually be written.
2013-02-11 - Initial submission
"""

import os
import re

luxsetpythonpaths_path = "DataProcessingFramework/Utilities/LUXSetPythonPaths.py"

function_dict = {}
install_action_list = []
global_dict = {}

def install_abort():
    os.sys.exit("\n** Installation aborted! No files have been written **")

def raw_input_with_default(prompt, default):
  resp = raw_input(prompt)
  if resp in ['q', 'Q', "exit", "EXIT", "Exit", "quit", "Quit", "QUIT"]:
    install_abort()
  if resp == "":
    return default
  return resp

def is_blank(resp):
  if resp in ["Y", "y", "Yes", "yes", "YES"]:
    return True
  return False

def is_yes(resp):
  resp = resp.replace('[', '').replace(']','')
  if resp in ["Y", "y", "Yes", "yes", "YES"]:
    return True
  return False


def preparatory_action(cl):
  install_action_list.append(cl())
  return cl

def install_action(cl):
  instance = cl()
  instance.global_dict = global_dict
  instance.active = False
  install_action_list.append(instance)
  return cl

@install_action
class check_for_mysqldb:
  def __init__(self):
    pass
  def pre_install_action(self):
    try:
      import MySQLdb
    except ImportError:
      print "It seems you do not have MySQL db installed. This is required"
      print "for LUG/Slow Control database communications. Try:"
      print "  ", "* apt-get install mysql-python (Debian, Ubuntu, etc.)"
      print "  ", "* yum install mysql-python (Red Hat, Cent OS, Fedora, etc.)"
      print "  ", "* Install mysql (www.mysql.com/) ",
      print "and try pip install MySQL-python (Mac)"
      print ""
      print """If when you run the command"""
      print '\t' + 'python -c "import MySQLdb"'
      print "it executes and quits without erroring, then you're ready"
      install_abort()
  def install_action(self):
    pass

@install_action
class get_dp_path:
  def __init__(self):
    self.prompt = "Enter full path to data processing directory [CWD]:"
  def pre_install_action(self):
    default_dp_path = os.getcwd()
    dp_path = raw_input_with_default(self.prompt, default_dp_path)
    self.global_dict["dp_path"] = dp_path
  def install_action(self):
    pass

@install_action
class create_luxsetpythonpaths:
  def __init__(self):
    pass
  def pre_install_action(self):
    self.active = True
  def install_action(self):
    print "Creating..."
    print "\tDataProcessing/" + luxsetpythonpaths_path
    file_path = self.global_dict['dp_path'] + "/" + luxsetpythonpaths_path
    self.global_dict['setpythonpaths'] = file_path
    f = open(file_path, "w")

    f.write("""
import sys
import os
from os.path import expanduser

dp_path = '%s'

# all paths need for python framework need to be added here
relative_paths = ['DataProcessingFramework',\
      '/',\
              'PythonModules/']

# put all full paths into list
paths = [dp_path + s for s in relative_paths]

# add all paths to python path
for p in paths:
  for root, dirs, files in os.walk(p):
    if not os.path.basename(root).startswith('.'):
      if root not in sys.path:
        sys.path.append(root)
        """ % self.global_dict['dp_path'])

  
@install_action
class make_cpp_modules:
  def __init__(self):
    self.prompt = "Do you want to compile the C++ modules? [Y]/n "
  def pre_install_action(self):
    resp = raw_input_with_default(self.prompt, "Y")
    self.active = True if is_yes(resp) else False
  def install_action(self):
    if self.active:
      bst_file = open("/tmp/boostcheck.cc","w")
      bst_file.write("#include <boost/property_tree/ptree.hpp>\nint main(){return 0;}")
      bst_file.close()
      ret = os.system("g++ $BOOST /tmp/boostcheck.cc -o /tmp/boostcheck")
      if ret != 0:
        print "** I can't find your boost libraries (www.boost.org). Check"
        print "   your installation and try again."
        install_abort()
      ret = os.system("g++ $BOOST /tmp/boostcheck.cc -o /tmp/boostcheck")
      if ret != 0:
        print "** Your version of g++ doesn't support c++11. Time for an update!"
        print "** MAC USERS: You may use clang (Apple's LLVM front end) instead"
        print "   of g++. Try replacing g++ with clang++ in all makefiles until"
        print "   this is automated."
        install_abort()
      os.system("make -C " + self.global_dict['dp_path'] + "/CppModules/ cpp")
      ret = os.system("rm /tmp/boostcheck.cc /tmp/boostcheck")

@install_action
class make_root_modules:
  def __init__(self):
    self.prompt = "Do you want to compile the ROOT modules? [Y]/n "
  def pre_install_action(self):
    resp = raw_input_with_default(self.prompt, "Y")
    self.active = True if is_yes(resp) else False
  def install_action(self):
    if self.active:
      os.system("make -C " + self.global_dict['dp_path'] + "/CppModules/ root")

@install_action
class setup_matlab_path:
  def __init__(self):
    self.prompt = "Do you want to setup the matlab path? You must do this\n"
    self.prompt += "\twhenever directory structures change. [Y]/n "
  def pre_install_action(self):
    resp = raw_input_with_default(self.prompt, "Y")
    self.active = True if is_yes(resp) else False
    if self.active:
      prompt = "Specify the location of your matlab binary. (e.g. /usr/local/bin/matlab) "
      self.matlab_path = raw_input(prompt)
      self.global_dict['matlab_path'] = self.matlab_path
  def install_action(self):
    if self.active:
      dp_path = self.global_dict['dp_path']
      comm = dp_path
      comm += "/MatlabModules/Utilities/setup_matlab_path.sh"
      comm += " " + dp_path + " " +  self.matlab_path
      os.system(comm)


@install_action
class update_run_dp_framework:
  def __init__(self):
    self.prompt = "Do you want to update run_dp_framework.py with your settings? y/[N]? "
  def pre_install_action(self):
    resp = raw_input_with_default(self.prompt, "N")
    self.active = True if is_yes(resp) else False
    if not self.active:
      return
    # First set the default xml dp settings
    prompt = "What xml settings file do you want to use (enter path relative"
    prompt += " to DataProcesssing directory)? "
    prompt += "[default_data_processing_settings.xml] "
    def_xml = 'default_data_processing_settings.xml'
    xml_path = raw_input_with_default(prompt, def_xml)
    dp_path = self.global_dict['dp_path']
    xml_path = dp_path + '/' + xml_path
    #print "\n\t** Ensure paths are correct set in your xml file!"
    #print "\n\t** This install will not do that for you, sadly."
    self.global_dict['xml_path'] = xml_path
  def install_action(self):
    dp_path = self.global_dict['dp_path']
    if not self.active:
      return
    run_dp_fw = open(dp_path + "/Examples/run_dp_framework.py").readlines()
    for i, line in enumerate(run_dp_fw):
      if "execfile" in line:
        line = "execfile('" + self.global_dict['setpythonpaths'] + "')\n"
        run_dp_fw[i] = line
      if "dp_settings_xml_path = " in line:
        line = "dp_settings_xml_path = '" + self.global_dict['xml_path'] + "'\n"
        run_dp_fw[i] = line
    open(dp_path + "/Examples/run_dp_framework.py", 'w').writelines(run_dp_fw)
    print "run_dp_framework.py edited..."

@install_action
class update_dp_settings:
  def __init__(self):
    # We don't look for the xml path here becuase __init__ is run before we
    # get any user input.
    pass
  def pre_install_action(self):
    if not self.global_dict.has_key('xml_path'):
      self.active = False
      return
    self.prompt = "Does your dp settings xml file need updating, too? [Y]/n? "
    resp = raw_input_with_default(self.prompt, "Y")
    self.active = True if is_yes(resp) else False
    if not self.active:
      return
  def install_action(self):
    if not self.active:
      return
    xml_path = self.global_dict['xml_path']
    xml = open(xml_path).read()
    try:
      matlab_path = self.global_dict['matlab_path']
      xml = re.sub(r'(<matlab_path.*>).*(</matlab_path>)', r'\1%s\2' % matlab_path, xml)
    except:
      pass
    open(xml_path, 'w').write(xml)
    print xml_path + " edited..."

@install_action
class mex_transparent_rubiks_cube:
  def __init__(self):
    self.prompt = "Do you want to mex TransparentRubiksCube? [Y]/n? "
  def pre_install_action(self):
    resp = raw_input_with_default(self.prompt, "Y")
    self.active = True if is_yes(resp) else False
    if not self.active:
      return
    prompt = "Specify full path to mex if needed [mex]: "
    self.mex_path = raw_input_with_default(prompt, "mex")
  def install_action(self):
    dp_path = self.global_dict['dp_path']
    comm = "cd " + dp_path + "/"
    comm += "MatlabModules/PulseFinder_TransparentRubiksCube/ "
    comm += " && "
    comm += self.mex_path + " perusePeeks.c"
    os.system(comm)

# Set up dictionary of methods and 
if __name__ == "__main__":
  
  print "Welcome to the Data Processing Framework installer!"
  print ""
  print "When prompted, the default response is in [brackets] and you can just"
  print "press enter to accept it. You can type exit at any prompt to quit."
  print ""
  print "_"*76
  print ""

  #create_luxsetpythonpaths("/Users/mwoods/School/Research/LUX/DataProcessing")
  #print "Install dict = ", install_action_list
  #c = create_luxsetpythonpaths()
  #print c
  #print "Install dict = ", install_action_list

  try:
    for func in install_action_list:
      func.pre_install_action()
  except KeyboardInterrupt:
    install_abort()
      
  for func in install_action_list:
    if func.active:
      print "v"* 79
      print " "*10, "Running", func.__class__.__name__
      func.install_action()
      print "^"* 79
      print ""
      print ""

  # Set functions that need to run.
