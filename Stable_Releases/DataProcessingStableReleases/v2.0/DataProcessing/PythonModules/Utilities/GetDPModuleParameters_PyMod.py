# GetDPModuleParameters_PyMod
# Tomasz Biesiadzinski
# July 07, 2014

from xml2dict import xml2dict					# reading/parsing xml files
import numpy as np								# mathematical operations
from collections import OrderedDict as odict	# ordered keys in a dictionary

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# Constant definitions

data_processing_dict_key = 'data_processing_settings'
generic_module_dict_key = 'module'
module_name_dict_key = 'module_name'
parameters_dict_key = 'parameters'

generic_global_dict_key = 'global'

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

def GetDPModuleParameters(data_processing_xml_path, \
	modules, return_parsed_xml = False):
	"""
	This function reads the data processing xml file and extract the parameters
	needed by a given module.
	
	_________________________
	Options:
		data_processing_xml_path
			The path to the data processing xml file.
		modules
			Name of the module settings in the xml file. Or a list of such
			module names.
	
	Optional
		return_parsed_xml = False
			If set to True, in addition to the normal return dictionary we
			will return the full xml dictionary than can be accessed arbitrarly.
	
	_________________________
	Returns:
		Two dictionaries will be returned. The second one will contain the
		entries within the <global> tag of the DP xml file. The first dictionary
		will be as follows:
		- If a single module name was provided, the function will return a
		  dictionary with the various parameters keyed to their names. By
		  default the function attempts to convert them to floats. If a value
		  cannot be converted to float it is returned as a string.
		- If multiple module names were provided, the return dictionary will
		  contain multiple dictionaries described above, each keyed to the
		  module name.
		
		Note that if return_parsed_xml is set to True, a third dictionary, 
		containing the entire DP xml file will be returned.
	"""
	# make sure module_names are a list
	if not isinstance(modules, list) and \
		not isinstance(modules, np.ndarray) and \
		not isinstance(modules, tuple):
		modules = [modules]
	#
	# Read in the xml file
	dp_settings_xml = xml2dict(data_processing_xml_path)
	# Get the list fo all module names
	module_names = [tempdict[module_name_dict_key] for tempdict in \
		dp_settings_xml[data_processing_dict_key ][generic_module_dict_key]]
	# Convert to a numpy array of strings so that we can search it easily
	module_names = np.array(module_names)
	# initialize the return dictionary
	retpars = odict()
	# step through the moduled desired by the user
	for module_name in modules:
		# locate the module settings index and check that it is found.
		index_module_arr = np.where(module_names == module_name)[0]
		if index_module_arr.size == 0:
			# The settings are missing.
			raise(DPSettingsMissingError(module_name, module_names))
		if index_module_arr.size > 1:
			# There are more than one version of the settings.
			raise(DPSettingsMultipleError(module_name, module_names))
		#
		# we can now get the index module
		index_module = index_module_arr[0]
		#
		# get the module parameters
		mymodule_settings = dp_settings_xml[data_processing_dict_key]\
			[generic_module_dict_key][index_module][parameters_dict_key]
		parameters = mymodule_settings.keys()	# available xml tags
		if len(modules) > 1:
			retpars[module_name] = odict()
			retparspointer = retpars[module_name]	# pointer to the module
		else:
			retparspointer = retpars	# pointer to the full return dictionary
		for parameter_name in parameters:
			tempstr = mymodule_settings[parameter_name]
			try:	# convert to a float
				retparspointer[parameter_name] = float(tempstr)
			except(ValueError):	# return the string instead
				retparspointer[parameter_name] = tempstr
	#
	# get the global dictionary
	if not dp_settings_xml[data_processing_dict_key].has_key( \
		generic_global_dict_key):
		raise(DPSettingsMissingError(generic_global, \
			dp_settings_xml[data_processing_dict_key].keys()))
	retparglob = dp_settings_xml[data_processing_dict_key]\
		[generic_global_dict_key]
	# Attempt to convert to floats
	for skey in retparglob.keys():
		if not isinstance(retparglob[skey], str): continue	# can't convert
		try: retparglob[skey] = float(retparglob[skey])	# try to convert
		except(ValueError):	pass	# don't do anything
	#
	# return the parameters
	if return_parsed_xml:	# also return the parsed dictionary for direct access
		return retpars, retparglob, dp_settings_xml[data_processing_dict_key]
	else: return retpars, retparglob


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
#		Exception definitions
#

class DPSettingsMissingError(Exception):
	"""
	Raised when the module settings are not found in the DP settings xml
	"""
	def __str__(self):
		messgstr = "\nSettings for %s are missing." %(self.args[0])
		messgstr = "Settings are only available for: "
		for item in self.args[1]: messgstr += "%s, " %(item)
		return messgstr

class DPSettingsMultipleError(Exception):
	"""
	Raised when the module settings are found multiple times.
	"""
	def __str__(self):
		messgstr = "\nSettings for %s occur multiple times in: " %(self.args[0])
		for item in self.args[1]: messgstr += "%s, " %(item)
		return messgstr


