# RQWriter_PyMod
# Tomasz Biesiadzinski
# 07/22/2014

#
# 20140722 TPB - Created
#

import numpy as np			# read/write binary, numerical operations
# the follwoing module is used to control endianness and define the numpy
# read variables with a convenient naming convention.
from binary_read_utilities_PyMod import dt	# endianness handling class
from binary_read_utilities_PyMod import UnknownDataTypeError
# the following function constructs XML string from a dictionary
from dict_to_xmlstring_PyMod import recursive_xml_string_constructor
import time					# for timestamp for temporary file naming
import os					# change the tempoaray file name to the true one
from collections import OrderedDict as odict	# ordered keys in a dictionary

def RQWriter(filename, rq_dict, file_dict, settings_dict, \
	livetime_dict, rq_tags = None, allow_not_per_event_rq=False):
	"""
	Write a binary RQ file using various input dictionaries.
	
	Note that by default we assume that RQ array esists for each event if
	the first dimmension of the numpy array provided is number_of_events.
	If there is only 1 dimmension and it is of length number_of_events the
	the code will assume that the RQ is a single value per each event.
	
	_________________________
	Options:
		filename
			Name of the RQ file to write
		rq_dict
			A dictionary containing the reduced quantities (RQs) to write.
		file_dict
			File header dictionary to write.
		settings_dict
			A dictionary containing the DAQ and event builder settings to write.
		livetime_dict
			A dictionary containing the livetime data. 
	
	Optional:
		rq_tags = None
			A list or array with the RQ names in the order desired. If not
			provided the RQ names in the binary file will be sorted
			alphabetically.
		
		allow_not_per_event_rq = False
			Normally we expect an RQ per each event. However, by settings this
			to True, arbitrarly dimensional RQs will be writter before the 
			normal ones.
	
	_________________________
	Returns:
		Status. Equals to 1 if write was successful. 
		A file named "filename" is written out. 
	
	"""
	# Open file for writing. Give it a different filename initially so that 
	# we retain the old RQ if this code crashes.
	# 13 character (YYMMDDTHHMM) time stamp for version control
	tarr = time.localtime()
	timestamp = "%02d%02d%02dT%02d%02d%02d" \
		%(tarr[0]-2000,tarr[1],tarr[2],tarr[3],tarr[4],tarr[5])
	f = open(filename+'.%s' %(timestamp),'wb')
	# write the endian
	np.array(0x01020304, dtype=np.uint32).tofile(f)
	#
	# Get the number of events in the file.
	events_number = file_dict['nb_evts_in_file'][0]
	#
	############################################################################
	# Settings
	#
	# There may be an issue with the settings_dict. Nominally, the settings XML
	# does not have a single, all-encomapssing tag but the dictionary may since
	# xml2dict (python) requires one to work. Deal with this here.
	initags = settings_dict.keys()	# Get the settings tags.
	if len(initags) == 1 and initags[0] == 'settings':
		settings_dict_use = settings_dict['settings']
	else: settings_dict_use = settings_dict
	# convert the settings dictionary to a string
	settings_string = recursive_xml_string_constructor(settings_dict_use)
	# get the size of the settings string
	settings_size = len(settings_string)
	# write the setting size to file
	np.array(settings_size, dtype=np.uint32).tofile(f)
	# write the settings string to file
	f.write(settings_string)
	############################################################################
	# File description header
	#
	# Write the livetime header sting
	new_file_header_string  = ''
	filekeys = file_dict.keys()
	for i in range(len(filekeys)):
		# get the shape of this file. Allow for either a numpy array or list
		# of arrays
		if isinstance(file_dict[filekeys[i]], np.ndarray):
			# get the variable data type
			variable_type, bytenum = \
				dt.get_rq_vartype(file_dict[filekeys[i]].dtype.name)
			# the file is stored in a numpy array
			fileshapeall = file_dict[filekeys[i]].shape
		else:
			# get the variable data type
			variable_type, bytenum = \
				dt.get_rq_vartype(file_dict[filekeys[i]][0].dtype.name)
			# The file is stored in a list of numpy arrays
			fileshapeall = (len(file_dict[filekeys[i]]),) + \
				file_dict[filekeys[i]][0].shape
		# Strings are different. They will be written as characters
		# and we need their length. 
		if variable_type == 'char': fileshapeall = fileshapeall + (bytenum,)
		fileshape = fileshapeall
		if len(fileshape) > 1 and fileshape[0] == 1:
			# Useless first dimmension. Drop it.
			fileshape = fileshape[1:]
		# update the variable name and type in the header
		new_file_header_string += filekeys[i] + ";" + variable_type + ";"
		# add the dimmension to the header string.
		for j in range(len(fileshape)): new_file_header_string += \
			"%d," %(fileshape[j])
		# truncate the last coma
		new_file_header_string = new_file_header_string[:-1]
		# last character
		new_file_header_string += ";"
	# Get the size of the header.
	new_file_header_size = len(new_file_header_string)
	# write the new file header size
	np.array([new_file_header_size], dtype=np.uint16).tofile(f)
	# write the new header header
	f.write(new_file_header_string)
	# Write the number of lines. It is 1 as far as I know. I don't know what
	# this is supposed to mean since binary files don't have line breaks. 
	np.array([1], dtype=np.int32).tofile(f)
	# Write the header.
	for skey in filekeys: file_dict[skey].tofile(f)
	############################################################################
	# RQs
	#
	######### Test mode
	###newrq = {'a':np.arange(events_number*10).reshape((events_number,10)),
	###'b':np.arange(events_number)+100.,
	###'c':[np.arange(10)+200 for j in range(events_number)],
	###'d':np.array(np.arange(events_number)),
	###'e':np.tile(['afafasfasfa'], events_number), 
	###'f':np.tile(['tge'], events_number*10).reshape((events_number,10)),
	###'g':np.array([4]), 
	###'h':np.array(5), }
	###for skey in newrq.keys(): rq_dict[skey] = newrq[skey]
	######### 
	# get a list of the RQs from the RQ ditionary
	rqkeys = rq_dict.keys()
	# note the mode with which we will save the data
	write_mode = odict()
	# The possible write_dims are:
	# 'array per event', 'value per event', 'event-length array',
	# 'arbitrary array', 'single value'
	#
	# The following list will store the names of the RQs that do not have the
	# event number as the first dimmension. These cannot be split by events
	# when RQs are writter. 
	non_evt_rq = []
	for i in range(len(rqkeys)):
		# get the shape of this RQ 
		if isinstance(rq_dict[rqkeys[i]], np.ndarray):
			# get the RQ data type
			variable_type, bytenum = \
				dt.get_rq_vartype(rq_dict[rqkeys[i]].dtype.name)
			# the RQ is stored in a numpy array
			rqshapeall = rq_dict[rqkeys[i]].shape
		else:
			# get the RQ data type
			variable_type, bytenum = \
				dt.get_rq_vartype(rq_dict[rqkeys[i]][0].dtype.name)
			# The RQ is stored in a list of numpy arrays
			rqshapeall = (len(rq_dict[rqkeys[i]]),) + \
				rq_dict[rqkeys[i]][0].shape
		# Strings are different. They will be written as characters and we
		# need their length. 
		if variable_type == 'char': rqshapeall = rqshapeall + (bytenum,)
		#
		# get the number of dimmensions
		dimnum = len(rqshapeall)
		if dimnum > 1 and rqshapeall[0] == events_number:
			# Most commong case. There are multiple values per each event
			rqshape = rqshapeall[1:]
			write_mode_str = 'array per event'
		elif dimnum == 1 and rqshapeall[0] == events_number:
			# Single value per event
			rqshape = (1,)
			write_mode_str = 'value per event'
		elif dimnum == 0:
			# A single value
			rqshape = (1,)
			write_mode_str = 'single value'
			non_evt_rq.append(rqkeys[i])
		else:
			# just a single RQ of an arbitrary shape. 
			rqshape = rqshapeall
			write_mode_str= 'arbitrary array'
			non_evt_rq.append(rqkeys[i])
		write_mode[rqkeys[i]] = (write_mode_str, variable_type, rqshape,bytenum)
	# sort the RQs by name if their order wasn't provided.
	if rq_tags == None:
		# The RQ tags were not provided. Get them from the RQ keys. 
		rq_tags = rqkeys
		# if the dictionary is not ordered, sort the RQs alphabetically
		if not isinstance(rq_dict, odict):
			rq_tags.sort()
	#
	# We have to deal with the case of the special RQs that have wrong
	# dimmensions
	if allow_not_per_event_rq and len(non_evt_rq) > 0:
		# don't bother with this if all of the RQs are "good" dimmensionally
		rq_tags_temp = list(rq_tags)	# get a copy of the current keys
		rq_tags = [None for i in range(len(non_evt_rq))]	# placeholder 
		i=0
		for key in rq_tags_temp:	# Step through the keys.
			if key in non_evt_rq:	# Place the special RQ tags in the
				# beggining but in the same order as they appeared in the
				# rq_tags in case the user already anticipated this problem
				# and aranged the keys properly.
				rq_tags[i] = key
				i+=1
			else: rq_tags.append(key)	# Store the RQ. 
	elif (not allow_not_per_event_rq) and len(non_evt_rq) > 0:
		# throw an exception
		rqtypes = []
		rqshapes = []
		rqdtypes = []
		for i in range(len(non_evt_rq)):
			rqtypes.append(write_mode[non_evt_rq[i]][0])
			rqdtypes.append(write_mode[non_evt_rq[i]][1])
			rqshapes.append(write_mode[non_evt_rq[i]][2])
		raise(RQNotPerEvent(non_evt_rq,rqtypes,rqdtypes,rqshapes, f))
	#
	# Write the RQ header sting
	new_rq_header_string  = ''
	for i in range(len(rq_tags)):
		# update the variable name and type in the header
		variable_type = write_mode[rq_tags[i]][1]
		new_rq_header_string += rq_tags[i] + ";" + variable_type + ";"
		# add the dimmension to the header string.
		rqshape = write_mode[rq_tags[i]][2]
		for j in range(len(rqshape)): new_rq_header_string += \
			"%d," %(rqshape[j])
		# truncate the last coma
		new_rq_header_string = new_rq_header_string[:-1]
		# last character
		new_rq_header_string += ";"
	#
	# Get the size of the RQs.
	new_rq_header_size = len(new_rq_header_string)
	# write the new RQ header size
	np.array([new_rq_header_size], dtype=np.uint16).tofile(f)
	# write the new RQ header
	f.write(new_rq_header_string)
	# Write the number of lines. It is 1 as far as I know. I don't know what
	# this is supposed to mean since binary files don't have line breaks. 
	np.array([1], dtype=np.int32).tofile(f)
	#
	if allow_not_per_event_rq and len(non_evt_rq) > 0:
		# Write the RQs. RQs are generally divided by events. But first, in
		# case we have any RQs who's first dimmension isn't event number,
		# we'll write those first.
		for skey in non_evt_rq: rq_dict[skey].tofile(f)
	#
	# Write the RQs
	for i in range(events_number):
		for skey in rq_tags:
			if skey in non_evt_rq: continue
			rq_dict[skey][i].tofile(f)
	############################################################################
	# Livetime
	#
	# Write the livetime header sting
	new_lt_header_string  = ''
	# get the livetime dictionary keys
	ltkeys = livetime_dict.keys()
	# compose the header string
	for i in range(len(ltkeys)):
		# get the shape of this lt 
		if isinstance(livetime_dict[ltkeys[i]], np.ndarray):
			# get the livetime data type
			try: variable_type, bytenum = \
				dt.get_rq_vartype(livetime_dict[ltkeys[i]].dtype.name)
			except(UnknownDataTypeError):
				# I run into an issue with numpy. In the RQ reader using
				# struct.unpack numpy.array() doesn't convert an array of
				# ulonglong to uint64 but instead to dtype of "object". 
				if livetime_dict[ltkeys[i]].max() > 2**63:	
					# 64 bit unsigned integer, at least.
					livetime_dict[ltkeys[i]] = \
						np.array(livetime_dict[ltkeys[i]], dtype=np.uint64)
					variable_type, bytenum = \
						dt.get_rq_vartype(livetime_dict[ltkeys[i]].dtype.name)
				else:
					raise(UnknownDataTypeError(\
						livetime_dict[ltkeys[i]].dtype.name))
			# the livetime is stored in a numpy array
			ltshapeall = livetime_dict[ltkeys[i]].shape
		else:
			# get the livetime data type
			try: variable_type, bytenum = \
				dt.get_rq_vartype(livetime_dict[ltkeys[i]][0].dtype.name)
			except(UnknownDataTypeError):
				# I run into an issue with numpy. In the RQ reader using
				# struct.unpack numpy.array() doesn't convert an array of
				# ulonglong to uint64 but instead to dtype of "object".
				if livetime_dict[ltkeys[i]][0].max() > 2**63:	
					# 64 bit unsigned integer, at least.
					livetime_dict[ltkeys[i]] = \
						np.array(livetime_dict[ltkeys[i]], dtype=np.uint64)
					variable_type, bytenum = \
						dt.get_rq_vartype(livetime_dict[ltkeys[i]].dtype.name)
				else:
					raise(UnknownDataTypeError(\
						livetime_dict[ltkeys[i]].dtype.name))
			# The lt is stored in a list of numpy arrays
			ltshapeall = (len(livetime_dict[ltkeys[i]]),) + \
				livetime_dict[ltkeys[i]][0].shape
		# Strings are different. They will be written as characters and we
		# need their length. 
		if variable_type == 'char': ltshapeall = ltshapeall + (bytenum,)
		ltshape = ltshapeall
		# update the variable name and type in the header
		new_lt_header_string += ltkeys[i] + ";" + variable_type + ";"
		# add the dimmension to the header string.
		for j in range(len(ltshape)): new_lt_header_string += \
			"%d," %(ltshape[j])
		# truncate the last coma
		new_lt_header_string = new_lt_header_string[:-1]
		# last character
		new_lt_header_string += ";"
	# Get the size of the livetimes.
	new_lt_header_size = len(new_lt_header_string)
	# write the new livetime header size
	np.array([new_lt_header_size], dtype=np.uint16).tofile(f)
	# write the new lt header
	f.write(new_lt_header_string)
	# Write the number of lines. It is 1 as far as I know. I don't know what
	# this is supposed to mean since binary files don't have line breaks. 
	np.array([1], dtype=np.int32).tofile(f)
	# Write the livetimes.
	for skey in ltkeys: livetime_dict[skey].tofile(f)
	############################################################################
	# close the file
	f.close()
	if os.path.isfile(filename):
		# if the RQ file already exists, remove it first.
		os.remove(filename)
	os.rename(filename+'.%s' %(timestamp), filename)
	return 1

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
#		Exception definitions
#

class RQNotPerEvent(Exception):
	"""
	Exception raised when given RQs are not provided per each event as expected
	by the RQ readers. 
	"""
	def __str__(self):
		messgstr = "\nThe following RQs are not defined per each event:\n"
		itemnum = len(self.args[0])
		for i in range(itemnum):
			messgstr += \
				"RQ %s is of type '%s'. It has a shape of %s and data " \
				%(self.args[0][i], self.args[1][i], \
				self.args[3][i].__repr__()) + "types of %s.\n" \
				%(self.args[2][i])
		# close the RQ file gracefully (essentially saves the content)
		self.args[4].close()
		return messgstr

"""
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# testing


import numpy as np
import matplotlib.pyplot as pyp

import read_root as rr

import RQReader_PyMod as rrq		# reading in RQ files
import RQWriter_PyMod as rqw



filename = 'lux10_20131218T1713_f000001_cp09205.rq'



reload(rrq)


rq_dicto, settings_dicto, livetime_dicto, file_dicto = rrq.ReadRQFile('lux10_20131218T1713_f000001_cp09205_original.rq')


number_of_events = rq_dicto['pulse_classification'].shape[0]


reload(rqw)
rqw.RQWriter(filename, rq_dicto, file_dicto, settings_dicto, livetime_dicto, rq_tags = None, allow_not_per_event_rq=False)





froot = rr.read_root('/home/tomaszbi/lux/lux_data/tritium/rq/lux10_20131218T1713_cp09205/lux10_20131218T1713_f000001_cp09205.rq.root')


# compare the input and output
rq_dict, settings_dict, livetime_dict, file_dict = rrq.ReadRQFile('lux10_20131218T1713_f000001_cp09205.rq')

for skey in rq_dicto.keys():
	if 'string' not in rq_dicto[skey].dtype.name:
		aa=rq_dicto[skey] - rq_dict[skey]
		print aa.mean(), aa.std(ddof=1), 
	else: print '\t',
	print '\t', skey


"""


