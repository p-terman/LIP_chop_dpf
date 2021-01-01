
# EvtReader_PyMod.py
# Tomasz Biesiadzinski
# 03/02/2014

import os, re
import numpy as n
import glob
import xml2dict

# the follwoing module is used to control endianness and define the numpy
# read variables with a convenient naming convention.
from binary_read_utilities_PyMod import dt	# endianness handling class

#try:
#	from multiprocessing import Process, Queue
#	import numprocs	# to count the cores
#	can_multiprocess = True
#except:
#	print 'Unable to import the multiprocessing module'
#	can_multiprocess = False
# Multiprocessing is not enabled yet


# Notes:
#
# PMTs 1-122 are in TPC, PMTs 129-136 are in water, 123 is the analog
# top PMT sum, 124 is the analog bottom PMT sum, 125 is empty, 126 is
# the TTE, 127 is empty, 128 is the Trigger. 

class EvtReader:
	def __init__(self, file_or_directory, filematch = '*.evt', \
		max_channel = None,
		recalculate_baseline = True, recalculate_baseline_samples = 22,
		apply_amp_gains=False, preamp_gain = None, postamp_gain = None, 
		ns_per_sample = 10., 
		debug=False):
		"""
		Initialize the EvtReader object.
		
		
		Indecies are relative to 0. That is, the first event is indexed 0
		when you need to retrieve it regardless of what it was labeled in the
		evt file. This is an issue if you've read files not in order or
		all of the preceeding ones. You can look up the event id from
		the file and then request the corresponding event index here.
		
		_________________________
		Options:
			file_or_directory
				An evt file, list of evt files or a directory containing
				evt files.
		
		Optional:
			filematch = '*.evt'
			The format of files we're going to use. Can use unix wildcards.
			
			max_channel = None
				If given, don't read any channels above it. Note, this is
				0 indexed which means that the TPC PMTs run from 0 to 121,
				NOT 1 to 122.
			
			recalculate_baseline = True
				The baseline obtained by the DAQ and saved in the EVT file is
				biased by (very) roughly 0.5ADU. If this option is set to True,
				the baseline will be recomputed using the first
				recalculate_baseline_samples of each POD waveform.
				Note that this option specified the default behavior of the
				reader. When self.get_event_data() is called the user
				can change this option.
		
			recalculate_baseline_samples = 22
				The number of samples that will be used to recalcuate the
				baseline if recalculate_baseline is set to True.
				Note that this option specified the default behavior of the
				reader. When self.get_event_data() is called the user
				can change this option.
			
			apply_amp_gains = False
				If True, preamp_gain and postamp_gain will be applied to the
				waveforms along witn ns_per_sample
			
			preamp_gain = None
				The pre-amplifier gain. Used only if apply_amp_gain is set to
				True. If set to None, it will be read from
				the xml settings within the EVT file. 
			
			postamp_gain = None
				The post-amplifier gain. Used only if apply_amp_gain is set to
				True. If set to None, it will be read from
				the xml settings within the EVT file.
			
			ns_per_sample = 10.
				For conversion from mV to mVns. Applied only if apply_amp_gains
				is set to True.
			
			debug = False
				If set to true diagnostic info will be printed out from
				various functions.
		
		_________________________
		Returns:
			An EvtReader class object. 
		
		"""
		# for debugging
		self.debug = debug
		# read the evt file headers
		self.filenum,self.filelist,self.xmlsettings,self.evt_number_per_file, \
			self.date_and_time, self.live_time_header, self.evt_ids, \
			self.evt_gids, self.evt_file_inds = \
				prepare_evt_info(file_or_directory, filematch = filematch)
		#_________________________
		# get the number of events
		self.evt_number = self.evt_gids.size
		if self.debug:
			print 'Debugging mode'
			print '%d files' %(len(self.filelist))
			print '%d events in files' %(self.evt_number)
		#_________________________
		# Currently all of the files are closed. However once we start addressing
		# them we should keep the current file open for speed. We'll initialize
		# the variables that will keep track of which file, if any, is open and 
		# the filehandle of that file.
		self.filehandle = None
		self.currently_opened_file_index = None
		# There is an issue with unhandled exceptions. If an error occurs during
		# a file read the file pointer will not be at the event beginning
		# location. Attempting to read again from the same object will then
		# throw EVTFilePositionError. In order to clean this up the file must be
		# advanced back to its proper position
		self.fileread_in_progress = False
		# Also keep track of which event we're on right now
		self.current_evt_ind = 0
		# store the maximum channel to be read
		self.max_channel = max_channel
		#_________________________
		# store the default behavior for how the baseline will be treated.
		self.recalculate_baseline = True
		self.recalculate_baseline_samples = 22
		#_________________________
		# find out whether the trigger info needs to be read for each event.
		# If it's not present that reading those bytes will completally mess up
		# our location in the file.
		if 'read_xlm' in self.xmlsettings['daq_settings']['sis3301']['global'].keys():
			# learned about this from Surge's dat code the hard way after
			# my code failed on some bad data.
			trig_info_flag = int(self.xmlsettings['daq_settings']['sis3301']['global']['read_xlm'][0])
			if trig_info_flag == 1:
				# read in the trigger information below
				self.read_trigger_info = True
			else:
				self.read_trigger_info = False
		else:
			raise(SettingMissingError('read_xlm',
				'Number of bytes to read cannot be determined'))
		#_________________________
		# read in the preamp and postamp gains
		if apply_amp_gains:
			if preamp_gain == None:
				if 'preamp' in self.xmlsettings['daq_settings']['global']:
					self.preamp_gain = float(self.xmlsettings['daq_settings']['global']['preamp'][0])
				else:
					raise(SettingMissingError('preamp',
						'Cannot locate the pre-amplifier gain in the file.'))
			else: self.preamp_gain = preamp_gain
			if postamp_gain == None:
				if 'postamp' in self.xmlsettings['daq_settings']['global']:
					self.postamp_gain = float(self.xmlsettings['daq_settings']['global']['postamp'][0])
				else:
					raise(SettingMissingError('postamp',
						'Cannot locate the post-amplifier gain in the file.'))
			else: self.postamp_gain = postamp_gain
			# save the temporal resolution
			self.ns_per_sample = ns_per_sample
		else:
			# don't apply the gains
			self.preamp_gain, self.postamp_gain, self.ns_per_sample = 1., 1., 1.
	
	def __getitem__(self, i):
		"""
		x.__getitem__(i) <==> x[i]
		"""
		if i < 0: # support for negative indecies
			start = self.evt_number+i
		else:
			start = i
		if self.debug: print '__getitem__ of %d' %(start)
		return self.get_event_data(start=start)
	
	def __getslice__(self, i, j):
		"""
		x.__getslice__(i, j) <==> x[i:j]
		"""
		if i < 0: # support for negative indecies
			start = self.evt_number+i
		else:
			start = i
		if j < 0: # support for negative indecies
			end = self.evt_number+j
		else:
			end = j
		if self.debug: print '__getslice__ from %d to %d' %(start,end)
		return self.get_event_data(start = start, end = end)
	
	def get_event_data(self, start = None, end = None, length = 1, 
		get_all=False, pulse_pointers_only=False, max_channel=None,
		recalculate_baseline = None, recalculate_baseline_samples = None):
		"""
		Get the raw data read from an evt file. It's not completally raw
		(no binary or anything) but it's not processed much.
		
		This function performs a lot of error checking on file location. It
		does not override existing file location if the current event is
		requested. In case we are at the wrong location an exception will
		be raised since something went wrong and the user should be aware
		of it.
		
		_________________________
		Options:
			None
		
		Optional:
			start = None
				Start event index. If None we'll start from the current
				event index.
			
			end = None
				The end event index to be read. If None, only one event
				will be read.
				Note that this is inclusive. That is, event with index
				end will be read in.
			
			length = 1
				The number of events to read. If end and length are provided
				length will be ignored.
			
			get_all = False
				If True, return all events. Overrides start, end and length.
			
			pulse_pointers_only = False
				If set to True, instead of loading PODs into 'pulse_data'
				we will load the file pointers to the PODs for each channel
				in each event.
			
			max_channel = None
				Read only up to this channel if set. Indexed from 0 so PMT
				122 would be max_channel of 121
			
			recalculate_baseline = None
				The baseline obtained by the DAQ and saved in the EVT file is
				biased by (very) roughly 0.5ADU. If this option is set to True,
				the baseline will be recomputed using the first
				recalculate_baseline_samples of each POD waveform.
				If this option is not specified (that is, set to None) the
				default defined when the class was instantiated will be used.
		
			recalculate_baseline_samples = None
				The number of samples that will be used to recalcuate the
				baseline if recalculate_baseline is set to True.
				If this option is not specified (that is, set to None) the
				default defined when the class was instantiated will be used.
		_________________________
		Returns:
			A dictionary with keys:
				'luxstamp_samples', 'active_channels', 'pulse_nums', 
				'pulse_starts', 'pulse_lengths', 'pulse_baselines',
				'pulse_data', 'pulse_time_arr'
			Each dictionary entry is a list of the specified items
		
		"""
		#______________________________________________
		if get_all:
			# get all events
			start = 0
			end = self.evt_number
		else:
			if start == None:
				# User did not provide a start event index. Read the current
				# event.
				start = self.current_evt_ind
			if end == None:
				# User did not provide an end event index. Read the current
				# event only.
				end = start+length
		#______________________________________________
		# check how the baseline will be treated
		if recalculate_baseline == None:
			# if the option is not provided, use the default
			recalculate_baseline = self.recalculate_baseline
		if recalculate_baseline_samples == None:
			# if the option is not provided, use the default
			recalculate_baseline_samples = self.recalculate_baseline_samples
		#______________________________________________
		# check if the user doesn't want too many events
		if end > self.evt_number:
			raise(InsufficientEventNumnerError(end, self.evt_number))
		if self.debug: print 'get_event_data from %d to %d' %(start,end)
		# check if there is a limit on the channels. If non specified by
		# the user, use the default object one.
		if max_channel == None: max_channel = self.max_channel
		#______________________________________________
		# go through the events requested by the user
		ordered_event_return_keys = \
			['luxstamp_samples', 'active_channels', 'pulse_nums', \
			'pulse_starts', 'pulse_lengths', 'pulse_baselines',
			'pulse_data', 'pulse_time_arr']
		# The above constiture the keys to the return dictionary initialized
		# below. The MUST be in the same order as the returned values of
		# read_evt_file_data().
		retdict = {}
		for skey in ordered_event_return_keys: retdict[skey] = []
		for i in range(start, end):
			# self.open() will take care of closing the previous file if
			# neccessary and open the new one unless it's already opened.
			self.open(self.evt_file_inds[i])
			if self.current_evt_ind != i:
				# Go to the file location if we're not already there
				if self.debug: print 'moving file pointer independently'
				self.filehandle.seek(self.evt_gids[i])
				# set the new location
				self.current_evt_ind = i
			elif self.fileread_in_progress:
				# If this was triggered there was a problem with a previous
				# read. We must therefore go to a correct file locaton.
				if self.debug:
					print 'moving file pointer after a bad previous read'
				self.filehandle.seek(self.evt_gids[i])
				self.fileread_in_progress = False
			# else:
				# Assume that we are at the correct file location if the
				# current index matches hence don't go there. Moving around a
				# file involves accessing the harddrive uneccessarily so we want
				# to avoid it if possible. But we can perform error checks on
				# the current file location.
			#_______________________
			# Check the file location. Raise an exception if it's not correct.
			floc = self.filehandle.tell()
			if floc != self.evt_gids[i]:
				if self.debug: print 'Error on the %dth event read' %(i-start+1)
				raise(EVTFilePositionError(floc, self.evt_gids[i]))
			#_______________________
			# Finally, read the event
			self.fileread_in_progress = True
			evtret = read_evt_file_data(self.filehandle,
				pulse_pointers_only=pulse_pointers_only, max_channel=max_channel,
				ns_per_sample = self.ns_per_sample, 
				read_trigger_info=self.read_trigger_info,
				recalculate_baseline = recalculate_baseline,
				recalculate_baseline_samples = recalculate_baseline_samples,
				preamp_gain = self.preamp_gain, postamp_gain = self.postamp_gain)
			self.fileread_in_progress = False
			# put the event data into a return dictionary
			for j in range(len(ordered_event_return_keys)):
				skey = ordered_event_return_keys[j]
				retdict[skey].append(evtret[j])
			#_______________________
			# reset the current event index
			self.current_evt_ind = i+1
		if self.debug: print 'Read %d events' %(i-start+1)
		return retdict
	
	def open(self, fileind, go_to_start=False):
		"""
		Open the file of index fileind
	
		_________________________
		Options: 
			fileind
				The index of the file to open
		
		Optional:
			go_to_start=False
				For a file that's already open this will go back
				to start.
	
		_________________________
		Returns: None
		"""
		# If there is a file currently opened, close it
		# unless it's the same file. 
		if self.currently_opened_file_index != fileind:
			if self.filehandle != None:
				# File is opened but it's not the right one.
				self.close()	# close it
			# open the file
			self.filehandle = open(self.filelist[fileind], 'rb')
			self.currently_opened_file_index = fileind
			# go to the first event
			self.current_evt_ind = \
				self.evt_number_per_file[:self.currently_opened_file_index].sum()
			self.filehandle.seek(self.evt_gids[self.current_evt_ind])
			if self.debug: print 'Opened file %d, %s' \
				%(self.currently_opened_file_index, \
				self.filelist[self.currently_opened_file_index])
		elif go_to_start:
			# the correct file is opened. However, we could be at
			# the wrong position within it. If an open was requested
			# we should assume the user wants to start at location 0.
			self.current_evt_ind = \
				self.evt_number_per_file[:self.currently_opened_file_index].sum()
			self.filehandle.seek(self.evt_gids[self.current_evt_ind])
			if self.debug: print 'File moved to start pointer'
	
	def close(self):
		"""
		Force the currently opened file closed.
	
		_________________________
		Options: None
	
		_________________________
		Returns: None
		"""
		if self.filehandle != None:
			# file is actually open
			self.filehandle.close()			# close the file
			if self.debug: print 'File %d has been closed' \
				%(self.currently_opened_file_index)
		self.filehandle = None			# set the flags
		self.currently_opened_file_index = None
	
	def __del__(self):
		if self.debug: print 'Object destructor initiated'
		self.close()	# close the file


def prepare_evt_info(file_or_directory, filematch = '*.evt'):
	"""
	We will open all of the files and get the headers from them. 
	Then we'll perform some checks to make sure that the headers
	agree and save the neccessary info.
	
	_________________________
	Options:
		file_or_directory
			An evt file, list of evt files or a directory containing
			evt files.
	
	Optional:
		filematch = '*.evt'
			The format of files we're going to use. Can use unix wildcards.
	
	_________________________
	Returns:
		filenum, filelist, xmlsettings, evt_number_per_file, \
		date_and_time, live_time_header, evt_ids, evt_gids, evt_file_inds
	"""
	#_________________________________________________
	# get a list of files to process
	if isinstance(file_or_directory, n.ndarray) or \
		isinstance(file_or_directory, list):
		# a list of filenames; just what we need
		filelist = file_or_directory
	elif os.path.isdir(file_or_directory):
		# directory containing evt files given
		filelist = glob.glob(os.path.join(file_or_directory, filematch))
	else:	# assume it's a string
		filelist = [file_or_directory]
	filenum = len(filelist)
	#_________________________________________________
	# initialize data storage lists
	evt_number_per_file = []		# number of events per each file
	evt_ids = []					# IDs of the events			
	evt_gids = []					# the file pointers of each event
	evt_file_inds = []				# the idecies of the file. Matched to evt_gids
	date_and_time = []				# date and time. I don't know what the format is
	#loc_exp_run_version = []		# don't bother
	live_time_header = {}			# livetime info
	for i in range(filenum):		# go through the files
		filename = filelist[i]		# current file
		f = open(filename, 'rb')	# Open the file
		#_______________________
		# read the file header
		tevt_number_per_file, tevt_ids, tevt_gids, txmlsettings, \
			tdate_and_time, tloc_exp_run_version, tlive_time_header = \
				read_evt_file_header(f)
		#_______________________
		if i == 0:
			# the settings should be the same for each file so store only once
			xmlsettings = txmlsettings
		#else:	# check if the settings are the same
			# 	don't check
		#_______________________
		# Perform an integrity check. The current filehandle location should
		# equal to the first evt_gids
		floc = f.tell()
		if floc != tevt_gids[0]:
			raise(EVTFilePositionError(floc, evt_gids[0]))
		#_______________________
		# save the neccessary data
		#
		# The number of events per file the date and time and the
		# livetime header are properties that belong to the file
		# hence we get one value per file.
		evt_number_per_file.append(tevt_number_per_file)
		date_and_time.append(tdate_and_time)
		for skey in tlive_time_header.keys():
			# store dictionary of lists, not vice versa
			if i == 0:	# initialize the list first time around
				live_time_header[skey] = [tlive_time_header[skey]]
			else:		# append to an existing list
				live_time_header[skey].append(tlive_time_header[skey])
		#
		# The event IDs, event gids and event file indeces are properties
		# of each event hence we get one value per each event.
		evt_ids.extend(tevt_ids)
		evt_gids.extend(tevt_gids)
		# The above evt_gids tell us where the events are in a file. But we
		# must also know which file each evt_gid is in. Hence evt_file_inds
		# below.
		evt_file_inds.extend(n.zeros(tevt_gids.size,dtype=n.int64)+i)
		#_______________________
		# close the file
		f.close()
	# convert some of the lists to numpy arrays
	evt_number_per_file = n.array(evt_number_per_file, dtype=n.uint32)
	date_and_time = n.array(date_and_time, dtype=n.uint32)
	evt_ids = n.array(evt_ids, dtype=n.uint32)
	evt_gids = n.array(evt_gids, dtype=n.uint32)
	evt_file_inds = n.array(evt_file_inds, dtype=n.int64)
	# Should we sort by evt_ids in case the user did not enter files
	# in order. In principle, since the events are read in one by one 
	# we want the events from a single file grouped together so that we
	# don't have to keep opening and closing files which means we want
	# them in order. But what if the user entered them in order and it's
	# the internal IDs that are out of order?
	# Nope. If we want to quickly read in order then the same order must
	# be maintained as during the original read. We cannot be jumping
	# around the file looking for different evt_gids.
	return filenum, filelist, xmlsettings, evt_number_per_file, \
		date_and_time, live_time_header, evt_ids, evt_gids, evt_file_inds


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
def read_evt_file_header(f):
	"""
	Read the file header. This is everything preceeding the event itself.
	
	_________________________
	Options:
		f
			An opened file handle
	
	_________________________
	Returns:
		evt_number_per_file, evt_ids, evt_gids, xmlsettings,
		date_and_time, loc_exp_run_version, live_time_header
	"""
	#__________________________ Initial Info
	# read the first 4 bytes. Is supposed to be 0x01020304 if correct endian
	checknum = n.fromfile(f,dtype=dt.uint32,count=1)[0]
	dt.set_endian(checknum)	# set endian if neccessarry
	#print '%x' %(checknum)
	# the next 4 bytes is the lengths of the settings strings
	settings_string_length = n.fromfile(f,dtype=dt.uint32,count=1)[0]
	#__________________________ read in the event builder and DAQ settings
	event_builder_daq_settings = \
		n.fromfile(f,dtype=dt.int8,count=settings_string_length)
	# build the settings string
	settings_xml = ''
	for i in range(event_builder_daq_settings.size):
		settings_xml += '%c' %(event_builder_daq_settings[i])
	## clean up the xml tag which is causing a lot of trouple.
	settings_xml = re.sub('<.+?xml.+?>\n?', '', settings_xml)
	# clean up waste Null characters
	settings_xml = re.sub('\x00', '', settings_xml)
	#print settings_xml
	# Grab the xml settings
	xmlsettings = xml2dict.xml2dict('<settings>\n'+settings_xml+'\n</settings>')['settings']
	#__________________________ File Header
	# the next 4 bytes tell us the header info
	checknum = n.fromfile(f,dtype=dt.uint32,count=1)[0]
	dt.set_endian(checknum)	# set endian if neccessarry
	date_and_time, loc_exp_run_version, evt_number_per_file = \
		n.fromfile(f,dtype=dt.uint32,count=3)
	#__________________________ 
	# now read in the event numbers (IDs) and byte location of event GIDs
	temparr = n.fromfile(f,dtype=dt.uint32,count=evt_number_per_file*2)
	evt_ids = temparr[0::2]
	evt_gids = temparr[1::2]		# file pointers to the events
	del temparr
	#__________________________ Live Time Header
	
	seqnum = n.fromfile(f,dtype=dt.uint16,count=1)[0]
	# read in timestamp latches and ends
	temparr = n.fromfile(f,dtype=dt.uint64,count=seqnum*2)
	live_time_header = {'seqnum':seqnum, 'timestamp_latch':temparr[0::2], \
		'timestamp_end':temparr[1::2]}
	del temparr, seqnum
	# return the file information
	return evt_number_per_file, evt_ids, evt_gids, xmlsettings, \
		date_and_time, loc_exp_run_version, live_time_header


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
def read_evt_file_data(f, pulse_pointers_only=False, max_channel=None,
	ns_per_sample = 10., read_trigger_info = True, 
	recalculate_baseline = True, recalculate_baseline_samples = 22,
	preamp_gain = 5., postamp_gain = 1.5):
	"""
	Read a single event data from an evt file.
	
	_________________________
	Options:
		f
			An opened file handle
	
	Optional:
		pulse_pointers_only = False
			If set to True, instead of loading PODs into 'pulse_data'
			we will load the file pointers to the PODs for each channel
			in each event.
		
		max_channel = None
			Read only up to this channel if set. Indexed from 0 so PMT
			122 would be max_channel of 121
		
		ns_per_sample = 10.
			Assume that each sample lasts 10ns
		
		read_trigger_info = True
			
		
		recalculate_baseline = True
			The baseline obtained by the DAQ and saved in the EVT file is biased
			by (very) roughly 0.5ADU. If this option is set to True, the
			baseline will be recomputed using the first
			recalculate_baseline_samples of each POD waveform.
		
		recalculate_baseline_samples = 22
			The number of samples that will be used to recalcuate the baseline
			if recalculate_baseline is set to True.
		
		preamp_gain = 5.
			The pre-amplifier gain
		
		postamp_gain = 1.5
			The post-amplifier gain
	_________________________
	Returns:
		luxstamp_samples, active_channels, voltres_voltoffset, \
		pulse_nums, pulse_starts_lengths_baselines, pulse_data
		
		luxstamp_samples is the trigger time stamp, a single value.
		active_channels is a numpy array of channels (PMTs) participating
			in this event. Indexed from 0
		voltres_voltoffset is a 2 by <channe number> array that contains
			the voltage resolution and voltage offset for each
			channel.
		pulse_nums is an array containing the number of pulses for each channel
		pulse_starts_lengths_baselines is a list, one entry per channel. Each
			entry is a 3 by <pulsenum> array containing the pulse start, length
			and basline for each POD.
		pulse_data is a list of lists (for each channel and pulse) containing
			the POD waveform arrays.
	
	"""
	#
	##__________________________ Event GID Header
	#gidheader = n.fromfile(f,dtype=dt.uint32,count=3)
	##[0] Date and time of run (filename prefix)
	##[1] Location, experiment, run version
	##[2] Event number
	#
	# go to event
	f.seek(12, 1)	
	# Jump to number of channels per event and then skip the event size entry
	channels_per_event = n.fromfile(f,dtype=dt.uint32,count=2)[0]
	##[4] Event size in bytes
	#__________________________
	if read_trigger_info:
		# read trigger information
		#ddc_trig_timestamp = n.fromfile(f,dtype=dt.uint64,count=1)[0]
		#ddc_trig_seqnum, max_filt_response = n.fromfile(f,dtype=dt.uint32,count=2)
		#max_channel_id = n.fromfile(f,dtype=dt.int8,count=1)[0]
		f.seek(17, 1)
		# we need the ddc board number to skip towards later
		ddc_board_num = n.fromfile(f,dtype=dt.uint16,count=1)[0]
		#
		#ddc_board_s1_vector = n.fromfile(f,dtype=dt.int8,count=ddc_board_num)
		#ddc_board_s2_vector = n.fromfile(f,dtype=dt.int8,count=ddc_board_num)
		#ddc_check_byte = n.fromfile(f,dtype=dt.int8,count=1)[0]
		f.seek(2*ddc_board_num+1, 1)
	## records. Not sure what these are
	#record_format, record_size = n.fromfile(f,dtype=dt.uint32,count=2)
	f.seek(8, 1)
	#__________________________
	# Grab the event timestamp
	luxstamp_samples = n.fromfile(f,dtype=dt.uint64,count=1)[0]
	# initialize lists that will store various information
	active_channels = []	# will keep a list of channels for which we have PODs
	pulse_nums = []			# will keep the number of PODs for each active channel
	pulse_starts = []
	pulse_lengths = []
	pulse_baselines = []
	# Will keep arrays of POD starts, lengths and baselines for each channel
	#
	pulse_data = []			# will keep lists of arrays containing POD data
	pulse_time_arr = []		# will store the pulse times
	#voltres = []
	#voltoffset = []
	# Will keep arrays of voltage resolutions and voltage offsets
	# for each channel
	#
	# Step through a list of active channels. Some of them may be empty
	for jj in range(channels_per_event):
		# not sure how this is used
		binary_datatype = n.fromfile(f,dtype=dt.int8,count=1)[0]
		tvoltres, tvoltoffset = n.fromfile(f,dtype=dt.double,count=2)
		# skip the time resolution since it is wrong
		f.seek(8, 1)
		temp_read_uint32 = n.fromfile(f,dtype=dt.uint32,count=5)
		if temp_read_uint32[-1] == 0:
			# Pulse number is 0. Skip this channel
			continue
		pulsenum = temp_read_uint32[4]
		# skip pretrigger, evt_size, pulse_pretrigger, pulse_posttrigger
		#
		# - Number of samples relative to trigger that pulse begins,
		# - Number of samples in the pulse (includes pulse detect
		# 	pretrigger and pulse end posttrigger)
		# - Pulse baseline in ADC counts
		tpulse_starts, tpulse_lengths, tpulse_baselines = \
			n.fromfile(f,dtype=dt.int32,count=pulsenum*3).reshape((3,pulsenum))
		tpulse_baselines = n.float32(tpulse_baselines)
		#__________________________
		# save the pulse data
		if pulse_pointers_only:	# don't read the POD waveforms
			# save only the file pointer to the PODs of that channel
			tpulse_data = f.tell()
			# advance the file
			f.seek(tpulse_lengths.sum()*2)
		else:
			#__________________________
			# number of bytes to read
			all_pods_length = tpulse_lengths.sum()
			# This is the "new" way of doing in. All PODs in a single channel are
			# read in as a continous 1D array and time array is created for each.
			# The data is also converted to mV
			#__________________________
			# read the PODs
			tpulse_data = n.float32(n.fromfile(f,dtype=dt.uint16,count=all_pods_length))
			#__________________________
			if tpulse_data.size != all_pods_length:
				# check if the POD rate agrees with what we were
				# supposed to read
				raise(PODLengthError(tpulse_data.size, all_pods_length))
			#__________________________
			if recalculate_baseline:
				# correct the POD digitization bias
				currentpos = 0
				for kk in range(pulsenum):
					tpulse_baselines[kk] = \
						tpulse_data[currentpos:currentpos + \
						tpulse_lengths[kk]][: \
						recalculate_baseline_samples].mean()
					# increment position counter
					currentpos += tpulse_lengths[kk]
			#__________________________
			# convert to mV*ns AND create a time array
			#
			# get the baseline to be in the same units
			tpulse_baselines = n.float32(tpulse_baselines) * tvoltres*1000
			# pre-define the time array
			tpulse_time_arr = n.arange(tpulse_data.size, dtype=n.int64)
			currentpos = 0
			for kk in range(pulsenum):
				tpulse_data[currentpos:currentpos+tpulse_lengths[kk]] = \
					tpulse_baselines[kk] - \
					tpulse_data[currentpos:currentpos + \
					tpulse_lengths[kk]] * tvoltres*1000 + tvoltoffset*1000
				# create the time array by setting the start of each pulse to 0
				# and then addting pulse_start
				tpulse_time_arr[currentpos:currentpos+tpulse_lengths[kk]] += \
					-tpulse_time_arr[currentpos] + tpulse_starts[kk]
				# increment position counter
				currentpos += tpulse_lengths[kk]
			# Apply the postamp and preamp gains
			tpulse_data *= ns_per_sample / preamp_gain / postamp_gain
		# end once we reach the maximum channel
		if max_channel != None and jj > max_channel: continue
		# otherwise save the data
		# 
		active_channels.append(jj)		# store the current channel
		# pulse number
		pulse_nums.append(pulsenum)
		# voltage resolution in this channel, voltage offset in this channel,
		# time resolution in this channel
		#voltres.append(tvoltres)
		#voltoffset.append(tvoltoffset)
		pulse_starts.append(tpulse_starts)
		pulse_lengths.append(tpulse_lengths)
		pulse_baselines.append(tpulse_baselines)
		# pulse data
		pulse_data.append(tpulse_data)
		# pulse time array
		pulse_time_arr.append(tpulse_time_arr)
	# compile the other arrays into the dictionary
	active_channels = n.array(active_channels)
	pulse_nums = n.array(pulse_nums)
	#voltres = n.array(voltres)
	#voltoffset = n.array(voltoffset)
	# return the event information
	#return luxstamp_samples, active_channels, voltres, voltoffset, \
	#	pulse_nums, pulse_starts, pulse_lengths, pulse_baselines, pulse_data
	return luxstamp_samples, active_channels, pulse_nums, \
		pulse_starts, pulse_lengths, pulse_baselines, pulse_data, pulse_time_arr
	# the pulse_data is a 1D list of numpy arrays.


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
#		Exception definitions
#

class InsufficientEventNumnerError(Exception):
	"""
	Raised when the user tried to get events past the number
	actually available.
	"""
	def __str__(self):
		messgstr = "\nRequested %d final event\n" \
			%(self.args[0])
		messgstr += "but only %d are available\n" \
			%(self.args[1])
		return messgstr

class EVTFilePositionError(Exception):
	"""
	Raised if the current location of the file pointer (presumably
	at the start of	an event) does not match the evt_gid. This 
	should not happen unles the file was read incorrectly (i.e. the
	wrong number of bytes was read) or written incorrectly.
	"""
	def __str__(self):
		messgstr = "\nThe current file pointer location %d does not match\n" \
			%(self.args[0])
		messgstr += "the expected location %d as read from file\n" \
			%(self.args[1])
		return messgstr


class PODLengthError(Exception):
	"""
	Raised if the POD length read from the file does not match the
	length of the POD data. This should not happen since the number
	we use the POD length to know how many bytes to read into the
	data POD. The only way this would happen if the file is in
	some way truncated or if the last POD length extends beyond the
	end of the file.
	"""
	def __str__(self):
		messgstr = "\nThe POD data length %d does not match\n" %(self.args[0])
		messgstr += "the expected length %d as read from file\n" \
			%(self.args[1])
		return messgstr


class SettingMissingError(Exception):
	"""
	Raised if an XML tag/setting is missing.
	"""
	def __str__(self):
		messgstr = "\nSetting %s was not found within the XML\n" %(self.args[0])
		messgstr += "information included in the EVT file.\n"
		# add in a string messege
		messgstr += "%s\n" %(self.args[1])
		return messgstr


