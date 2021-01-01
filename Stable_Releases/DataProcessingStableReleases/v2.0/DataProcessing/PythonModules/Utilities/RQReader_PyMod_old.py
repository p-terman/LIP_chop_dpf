# RQReader_PyMod
#
# You can load an entire .rq file with ReadBinaryBlock:
#  rq,settings,livetime = RQReader_PyMod.ReadBinaryBlock('lux10_20130505T0712_f000001_cp05457.rq.hdf5')
#
# Optional boolean argument 'verbosity_flag' for printing each RQ loaded with dimensions
#
# 20140625 CHF - Created
# 20140703 CHF - Added int32, uint8 to formats list
#
from struct import unpack, calcsize
from xml2dict import xml2dict
import numpy as np
from collections import OrderedDict as odict	# ordered keys in a dictionary

# Define mapping from binary format name to struct.unpack codes
fmts = {'char':'c', 'short':'h', 'ushort':'H', 'int':'i', 'uint':'I',\
        'long':'l', 'ulong':'L', 'longlong':'q', 'ulonglong':'Q', 'double':'d',\
        'bool':'b', 'int32':'i', 'single':'f',\
        'int32':'i','uint32':'I','uint64':'Q','uint8':'B','int8':'b'}

def ReadRQFile(filename,verbosity_flag=0):
	"""
    This function reads an entire .rq file. 

    It outputs three dictionaries: rq, settings, livetime
    _________________________
    Inputs:
        filename
            Filename (with full path if not in same folder as this function)
        verbosity_flag
        	Boolean flag for displaying text along the way (default:0)
    Returns:
		rq_dict
			Dictionary with rqs as keys and np.array for values
		settings_dict
			XML settings turned into a dictionary
		livetime_dict
			Dictionary with livetime_latch_samples and livetime_end_samples keys and np.array for values.
    Versioning:
        20140625 CHF - Created

    To do:
    	Handling .rq.gz files
    	Error catching
    _________________________
    """
	# Read file
	f = open(filename,'rb')

	##########################
	###### BLOCK 0: XML ######
	##########################
	# Start off by reading xml header
	f.seek(0)
	a = unpack('<i',f.read(calcsize('i')))
	size_of_xml = unpack('<i',f.read(calcsize('i')))
	xml_string = f.read(size_of_xml[0])

	# Format xml header (it's broken)
	last_tag = '</daq_settings>'
	del_ind = xml_string.rfind(last_tag)
	xml_string_sanitize = xml_string[0:del_ind+len(last_tag)]
	xml_string_sanitize = '<settings>' + xml_string_sanitize.replace('\n','') + '</settings>'
	settings_dict = xml2dict(xml_string_sanitize)

	##################################
	###### BLOCK 1: Description ######
	##################################
	# Read header block
	var_names,var_types,var_len_str = ReadHeaderBlock(f)
	var_len = [int(x) for x in var_len_str]

	# This seems useless, but has to be read to advance
	other_stuff = unpack('<i',f.read(calcsize('i')))

	# Make dictionaries to give ReadDataBlock as input
	var_types_dict = dict(zip(var_names, var_types))
	var_len_dict = dict(zip(var_names, var_len))

	# Loop per variable names and read binary block for each
	header = dict()
	for var in var_names:
		header[var] = ReadDataBlock(f,var_types_dict[var],var_len_dict[var])

	# If dataset_name is a key, and it starts with a dot
	if header.get('dataset_name'):
		if header['dataset_name'][0] is '.':
			# then remove that pesky '.' before filename_prefix
			header['dataset_name'] = header['dataset_name'][1:] 

	number_of_events_in_file = int(header['nb_evts_in_file'][0])

	##########################
	###### BLOCK 2: RQs ######
	##########################
	# Read header block
	rq_var_names,rq_var_types,rq_var_len_str = ReadHeaderBlock(f)

	# Convert RQ dimension list from comma-separated str to int
	rq_var_len = list()
	for rq in rq_var_len_str:
		var_temp = [int(x) for x in rq.split(',')]
		rq_var_len.append(var_temp)

	# Make dictionaries to give ReadDataBlock as input
	rq_var_types_dict = dict(zip(rq_var_names, rq_var_types))
	rq_var_len_dict = dict(zip(rq_var_names, rq_var_len))

	# This seems useless, but has to be read to advance
	other_stuff = unpack('<i',f.read(calcsize('i')))

	# Find the number of pulses we are keeping
	# Using pulse_area_phe as template since it's a very basic rq
	# that should be there
	num_pulses = int(rq_var_len_dict['pulse_area_phe'][0])

	# Initialize output dictinary with numpy arrays
	rq_dict = dict()
	for rq in rq_var_names:
		# If it has dimensions of [1], it's one per event.
		# Also, if it has dimensions of total event numbers, it's not iterable, but
		# the size should also be number_of_events_in_file
		if (rq_var_len_dict[rq][0] == 1) or (rq_var_len_dict[rq][0] == number_of_events_in_file):
			total_dimensions = number_of_events_in_file
		# If it has dimensions of num pulses, then concatenate event number dimension
		elif rq_var_len_dict[rq][0] == num_pulses:
			total_dimensions = [number_of_events_in_file]+rq_var_len_dict[rq]
		else:
			print '***ERROR: Funny dimensions found for %s' % rq

		# Initialize with nans (they will all get populated, don't worry)
		rq_dict[rq] = np.empty(total_dimensions)
		rq_dict[rq].fill(np.nan)

	# Loop per event number and read each RQ at a time (not optimal, but does the job)
	for ee in xrange(number_of_events_in_file):
		
		if verbosity_flag:
			print "*** Loading event %d ***" % ee
		else:
			print '.',

		# This ensures they are read in order of appearance in file!
		for rq in rq_var_names:

			# Read the rq value
			output = ReadDataBlock(f,rq_var_types_dict[rq],rq_var_len_dict[rq])

			# Reshape the output if it's supposed to have multiple dimensions
			num_dims = len(rq_var_len_dict[rq])

			# If it's one-per-event
			if (num_dims == 1) and (rq_var_len_dict[rq][0] == 1):
				if verbosity_flag:
					print '%s: dimension %d' % (rq,1)
				rq_dict[rq][ee] = np.array(output)
			# If it's one-per-pulse
			elif (num_dims == 1) and (rq_var_len_dict[rq][0] == num_pulses):
				if verbosity_flag:
					print '%s: dimension %d' % (rq,num_pulses)
				rq_dict[rq][ee,:] = np.array(output)
			# If it's per-PMT (and, of course, per-pulse)
			elif num_dims == 2:
				if verbosity_flag:
					print '%s: dimensions %dx%d' % (rq,rq_var_len_dict[rq][0],rq_var_len_dict[rq][1])
				nparr = np.array(output)
				nparr_reshaped = nparr.reshape(rq_var_len_dict[rq])
				rq_dict[rq][ee,:,:] = nparr_reshaped


	###############################
	###### BLOCK 3: Livetime ######
	###############################

	# Read header
	livetime_var_names,livetime_var_types,livetime_var_len_str = ReadHeaderBlock(f)
	livetime_var_len = [int(x) for x in livetime_var_len_str]

	livetime_var_types_dict = dict(zip(livetime_var_names, livetime_var_types))
	livetime_var_len_dict = dict(zip(livetime_var_names, livetime_var_len))

	# Read binary data
	livetime = dict()
	for var in livetime_var_names:
		livetime[var] = np.array(ReadDataBlock(f,livetime_var_types_dict[var],livetime_var_len_dict[var]))

	print '\nFinished loading %s' % filename

	# Close file
	f.close()

	# Return data
	return rq_dict, settings_dict, livetime


def ReadHeaderBlock(f):
	"""
    This function reads a header block in a LUX .rq file. This should only be called
    when the file read location is the header string length.

    It does the following:
        1) Reads the size of the header string length
	    2) Reads the header string
	    	(e.g.'livetime_latch_samples;uint64;20;livetime_end_samples;uint64;20;') 
	    3) It parses the values into var_names, var_types and var_len_str lists
    _________________________
    Inputs:
        f
            Opened .rq file object [e.g. f = open(filename.rq,'rb')]
    Returns:
        var_names
            List of variable names
        var_types
            List of variable types, order corresponding to var_names
        var_len_str
        	List of variable lengths (strings!), order corresponding to var_names
    Versioning:
        20140625 CHF - Created
 	
    To do:
    	Add option to read this by providing a file location value (fseek)
    	Make output a list of ints instead of strings
    	Error catching
    _________________________
    """
	# Read header size
	header_size = unpack('<H',f.read(calcsize('H')))
	header_string = f.read(header_size[0])

	# Remove final ';' to avoid additional empty list entry
	if header_string[-1]==';':
		header_string = header_string[:-1]

	# Asign outputs
	var_names = header_string.split(';')[0::3]
	var_types = header_string.split(';')[1::3]
	var_len_str = header_string.split(';')[2::3]

	return var_names,var_types,var_len_str

def ReadDataBlock(f,var_types,var_lengths):
	"""
	This function reads a binary data block in a LUX .rq file. This should only be called
	after the header block has been read.

	It outputs a 1-D list of values corresponding to the data read.
	For multi-dimensional data e.g. (10,122), it reads 10x122=1220 and still outputs a 1-D list.
	Convert to np.array and perform a reshape to get it to the right shape.
	_________________________
	Inputs:
	    f
	        Opened .rq file object [e.g. f = open(filename.rq,'rb')]
	    var_types
	    	Variable type to read (e.g. 'char', 'uint32', etc/)
	    var_lengths
	    	Variable length (e.g. 10, or [10,122], etc.)
	    	The variable lenghts should be integers! (ReadHeaderBlock currently outputs string)
	Returns:
		output
			1-D list with the data in the specified format
			Convert to np.array and reshape for multi-dimensional data
	Versioning:
	    20140625 CHF - Created

	To do:
		Add option to read this by providing a file location value (fseek)
		Error catching
	_________________________
	"""
	# Total is the product of all dimensions
	total_var_length = np.prod(var_lengths)

	# Establish format from lookup dictionary and unpack
	fmt = fmts[var_types]*total_var_length
	temp_value = unpack(fmt,f.read(calcsize(fmt)))

	if var_types == 'char':
		# The if statement prevents empty line characters '\x00' from being used as a string
		output = ''.join(elem[0] if elem is not '\x00' else '' for elem in temp_value )
	else:
		output = temp_value
	return output

