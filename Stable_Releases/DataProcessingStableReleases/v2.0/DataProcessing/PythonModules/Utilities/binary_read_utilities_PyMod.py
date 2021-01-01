# binary_read_utilities_PyMod
# Tomasz Biesiadzinski
# 07/22/2014

#
# 20140722 TPB - Created

import re
import sys

#-----------------------------------------------------------------------------
class dtypes:
	def __init__(self, endian = None):
		"""
		A convenient object for handling binary reads/writes. It sets the 
		numpy datatypes to the proper endian. It also returns the numpy 
		datatype for an RQ variable type and the RQ variable type for a numpy
		datatype.
	
		_________________________
		Optional:
			endian = None
				Allows the user to set the initial endian. Otherwise the system
				value will be used.
	
		_________________________
		Returns:
			An instance of the dtypes object.
	
		"""
		#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
		# Here we define the data types with default endianness
		self.system_endian = sys.byteorder	# default endian of the system
		#
		if endian == None:
			# current endian of the object
			self.current_endian = self.system_endian
		elif endian == 'little' or endian == 'l' or endian == 'L' or \
			endian == '<':
			self.current_endian = 'little'
		elif endian == 'big' or endian == 'b' or endian == 'B' or \
			endian == '>':
			self.current_endian = 'big'
		else:
			raise(UnknownEndianError(endian))
			
		#
		# set the default data parameters
		if self.current_endian == 'little':
			# little endian
			self.uint8 = '<u1'
			self.uint16 = '<u2'
			self.uint32 = '<u4'
			self.uint64 = '<u8'
			self.int8 = '<i1'
			self.int16 = '<i2'
			self.int32 = '<i4'
			self.int64 = '<i8'
			self.single = '<f4'
			self.double = '<f8'
		else: # big endian
			self.uint8 = '>u1'
			self.uint16 = '>u2'
			self.uint32 = '>u4'
			self.uint64 = '>u8'
			self.int8 = '>i1'
			self.int16 = '>i2'
			self.int32 = '>i4'
			self.int64 = '>i8'
			self.single = '>f4'
			self.double = '>f8'
		# put the formats into a list
		self.formatdefs = [self.uint8, self.uint16, self.uint32, self.uint64, \
			self.int8, self.int16, self.int32, self.int64, self.single, \
			self.double]
	#
	def set_endian(self, checknum, allow_endian_override=False):
		"""
		Flip the endian depending on the value of a check number read
		from the file.
	
		_________________________
		Options:
			checknum
				A hex number whos value is supposed to be 0x01020304
				when read correctly. 
	
		Optional:
			allow_endian_override=False
				I noticed in some of the dat files that checknum is just
				full of zeros. This keeps the endiannes the same if
				checknum is zero.
	
		_________________________
		Returns:
			Nothing. Datatype defintions used by numpy.fromfile() are
			altered globally without being returned. 
		"""
		if checknum == 0x01020304:
			# the current endian appears to be correct
			pass
			# good. do nothing
		elif checknum == 0x04030201:
			# flip endian
			self.swap_endian()
		elif checknum == 0x00000000 and allow_endian_override:
			# I found that in the first .dat file I was using the endian
			# check wasn't there at the beggining of the file and
			# instead the first value was 0. Hence this option was added.
			pass
		else:
			# something went terribly wrong
			raise(BadEndianIndicatorError(checknum))
	#
	def swap_endian(self):
		for i in range(len(self.formatdefs)):
			if self.formatdefs[i][0] == '<':
				self.formatdefs[i] = '>' + self.formatdefs[i][1:]
				# update current endian to big
				self.current_endian = 'big'
			else:
				self.formatdefs[i] = '<' + self.formatdefs[i][1:]
				# update current endian to little
				self.current_endian = 'little'
		# Reset the data types definitions
		self.uint8, self.uint16, self.uint32, self.uint64, \
			self.int8, self.int16, self.int32, self.int64, self.single, \
			self.double = self.formatdefs
	#
	def set_ndarray_endian(self, array):
		"""
		Reverse the byte order of a numpy array based on the endian value.
		
		Native numpy datatypes are always defined as native ("=") which is
		hopefully the same as the system endian.
		
		_________________________
		Options:
			array
				The numpy array to byteswap.
		
		_________________________
		Returns:
			 byteswapped array
		"""
		if self.current_endian != self.system_endian:
			# swap the array endian if the current object endian does not agree
			# with the system (native, hopefully) endian.
			return array.byteswap()
	#
	def get_datatype(self, variable_type):
		"""
		Convert a string descriptor into a datatype understood by python.
	
		Note that normally making a dictionary with string descriptors as keys
		and datatypes as values would be the right way to accomplish this. But,
		storing the datatype in a dictionary appears to make a copy of the
		variable which means that set_endian() would not be able to change
		the endianness of that datatype. So we're using a function instead.
		
		_________________________
		Options:
			variable_type
				Type of a variable written in the RQ binary file.
		
		_________________________
		Returns:
			 Two values: The numpy datatype and the number of bytes per datatype
		"""
		if variable_type == 'char':
			return 'c', 1	# 1 byte, 8 bits per character
		elif variable_type == 'short' or variable_type == 'int16':
			return self.int16, 2		# 2 bytes
		elif variable_type == 'int' or variable_type == 'long' or \
			variable_type == 'int32':
			return self.int32, 4
		elif variable_type == 'longlong' or variable_type == 'int64':
			return self.int64, 8
		elif variable_type == 'ushort' or variable_type == 'uint16':
			return self.uint16, 2
		elif variable_type == 'uint' or variable_type == 'ulong' or \
			variable_type == 'uint32':
			return self.uint32, 4
		elif variable_type == 'ulonglong' or variable_type == 'uint64':
			return self.uint64, 8
		elif variable_type == 'single' or variable_type == 'float32':
			return self.single, 4
		elif variable_type == 'double' or variable_type == 'float64':
			return self.double, 8
		elif variable_type == 'bool' or variable_type == 'uint8':
			return self.uint8, 1
		elif variable_type == 'int8':
			return self.int8, 1
		else:
			raise(UnknownDataTypeError(variable_type))
	#
	def get_rq_vartype(self, dtype):
		"""
		Convert a string datatype (from numpy) into a string for the RQ file.
		_________________________
		Options:
			dtype
				A numpy datatype
		
		_________________________
		Returns:
			 Two values: RQ variable type and the number of bytes per datatype.
			 Note that this number is different for characters. Though a
			 character is only 1 byte, this function returns the number of
			 bytes per the entire string. So "hello" would return
			 "char", 5
		"""
		if dtype == 'float64':
			variable_type, bytenum = 'double', 8
		elif dtype == 'float32':
			variable_type, bytenum = 'single', 4
		elif dtype == 'uint8':
			variable_type, bytenum = 'uint8', 1
		elif dtype == 'uint16':
			variable_type, bytenum = 'ushort', 2
		elif dtype == 'uint32':
			variable_type, bytenum = 'uint32', 4
		elif dtype == 'uint64':
			variable_type, bytenum = 'ulonglong', 8
		elif dtype == 'int8':
			variable_type, bytenum = 'int8', 1
		elif dtype == 'int16':
			variable_type, bytenum = 'short', 2
		elif dtype == 'int32':
			variable_type, bytenum = 'int32', 4
		elif dtype == 'int64':
			variable_type, bytenum = 'longlong', 8
		elif dtype == 'bool':
			variable_type, bytenum = 'bool', 1
		elif 'string' in dtype:
			# character array
			variable_type = 'char'
			bytenum = int(re.split('\D+', dtype)[-1])/8
		else:
			raise(UnknownDataTypeError(dtype))
		#
		return variable_type, bytenum


#-----------------------------------------------------------------------------
######## Create an instance of the object that can be used right away. Of
# course the user is free to create their own instance.
dt = dtypes()


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
#		Exception definitions
#
class UnknownDataTypeError(Exception):
	"""
	Exception raised when the numpy datatype or RQ variable type are not 
	recognized.
	"""
	def __str__(self):
		messgstr = "\nUnknown variable type: %s" %(self.args[0])
		return messgstr

class UnknownEndianError(Exception):
	"""
	Raised when the endian name (little or big) is not recognized.
	"""
	def __str__(self):
		messgstr = "\nUnknown endian: %s" %(self.args[0])
		return messgstr

class BadEndianIndicatorError(Exception):
	"""
	Raised when the 4 bytes given as the endian check value do
	not match expectations.
	"""
	def __str__(self):
		messgstr = \
			"\nThe number %x is not understood as an endianness indicator" \
			%(self.args[0])
		return messgstr

