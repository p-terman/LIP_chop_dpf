"""
AdditionalFileFormat_SaveHDF5_PyMod.py

This module reads a .rq binary file and converts it to a .rq.hdf5 file.
It uses the HDF5 file format as described here:
   http://teacher.pas.rochester.edu:8080/wiki/bin/view/LUX/FileFormatDocumentation_RQHDF5

It uses RQReader_PyMod.ReadRQFile for reading the .rq file
then uses h5lux_PyMod.convert_rq2hdf5 for HDF5 conversion

20130703 CHF - Created
"""

from h5lux_PyMod import convert_rq2hdf5
import os.path

def AdditionalFileFormat_SaveHDF5(evt_file_name, evt_file_path, rq_file_name, rq_file_path, run_module_number, dp_settings_xml_path, lug_iqs_xml_path):
	"""
	This function converts .rq
	_________________________
	Inputs are self-explanatory according to DP standards for any module

	No outputs

	Versioning:
	    20130703 CHF - Created

	To do:
	_________________________
	"""
	# Bookkeeping
	rq_file_fullpath = rq_file_path + '/' + rq_file_name
	hdf5_output_dir = rq_file_path + '/' + 'hdf5'
	hdf5_file = rq_file_name.replace('.rq','.rq.hdf5')
	hdf5_file_fullpath = hdf5_output_dir + '/' + hdf5_file

	# Check that target file for conversion exists
	if not os.path.isfile(rq_file_fullpath):
		print "ERROR: File not found: %s" % rq_file_fullpath

	# If .../hdf5 subfolder doesn't exist, make it
	if not os.path.isdir(hdf5_output_dir):
		try:
			os.mkdir(hdf5_output_dir)
			print 'Created directory %s' % hdf5_output_dir
		except:
			print 'ERROR: Could not create directory %s' % hdf5_output_dir

	# Convert .rq to .rq.hdf5
	# It uses RQReader_PyMod.ReadRQFile for reading the .rq file
	# then uses h5lux_PyMod.convert_rq2hdf5 for HDF5 conversion
	try:
		convert_rq2hdf5(rq_file_fullpath,hdf5_file_fullpath)
		print 'Successfully converted %s to %s' % (rq_file_name,hdf5_file)
	except:
		print 'ERROR: Could not convert file'