# EnergyReconstruction_Impulsive_PyMod
# Tomasz Biesiadzinski (TPB)
# July 07, 2014
#
# This is a conversion of the matlab code EnergyReconstruction_Naive.m by 
# Patrick Phelps (PHP).
# The code functionality was coppied essentially line by line with
# small changes and bug fixes.


from RQReader_PyMod import ReadRQFile			# reading in RQ files
from xml2dict import xml2dict					# reading/parsing xml files
# function for Data Processing (DP) settings
from GetDPModuleParameters_PyMod import GetDPModuleParameters
# function for getting IQs
from GetRecursiveIQ_PyMod import GetRecursiveIQ
import numpy as np								# mathematical operations
import os										# file location handling
import sys								# output of exit status to the system
import RQWriter_PyMod as rqw					# RQ file writing
from collections import OrderedDict as odict	# ordered keys in a dictionary

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# Constant definitions

constant_maxnumber_dict_key = 'max_num_pulses'

parameter_dict_key_NR_conversion_A	= 'NR_conversion_A'
parameter_dict_key_NR_conversion_B	= 'NR_conversion_B'
parameter_dict_key_extractioneff	= 'extractioneff'
parameter_dict_key_W				= 'W'
parameter_dict_key_pde				= 'pde'
parameter_dict_key_single_e_area_bot= 'single_e_area_bot'
parameter_dict_key_single_e_area	= 'single_e_area'


#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

#_______________________________________________________________________________
def EnergyReconstruction_Impulsive(mymodule_settings, global_dp_tags, \
			iqlist, rq_dict, module_name = 'EnergyReconstruction_Impulsive', \
			verbose=True):
	"""
	This module computes the estimate of event-generated quanta; the number of
	photons (corrected for photon detection efficienty) and the number of
	electrons (corrected for extraction efficiency). It also computes the
	estimated energy for single scatter events. 

	According to the original Matlab version of this module, 
	EnergyReconstruction_Naive.m:
		This module will do a very naive energy reconstruction
		1) Calibrate the n_gamma in the S1 via the g1 / PDE found off the
		164 keV Xenon activation line (If the fields change this needs
		a recalibration)
		2) Using an estimate of the 1eS2 size and the extraction efficiency
		estimate the S2 n_e
	
	_________________________
	Options:
		mymodule_settings
		
		dp_settings_xml_dp
		
		iqlist
		
		rq_dict
			A dictionary of reduced quantity (RQ) arrays.
			Required RQs:
			* pulse_classification - Need to know which pulsea are S1 and S2s
			* xyz_corrected_pulse_area_all_phe
			* xyz_corrected_pulse_area_bot_phe
			
	Optional
		module_name = 'EnergyReconstruction_Impulsive'
			The module name that will be printed.
	
	_________________________
	Returns:
		newrq
			A dictionary with newly created reduced quantities (RQ) arrays.
		
			According to the original Matlab version of this module
			(EnergyReconstruction_Naive.m) the RQs created are:
			* energy_found_1eS2_flag (binary) - Did this dataset find a 1eS2
			 size out of the LUG (1) or did it use the default numbers (0).
			* energy_bottom_poor_fit (binary) - Was the bottom array 1eS2
			 considered a poor fit (1) or not (0) based on was it within a
			 window of 10.5 phe the window is defined based on the difference
			 between 1eS2 area (all) and average 1eS2 area (24.5 phe)
			* num_electrons (num_pulses,num_evts). Number of liquid electrons
			 in all S2 pulses (as defined by pulse_classification rq).
			 Outputs 0 for all non S2 pulses.  Done over all PMTs.
			* num_electrons_bot (num_pulses,num_evts). Number of liquid
			 electrons in all S2 pulses (as defined by pulse_classification
			 rq).  Outputs 0 for all non S2 pulses.  Done over bot PMTs.
			* num_photons (num_pulses,num_evts).  Number of produced photons
			 in all S1 pulses (as defined by pulse_classificaiton rq.) Outputs
			 0 for all non S1 pulses. Done over all PMTs.
			* energy_keVee_all (num_evts).  Energy reconstructed for single
			 scatter events from all PMTs (keVee).
			* energy_keVee_bot (num_evts).  Energy reconstructed for single
			 scatter events from just bottom PMTs (keVee).
			* energy_keVnr_all (num_evts). Energy reconstructed for single
			 scatter events from all PMTs (keVnr).
			* energy_keVnr_bot (num_evts). Energy reconstructed for single
			 scatter events from just the bottom PMTs (keVnr).
			* refined_energy_all (num_evts). Refined energy reconstructed for
			 single scatter events from all PMTS (keVee)
			* refined_energy_bot (num_evts). Refined energy reconstructed for
			 single scatter events from bot PMTS (keVee)
	
	_________________________
	History:
		
		v1.0 20140708 TPB - Created from EnergyReconstruction_Naive.m
			Fixed a bug with the keVnr energy where energy_keVnr_all RQ was
			being overwritten by energy_keVnr_bot and the energy_keVnr_bot
			RQ was then not set.
		
		Original EnergyReconstruction_Naive.m:
			v2.0 20130906 PHP - Add refined energy and work towards 1eS2
				size = IQ
			v1.0 20130702 PHP - Created
	
		RQ versions:
			v1.0 num_electrons
			v1.0 num_photons
			v1.0 energy_keVee_all
			v1.0 energy_keVee_bot
			v1.0 energy_keVnr_all
			v1.0 energy_keVnr_bot
	"""
	if verbose: print '\n\n *** Starting module %s\n' %(module_name)
	# initialize the ordered dictionary variable that will hold the new RQs
	newrq = odict()
	# get the maximum number of pulses from the global tag of the DP
	# settings xml file
	max_num_pulses = int(global_dp_tags[constant_maxnumber_dict_key])
	#
	# Count the number of events and pulses
	pulse_class_shape = rq_dict['pulse_classification'].shape
	evtnum, pulsenum = pulse_class_shape
	## In principle, rq_dict['event_number'].shape should be 1D and equal 
	## to evtnum. But check just in case. ---- Nope. Don't bother.
	#if evtnum != rq_dict['event_number'].shape[0]:
	#	raise(EventNumberMissMatchError('pulse_classification', \
	#		'event_number', rq_dict['pulse_classification'].shape.__repr__(), \
	#		rq_dict['event_number'].shape.__repr__()))
	#
	# Assign the IQs to various variables.
	# require that all of the IQs reqeusted are found
	if iqlist[0] != None and iqlist[1] != None and iqlist[2] != None:
		single_e_area_bot = iqlist[0]['mean']
		if iqlist[0].has_key('adjrsquare'): adjRsqr = iqlist[0]['adjrsquare']
		else: adjRsqr = 1.	# fit quality undetermined.
		single_e_area = iqlist[1]['mean']
		filename_prefixs = iqlist[2]
		found_it = True		# flag that IQs were found
		newrq['energy_found_1eS2_flag'] = \
			np.zeros(evtnum, dtype=np.bool) + True
		newrq['energy_bottom_poor_fit'] = np.zeros(evtnum, dtype=np.bool)
	else: found_it = False	# glag that IQs were not found
	#
	# Check things are kosher
	if not found_it:
		# We didn't find our 1eS2 size
		newrq['energy_found_1eS2_flag'] = np.zeros(evtnum, dtype=np.bool)
		newrq['energy_bottom_poor_fit'] = np.zeros(evtnum, dtype=np.bool)
		single_e_area = \
			mymodule_settings[parameter_dict_key_single_e_area]
		single_e_area_bot = \
			mymodule_settings[parameter_dict_key_single_e_area_bot]
	elif adjRsqr < 0.97:
		single_e_area_bot = \
			mymodule_settings[parameter_dict_key_single_e_area_bot]
		newrq['energy_bottom_poor_fit'] = \
			np.zeros(evtnum, dtype=np.bool) + True
	#
	# Identify single scatter events
	s1flags = np.where(rq_dict['pulse_classification'] == 1, 1., 0.)
	s2flags = np.where(rq_dict['pulse_classification'] == 2, 1., 0.)
	single_scatter = \
		np.where((s2flags.sum(axis=1) * s1flags.sum(axis=1)) == 1,1,0)
	## assign the truth value to all pulses in an event
	#single_scatter = \
	#	np.repeat(single_scatter,pulsenum).reshape(pulse_class_shape)
	#
	# Calculate Quanta
	newrq['num_photons'] = \
		s1flags * rq_dict['xyz_corrected_pulse_area_all_phe'] / \
			mymodule_settings[parameter_dict_key_pde]
	newrq['num_electrons'] = \
		s2flags * rq_dict['xyz_corrected_pulse_area_all_phe'] / \
			(single_e_area*mymodule_settings[parameter_dict_key_extractioneff])
	newrq['num_electrons_bot'] = \
		s2flags * rq_dict['xyz_corrected_pulse_area_bot_phe'] / \
		(single_e_area_bot*mymodule_settings[parameter_dict_key_extractioneff])
	# Also compute the totals for each event
	total_num_photons = newrq['num_photons'].sum(axis=1)
	total_num_electrons = newrq['num_electrons'].sum(axis=1)
	total_num_electrons_bot = newrq['num_electrons_bot'].sum(axis=1)
	#
	# Calculate Energy
	newrq['energy_keVee_all'] = (total_num_electrons+total_num_photons) * \
		single_scatter * mymodule_settings[parameter_dict_key_W]
	newrq['energy_keVee_bot'] = (total_num_electrons_bot+total_num_photons) * \
		single_scatter * mymodule_settings[parameter_dict_key_W]
	temp_exponent = 1. / mymodule_settings[parameter_dict_key_NR_conversion_B]
	newrq['energy_keVnr_all'] = \
		(newrq['energy_keVee_all'] / \
			mymodule_settings[parameter_dict_key_NR_conversion_A]) ** \
			temp_exponent
	newrq['energy_keVnr_bot'] = \
		(newrq['energy_keVee_bot'] / \
			mymodule_settings[parameter_dict_key_NR_conversion_A]) ** \
			temp_exponent
	#
	# Refine Energy (keVee)
	newrq['refined_energy_all'] = (Beta(newrq['energy_keVee_all']) * \
		total_num_photons + total_num_electrons) * single_scatter * \
			mymodule_settings[parameter_dict_key_W] * \
			Alpha(newrq['energy_keVee_all'])
	newrq['refined_energy_bot'] = (Beta(newrq['energy_keVee_bot']) * \
		total_num_photons + total_num_electrons_bot) * single_scatter * \
			mymodule_settings[parameter_dict_key_W] * \
			Alpha(newrq['energy_keVee_bot'])
	#
	# make sure all of the RQs are doubles
	for skey in newrq.keys(): newrq[skey] = np.float64(newrq[skey])
	# Return the newly computed RQs
	return newrq

#_______________________________________________________________________________
def Beta(raw_energy):
	"""
	A helper function for refined_energy_all and refined_energy_bot.
	
	_________________________
	Options:
		raw_energy
			The previously computed energy.
	
	Optional
	_________________________
	Returns:
		An array of values of the same dimmensions as raw_energy.
	"""
	return 0.83827 - 0.83336 * np.exp(-0.46154 * raw_energy)

#_______________________________________________________________________________
def Alpha(raw_energy):
	"""
	A helper function for refined_energy_all and refined_energy_bot.
	
	_________________________
	Options:
		raw_energy
			The previously computed energy.
	
	Optional
	_________________________
	Returns:
		An array of values of the same dimmensions as raw_energy.
	"""
	beta_raw_energy = Beta(raw_energy)	# pre-compute since we'll use it a
										# few times.
	beta_raw_energy_sq = beta_raw_energy**2.	# pre-compute to pottentially
												# make this faster								
	return 1.0504 + 1.0054*beta_raw_energy - 1.8906*beta_raw_energy_sq + \
		1.277*beta_raw_energy*beta_raw_energy_sq - \
		0.41817*beta_raw_energy_sq**2.

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
#		Exception definitions
#
class EventNumberMissMatchError(Exception):
	"""
	Raised when the numbers of events in an RQ file obtained from two different
	RQs disagree. This should never happen.
	This exception is here only because I explicitly force the event numbers to
	be the same while the matlab version initialized variables with two sets
	of event numbers. That being said, the matlab version then did math with
	these so the dimmensions were implicitly assumed to be the same.
	"""
	def __str__(self):
		messgstr = "\n%s and %s RQs have different event numbers." \
			%(self.args[0], self.args[1])
		messgstr = "The shape of %s is %s and the shape of %s is %s" \
			%(self.args[0], self.args[2], self.args[1], self.args[3])
		return messgstr


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

# main "function"
if __name__ == "__main__":
	# check the arguments
	if len(sys.argv) == 1:
		print 'Usage: python EnergyReconstruction_Impulsive_PyMod.py',
		print 'evt_file_name evt_data_path rq_file_name rq_data_path',
		print 'data_processing_xml_path iq_xml_path'
		sys.exit(0)	# 0 means failure
	if len(sys.argv) < 7:
		print "Not all of the required arguments were provided."
		print 'Usage: python EnergyReconstruction_Impulsive_PyMod.py',
		print 'evt_file_name evt_data_path rq_file_name rq_data_path',
		print 'data_processing_xml_path iq_xml_path'
		sys.exit(0) # 0 means failure
	else:
		# parse the command line arguments into filenames and paths.
		filename_evt,data_path_evt,filename_rq,data_path_rq, \
			data_processing_xml_path,iq_xml_path = sys.argv[1:]
		# filename_evt and data_path_evt are not used but are expected for
		# compatibility issues.
	#
	# Define the module name
	#module_name = 'EnergyReconstruction_Impulsive'
	module_name = 'EnergyReconstruction_Naive'		# temporary so that we can
				# use the existing DP settings file
	#
	# get the DP parameters. Also return the global setting tag content
	# in addition to the module settings this function
	mymodule_settings, global_dp_tags = \
		GetDPModuleParameters(data_processing_xml_path, module_name)
	#
	# Get the IQs
	iqlist = GetRecursiveIQ(iq_xml_path, [('fit', 'botpmt'), ('fit','allpmt')])
	#
	# Read in the existing RQ file
	filename = os.path.join(data_path_rq,filename_rq)
	rq_dict, settings_dict, livetime_dict, file_dict = ReadRQFile(filename, \
		num_pulses=global_dp_tags['max_num_pulses'])
	#
	# calculate the quanta and energies
	newrq = EnergyReconstruction_Impulsive(mymodule_settings, \
		global_dp_tags, iqlist, rq_dict)#, module_name=module_name)
	# add the new RQs to the old ones
	for skey in newrq.keys(): rq_dict[skey] = newrq[skey]
	#
	# Write Output File
	status = rqw.RQWriter(filename_rq, rq_dict, file_dict, settings_dict, \
		livetime_dict)
	sys.exit(status)


"""
## command line
#python EnergyReconstruction_Impulsive_PyMod.py 0 0 lux10_20131218T1713_f000001_cp09205.rq . dpcp9205.xml iqtemp.xml

# --------------------------------------
# test in shell

import RQReader_PyMod as rrq
import GetDPModuleParameters_PyMod as gdp
import GetRecursiveIQ_PyMod as griq
import EnergyReconstruction_Impulsive_PyMod as eri
import os

# temporary settings
data_path_rq = '/home/tomaszbi/lux/lux_data/tritium/rq/lux10_20131218T1713_cp09205/binary'
filename_rq = 'lux10_20131218T1713_f000001_cp09205.rq'
data_processing_xml_path = 'dpcp9205.xml'
iq_xml_path = 'iqtemp.xml'
module_name = 'EnergyReconstruction_Naive'

reload(griq)
iqlist = griq.GetRecursiveIQ(iq_xml_path, [('fit', 'botpmt'), ('fit','allpmt')])

reload(gdp)
mymodule_settings, global_dp_tags = gdp.GetDPModuleParameters(data_processing_xml_path, module_name)

reload(rrq)
rq_dict, settings_dict, livetime_dict, header_dict = rrq.ReadRQFile(os.path.join(data_path_rq,filename_rq))

reload(eri)
newrq = eri.EnergyReconstruction_Impulsive(mymodule_settings, global_dp_tags, iqlist, rq_dict)#, module_name=module_name)


# test this
rqnames = newrq.keys()

print (rq_dict['num_photons'] - newrq['num_photons'])/rq_dict['num_photons']

import matplotlib.pyplot as pyp

a=pyp.plot((rq_dict['num_photons'] - newrq['num_photons'])/rq_dict['num_photons'],'.')
pyp.show()

a=pyp.plot((rq_dict['num_electrons'] - newrq['num_electrons'])/rq_dict['num_electrons'],'.')
pyp.show()

a=pyp.plot((rq_dict['num_electrons_bot'] - newrq['num_electrons_bot'])/rq_dict['num_electrons_bot'],'.')
pyp.show()

a=pyp.plot((rq_dict['energy_keVee_all'] - newrq['energy_keVee_all'])/rq_dict['energy_keVee_all'], '.')
pyp.show()

#rq_dict['energy_keVnr_all'] is actually rq_dict['energy_keVnr_bot'] because of a bug in EnergyReconstruction_Naive.m
a=pyp.plot((rq_dict['energy_keVnr_all'] - newrq['energy_keVnr_bot'])/rq_dict['energy_keVnr_all'], '.')
pyp.show()

#rq_dict['energy_keVnr_bot'] is 0 because of a bug in EnergyReconstruction_Naive.m
#a=pyp.plot((rq_dict['energy_keVnr_bot'] - newrq['energy_keVnr_bot'])/rq_dict['energy_keVnr_bot'], '.')
#pyp.show()


a=pyp.plot((rq_dict['refined_energy_all'] - newrq['refined_energy_all'])/rq_dict['refined_energy_all'], '.')
pyp.show()

a=pyp.plot((rq_dict['refined_energy_bot'] - newrq['refined_energy_bot'])/rq_dict['refined_energy_bot'], '.')
pyp.show()





"""

