# GetRecursiveIQ_PyMod.py
# Tomasz Biesiadzinski
# July 07, 2014

from xml2dict import xml2dict				# reading/parsing xml files

temp_tag = \
	"this_is_a_complicated_tag_that_is_unlikely_to_be_used_by_anything_else6724"

def GetRecursiveIQ(iq_xml_path, descending_xml_tags):
	"""
	Get the IQs in a generic fashion.
	
	_________________________
	Options:
		iq_xml_path
			The path to the constructed IQ file that contains all of the IQs
			needed for the processing. 
		descending_iq_tags
			List of tupples that contain the XML tags needed to get an IQ. For
			instance, the bottom PMT-derived single electron parameters are
			contained in lug_iqs_xml['iq'][i]['fit']['botpmt']. To get them,
			along with all PMTs, use:
			descending_iq_tags = [('fit','botpmt'), ('fit','allpmt')]
			The returned list will then contain two objects, one for the 
			bottom PMTs and one for all PMTs. In this specific example, each
			object will be a dictionary. 
	
	Optional
		
	
	_________________________
	Returns:
		A list with an item for each tupple in descending_iq_tags.
		If the content of the tags in descending_iq_tags is a value, then that
		value will be returned in the list. If that content is a dictionary,
		then that dictionary will be in that list. Finally, if multiple IQs
		match the tupple in descending_iq_tags, crash and burn.	
		
		We will attempt to convert all strings to floats.
		
		Note that an exception will NOT be raised if a tupple is not found. The
		returned value (in the list) for such an item will be None.
	"""
	# get the number of IQs requsted
	list_len = len(descending_xml_tags)
	# Check that the descending_xml_tags contains tupples (or sub-lists)
	# of tags.
	for j in range(list_len):
		# First check that each object in the list of XML tupples is itterable.
		# This is not a full check since strings will pass this test.
		try:
			nothing = len(descending_xml_tags[j])
		except(TypeError):
			print "%s is not itterable" %(descending_xml_tags[j].__repr__())
			raise(BadTagContainerError(descending_xml_tags[j].__repr__()))
		# Now check that it's not a string
		if isinstance(descending_xml_tags[j], str):
			print "%s is a string" %(descending_xml_tags[j].__repr__())
			raise(BadTagContainerError(descending_xml_tags[j].__repr__()))
	#
	# Load the IQ xml file. We need to prepend and append tags to make xml2dict
	# work.
	f=open(iq_xml_path)
	ntxt = "<%s>" %(temp_tag) + f.read() + "</%s>" %(temp_tag)
	f.close()
	lug_iqs_xml = xml2dict(ntxt)[temp_tag]
	del ntxt
	num_iqs = len(lug_iqs_xml['iq'])	# IQ number
	# initialize return list
	retlist = [None for i in range(list_len)]
	#
	# go through all of the IQs
	for i in range(num_iqs):
		for j in range(list_len):
			iq_not_here = False	# Initialize a flag that will keep track of
			# whether the desired tupple is present in the IQ number i.
			tag_not_here = False
			#
			# Current dictionary which we will address recursively.
			cdict = lug_iqs_xml['iq'][i]
			if not isinstance(cdict, dict):
				# not a dictionary. Means this is an empty IQ. 
				iq_not_here = True
				break
			for skey in descending_xml_tags[j]:
				try:
					if not cdict.has_key(skey):
						# the desired IQ is not here. We'll exit the loop over
						# the descending_xml_tags hence iq_not_here = True
						iq_not_here = True
						break
					cdict = cdict[skey]	# reset the dictionary to the next level
				except(AttributeError):
					# not a dictionary. Means this is an empty IQ.
					tag_not_here = True	# Don't exit the loop over the
										# the descending_xml_tags but don't save
					break
			if iq_not_here:	# the deisred IQ is not in this one
				continue
			# else: The tupple is here
			if retlist[j] != None:
				# this value was already found
				raise(MultipleIQsError('%s'\
					%(descending_xml_tags[j].__repr__())))
			#
			# We've reached the bottom of the dictionary tree and can now
			# save the IQ value or dictionary. Try to convert to a float if
			# not a dictionary.
			if not tag_not_here:
				if isinstance(cdict, dict):
					retlist[j] = cdict	# a dictionary
					# attempt to convert all to 
					tdict = retlist[j]
					for sskey in tdict.keys():
						if isinstance(tdict[sskey], str):
							# this entry is a string
							try:	# convert to a float
								tdict[sskey] = float(tdict[sskey])
							except(ValueError):	# return the string instead
								pass
							# done with this one
							continue
						else:	# continues to be a dictionary
							tdict = tdict[sskey]
				else:	# not a dictionary
					try:	# convert to a float
						retlist[j] = float(cdict)
					except(ValueError):	# return the string instead
						retlist[j] = cdict
	#
	# Done. Return the list
	return retlist


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
#		Exception definitions
#

class BadTagContainerError(Exception):
	"""
	Raised when the container holding xml tags is not itterable (and not a
	string). Possible if someone forgot to put the tag container into a higher
	level list.
	"""
	def __str__(self):
		messgstr = "\n%s is not an XML tag container." %(self.args[0])
		return messgstr

class MultipleIQsError(Exception):
	"""
	Raised when there are multiple IQs matching the same requested tupple.
	"""
	def __str__(self):
		messgstr = "\nMultiple values found for tags %s" %(self.args[0])
		return messgstr

