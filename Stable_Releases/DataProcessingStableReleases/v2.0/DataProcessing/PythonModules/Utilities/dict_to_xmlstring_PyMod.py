# dict_to_xmlstring_PyMod
# Tomasz Biesiadzinski (TPB)
# July 20, 2014
#
# This is a utility module. recursive_string_constructor(cdict) returns a string
# with xml tags constructed from the input dictionary's keys and values from
# the dictionary.
#
# 20140722 TPB - Created

def recursive_xml_string_constructor(cdict):
	"""
	Construct an XML string from a dictionary.
	
	_________________________
	Options:
		cdict
			A dictionary that will be turned into a string with xml tags.
			The keys will become the openning and closing tags. The 
			dictionary values become the xml entries. 
			Can be a dictionary or ordered dictionary which extends dictionary.
	
	_________________________
	Returns:
		A string of xml values
	
	_________________________
	History:
		v1.0 20140720 TPB - Created.
	
	"""
	clist = []
	# Construct a list of tag and value strings
	recursive_xml_list_constructor(cdict,clist)
	# Convert the list into a string.
	settings_string = ''.join(clist)
	return settings_string

def recursive_xml_list_constructor(cdict, xml_entries_list):
	"""
	Construct an XML list from a dictionary.
	
	Have to create a list because strings are imutable. That is, you
	cannot extend one. Even using "+=" you'll get a copy instead of the
	original reference. I *could* pass a list with a single element, the
	string. Then that list would be passed as a refence and each time I
	made a new string it would be preserved. But that means copying the same
	list over and over again while adding a word to it each time. This isn't
	eficent. Instead I'll construct a list of words and then they can be 
	joined once.
	_________________________
	Options:
		cdict
			A dictionary that will be turned into a list with xml tags.
			The keys will become the openning and closing tags. The 
			dictionary values become the xml entries. Each list element is
			either an xml tag or entry.
			Can be a dictionary or ordered dictionary which extends dictionary.
		
		xml_entries_list
			A list that will have the xml data appended to it. It should be
			an empty list when the function is called.
	
	_________________________
	Returns:
		Function returns None. The constructed xml tags and values are appened
		to xml_entries_list that the user can use.
	
	_________________________
	History:
		v1.0 20140720 TPB - Created.
	
	"""
	if isinstance(cdict, dict):
		# the current entry is a dictionary
		initags = cdict.keys()
		for skey in initags:
			xml_entries_list.append("<%s>" %(skey))
			#settings_string += \
			recursive_xml_list_constructor(cdict[skey], \
				xml_entries_list)
			xml_entries_list.append("</%s>" %(skey))
	else:
		# The current entry is no longer a dictionary. It should be a string.
		xml_entries_list.append("%s" %(cdict))



