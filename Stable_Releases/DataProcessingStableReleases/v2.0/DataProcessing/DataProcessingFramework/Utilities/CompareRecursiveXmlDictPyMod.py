






import xml2dict

def CompareRecursiveXmlDict(a,b,different_flag = 0):  
    """ 
    
    This function will go through each level of an XmlDictObject (created by xml2dict)
    and will stop if any dictionary key or value is different between a and b.
    The function is recursive by calling itself over and over if there are dictionaries or lists
    until it reaches a value.
    
    If every element is identical, the value of different_flag will remain at 0. 
    
    20130120 CHF - Created
    
     """
    
    # It's identical until proven otherwise.
    
    # Set verbosity
    verbose = 1
    
    # Make sure this level has the same number of keys.
    # This way, we can just pick the keys from one, and compare to the other ones
    # without needing to do the reverse as well
    
    
    a_has_keys = hasattr(a, 'keys')
    b_has_keys = hasattr(b, 'keys')

    same_properties = a_has_keys * b_has_keys
    if not same_properties:
        if verbose == 1:
            print '*** Different key properties:'
            print 'a has keys? %d, b has keys? %d' % (a_has_keys,b_has_keys)
        different_flag = different_flag+1
        return different_flag
            
    if len(a.keys()) != len(b.keys()):
        if verbose == 1:
            print '*** Different number of keys:'
            print '%s vs %s' % (a.keys(),b.keys())
        different_flag = different_flag+1
        return different_flag
    
    # For every key
    for i in a.keys():
        if verbose == 1:
            print '%s |' % i,
        # Check if it also exists in b
        if i in b.keys():
            # Good, it exists. Now go recursive!
            # If this is an xml dictionary, then go deeper
            if isinstance(a[i],xml2dict.XmlDictObject):
                different_flag = CompareRecursiveXmlDict(a[i],b[i],different_flag) 
            # if this is a list, repeat for each element
            elif type(a[i]) is list:
                
                # check if the lists have the same lengths!
                if len(a[i]) != len(b[i]):
                    if verbose == 1:
                        print '*** Different list length:'
                        print '%s vs %s' % (len(a[i]),len(b[i]))
                    different_flag = different_flag+1
                    return different_flag
                
                if verbose == 1:
                    print '*** looping through module list ***'
                for qq in range(len(a[i])):
                    if different_flag == 0: # only continue if it hasn't gotten a discrepancy
                        if verbose == 1:
                            print '(%d) ' % qq,
                        different_flag = CompareRecursiveXmlDict(a[i][qq],b[i][qq],different_flag)
            else: # then it's holding a value. Compare value, which should be the same
                if a[i] != b[i]:
                    if different_flag == 0:
                        if verbose == 1:
                            print '*** Different value in key "%s" ***' % i
                            print '"%s" vs "%s"' % (a[i],b[i]) 
                        different_flag = different_flag+1
                        return different_flag
                else:
                    if verbose == 1:
                        print 'value: %s' % a[i]
        else:
            # If it doesn't exist, it's not the same!
            if verbose == 1:
                print '*** found discrepancy in %s ***' % i
            different_flag = different_flag + 1
            return different_flag
            
    return different_flag
    
    
    