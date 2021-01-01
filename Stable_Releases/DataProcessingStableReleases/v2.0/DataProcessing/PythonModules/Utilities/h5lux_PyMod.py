# h5lux_PyMod
# CH Faham, James Verbus
# April 9, 2014

# You can load an entire .rq.hdf5 file with load_rq_file:
#  rq,settings,livetime = h5lux_PyMod.load_rq_file('lux10_20130505T0712_f000001_cp05457.rq.hdf5')
#
# You can write an additional rq variable to the .rq.hdf5 file:
#  new_rqs = { 'testing1':np.zeros((5,100)), 'testing2':np.ones((10,50)) }
#  h5lux_PyMod.write_rqs('lux10_20130505T0712_f000003_cp05457.rq.hdf5',new_rqs)
#
# You can convert a .rq.mat file to a .rq.hdf5 file:
#  h5lux_PyMod.convert_rq_mat2hdf5('lux10_20130505T0712_f000001_cp05457.rq.mat'
#
#
import h5py
import numpy as np
import os.path

#______________________________________________________________________________
#******************************************************************************
#                           FILE LOADING FUNCTIONS
#******************************************************************************
#
# load_evt_file
# load_rq_file

def load_evt_file(hdf5_filename):
    """
    This function loads all of the contents in a .evt.hdf5 file.
    _________________________
    Inputs:
        hdf5_filename
            The .evt.hdf5 file name
    Returns:
        evt
            Dictionary with evt fields
        settings
            Dictionary with evt,daq,trigger settings

    Versioning:
        20140505 JRV - Created

    To do:
    _________________________
    """

    # Open file
    f = h5py.File(hdf5_filename, 'r')

    # Get settings
    settings = load_attributes_recursive(f['/settings'])

    # Get livetime
    livetime = load_livetime(f)

    # Make a dictionary of all the field names
    evt_names_keys = dict(f['/'].items()).keys()

    # Remove the settings names so loop skips them
    evt_names_keys.remove('evt_settings')
    evt_names_keys.remove('daq_settings')
    evt_names_keys.remove('trigger_settings')
    evt_names_keys.remove('livetime')

    # Loop for each evt field
    evt = dict()
    for evt_name in evt_names_keys:
        evt[evt_name] = load_datasets_recursive(f[evt_name])

    f.close()

    return evt, settings, livetime

#______________________________________________________________________________


def load_rq_file(hdf5_filename):
    """
    This function loads all of the contents in a .rq.hdf5 file.
    _________________________
    Inputs:
        hdf5_filename
            The .rq.hdf5 file name
    Returns:
        rq
            Dictionary with rqs
        settings
            Dictionary with evt,daq,trigger settings

    Versioning:
        20140409 CHF - Created
        20140617 CHF - Changed to use /rqs/ and /settings/ groups

    To do:
        Error handling
        Load individual rqs from list
        Optional flag to load settings or not
    _________________________
    """

    # Open file
    f = h5py.File(hdf5_filename, 'r')

    # Get settings
    settings = load_attributes_recursive(f['/settings'])

    # Get livetime
    livetime = load_livetime(f)

    # Make a dictionary of all the RQ names
    rq_names_keys = dict(f['/rqs'].items()).keys()

    # Loop for each RQ
    rq = dict()
    for rq_name in rq_names_keys:
        rq[rq_name] = f['/rqs/%s'%rq_name][...]

    f.close()

    return rq, settings, livetime

#______________________________________________________________________________

#******************************************************************************
#                          FILE CONVERSION FUNCTIONS
#******************************************************************************
#
# convert_rq2hdf5
# convert_rq_mat2hdf5

def convert_rq2hdf5(rq_filename, hdf5_filename=None):
    """
    This function converts a .rq binary file into a .rq.hdf5 file.
    _________________________
    Inputs:
        rq_filename
            The .rq file name. Can be .rq.gz
    Optional arguments:
        hdf5_filename
            The .rq.hdf5 file name to save.
            If empty, default will be rq_filename with .rq.hdf5 extension

    Returns:
        None

    Versioning:
        20140630 CHF - Created
        20140723 CHF - Changed ReadRQFile line to include 'header' output as per 
                        Tomasz and Wing's updates to rq reader outputs

    To do:
        Error handling, verbosity arguments
    _________________________
    """
    # Need this to read .rq file.
    from RQReader_PyMod import ReadRQFile

    # Load .rq.mat file
    rq,settings_dict,livetime,header = ReadRQFile(rq_filename)

    # Create output file. If no output name provided, use mat file with hdf5
    if not hdf5_filename:
        hdf5_filename = rq_filename.replace('.rq', '.rq.hdf5')

    if os.path.exists(hdf5_filename):
        os.remove(hdf5_filename)
        print '*******************************************************'
        print '*** File already exists, deleting to create new one ***'
        print '*******************************************************'

    f = h5py.File(hdf5_filename)

    # Create a groupfor settings
    # Loop for each attribute and write to group
    settings_group = f.create_group("settings")

    write_attributes_recursive(settings_dict['settings'], settings_group)

    # Write livetime
    livetime_group = f.create_group("livetime")
    livetime_group['livetime_latch_samples'] = livetime['livetime_latch_samples']
    livetime_group['livetime_end_samples'] = livetime['livetime_end_samples']

    # Create list of rq names, remove settings and other __bookkeeping_
    rq_names_keys = rq.keys()
    
    # We don't want livetime in two places (it's its own group, not an RQ!)
    if 'livetime_latch_samples' in rq_names_keys:
        rq_names_keys.remove('livetime_latch_samples')
    if 'livetime_end_samples' in rq_names_keys:
        rq_names_keys.remove('livetime_end_samples')

    # Create RQs group
    rqs_group = f.create_group('rqs')

    # Loop for each rq
    for rq_name in rq_names_keys:

        rq_data = rq[rq_name]

        dset = rqs_group.create_dataset(rq_name, data=rq_data, compression=9, fletcher32=True)            

        print 'Done with %s' % rq_name

    print 'Finished converting file. Output file produced: %s' % hdf5_filename

    f.flush()
    f.close()

#______________________________________________________________________________

def convert_rq_mat2hdf5(mat_filename, hdf5_filename=None):
    """
    This function converts a .rq.mat file into a .rq.hdf5 file.
    _________________________
    Inputs:
        mat_filename
            The .rq.mat file name.
    Optional arguments:
        hdf5_filename
            The .rq.hdf5 file name to save.
            If empty, default will be mat_filename with .rq.hdf5 extension

    Returns:
        None

    Versioning:
        20140409 CHF - Created
        20140617 CHF - Changed to use /settings/ group
        20140620 CHF - Now uses /rqs, /livetime groups

    To do:
        Error handling, verbosity arguments, check if file exists,
        rewrite option
    _________________________
    """
    # Need this to read .rq.mat file. Only pre-v7.3 mat files supported
    from scipy.io import matlab

    # Load .rq.mat file
    m = matlab.loadmat(mat_filename, squeeze_me=True, struct_as_record=False)

    # Create output file. If no output name provided, use mat file with hdf5
    if not hdf5_filename:
        hdf5_filename = mat_filename.replace('mat', 'hdf5')

    if os.path.exists(hdf5_filename):
        os.remove(hdf5_filename)
        print '*******************************************************'
        print '*** File already exists, deleting to create new one ***'
        print '*******************************************************'

    f = h5py.File(hdf5_filename)

    # Get settings
    admin = m['admin']

    # Now do the settings
    # Create a group, one for each settings type
    settings_group = f.create_group("settings")
    evt_settings_group = settings_group.create_group("evt_settings")
    daq_settings_group = settings_group.create_group("daq_settings")
    trigger_settings_group = settings_group.create_group("trigger_settings")
    daq_settings_global_group = daq_settings_group.create_group("global")
    daq_settings_sis3301_group = daq_settings_group.create_group("sis3301")

    # Need to do this since 'global' is a reserved name (FAIL).
    daq_temp1 = getattr(admin.daq_settings, 'global')  # safe way to call 'global'
    daq_temp2 = getattr(admin.daq_settings.sis3301, 'global')  # safe way to call 'global'

    # Convert settings from mat structure format to dictionary
    evt_settings_dict = admin.evt_settings.__dict__
    daq_settings_global_dict = daq_temp1.__dict__
    daq_settings_sis3301_dict = daq_temp2.__dict__
    trigger_settings_dict_temp = admin.daq_settings.LUXTriggerSettings.__dict__

    if '_fieldnames' in evt_settings_dict:
        evt_settings_dict.pop('_fieldnames')
    if '_fieldnames' in daq_settings_global_dict:
        daq_settings_global_dict.pop('_fieldnames')
    if '_fieldnames' in daq_settings_sis3301_dict:
        daq_settings_sis3301_dict.pop('_fieldnames')
    if '_fieldnames' in trigger_settings_dict_temp:
        trigger_settings_dict_temp.pop('_fieldnames')

    # For some reason, the LUXTriggerSettings data is still in XML format (not dictionary).
    # Check if it's the case, and convert to dictionary if so.
    cl = trigger_settings_dict_temp['TriggerBuilder'].__class__
    if (cl == 'unicode') or (cl == 'str'):
        # Want to keep this to be loaded only if needed. Otherwise path not found will start giving errors 
        # and we otherwise don't need this!
        from xml2dict import xml2dict
        new_dict = xml2dict(str('<TriggerBuilder>'+trigger_settings_dict_temp['TriggerBuilder']+'</TriggerBuilder>'))
        # Yes, for some reason xml2dict requires a high-level wrap-around tag (TriggerBuilder)
        # We add it and then remove it in this next line
        trigger_settings_dict = new_dict['TriggerBuilder']
    else:
        trigger_settings_dict = trigger_settings_dict_temp['TriggerBuilder']

    settings_dict = admin.__dict__
    if '_fieldnames' in settings_dict:
        settings_dict.remove('_fieldnames')

    # Loop for each attribute and write to group
    write_attributes_recursive(settings_dict, settings_group)

    # Write livetime
    livetime_group = f.create_group("livetime")
    livetime_group['livetime_latch_samples'] = m['livetime_latch_samples']
    livetime_group['livetime_end_samples'] = m['livetime_end_samples']

    # Create list of rq names, remove settings and other __bookkeeping_
    rq_names_keys = m.keys()
    keys_to_remove = ['admin','source_filename','__header__','__version__','__globals__']
    for k in keys_to_remove:
        if k in rq_names_keys:
            rq_names_keys.remove(k)
    
    # We don't want livetime in two places (it's its own group, not an RQ!)
    if 'livetime_latch_samples' in rq_names_keys:
        rq_names_keys.remove('livetime_latch_samples')
    if 'livetime_end_samples' in rq_names_keys:
        rq_names_keys.remove('livetime_end_samples')

    # Create RQs group
    rqs_group = f.create_group("rqs")

    # Loop for each rq
    for rq_name in rq_names_keys:

        rq_data = m[rq_name]

        dset = rqs_group.create_dataset(rq_name, data=rq_data, compression=9, fletcher32=True)            

        print 'Done with %s' % rq_name


    print 'Finished converting file. Output file produced: %s' % hdf5_filename

    f.flush()
    f.close()

#______________________________________________________________________________

#******************************************************************************
#                              AUXILIARY FUNCTIONS
#******************************************************************************
#
# load_livetime
# write_rqs
# load_attributes_recursive
# write_attributes_recursive
# load_datasets_recursive
# nested_set

def load_livetime(f):
    """
    This function loads the livetime from any .hdf5 file
    _________________________
    Inputs:
        f
            Object containing open HDF5 file, e.g. f = h5py.File(filename)
    Returns:
        livetime
            Dictionary with latch_samples, end_samples for livetime calculation

    Versioning:
        20140415 CHF - Created
        20140702 CHF - Changed 'end_samples' to 'livetime_end_samples' and
                        'latch_samples' to 'livetime_latch_samples'

    To do:
        Error handling. Should give option to use filename instead of file
        object for stand-alone use
    _________________________
    """

    # Read settings
    livetime = dict()

    livetime_group = f['/livetime']
    livetime['livetime_latch_samples'] = livetime_group['livetime_latch_samples'][...]
    livetime['livetime_end_samples'] = livetime_group['livetime_end_samples'][...]

    return livetime

#______________________________________________________________________________

def write_rqs(hdf5_filename, rq_dict):
    """
    This function writes a single rq to a .rq.hdf5 file
    _________________________
    Inputs:
        hdf5_filename
            The .rq.hdf5 file name to write to
        rq_dict
            Dictionary containing new rqs to be written
    Returns:
        None

    Versioning:
        20140409 CHF - Created

    To do:
        Error handling, better feedback for a successful write
    _________________________
    """

    # Open file
    f = h5py.File(hdf5_filename)

    # If rqs group does not exist, make it
    if not f.get('rqs'):
        rqs_group = f.create_group('rqs')
    else:
        rqs_group = f['/rqs']

    # Write rq as a dataset. Loop for each key
    for rq_name in rq_dict.keys():

        rq_data = rq_dict[rq_name]

        # Default is to use max gzip compression (9) and checksum (fletcher32).
        dset = rqs_group['/rqs'].create_dataset(rq_name, data=rq_data, compression=9, fletcher32=True)

        print 'Wrote %s in %s' % (rq_name, hdf5_filename)

    f.flush()
    f.close()

#______________________________________________________________________________


def load_attributes_recursive(settings_group, level=0):
    """
    This function recursively traverses an hdf5 file's group structure
    and returns the contents as a dictionary.
    _________________________
    Inputs:
        settings_group
            The .hdf5 settings group, f['/settings']
        level
            Recursion level, leave blank
    Returns:
        output
            The output Python dictionary with the settings structure and values

    Versioning:
        20140701 CHF - Created

    To do:
    _________________________
    """
    # Initialize output (this is esp. needed in recursion so that same dict is not being nested)
    output = dict()

    # Grab fieldnames and group references, turn into dictionary for easy access
    names = [s[0] for s in settings_group.items()]
    groups = [g[1] for g in settings_group.items()]
    groups_dict = dict(zip(names,groups))

    level+=1
    for var in names:
        # If it's got more groups underneath... recurse
        print '\n' + '.'*level + '%s:' % var,
        if len(groups_dict[var].items()) > 0:
            output[var] = load_attributes_recursive(groups_dict[var], level)
        # Else it only has attributes inside!
        else:
            # Grab attribute names and values, turn into dictionary for easy access
            attr_names = [a[0] for a in groups_dict[var].attrs.items()]
            attr_values = [a[1] for a in groups_dict[var].attrs.items()]
            attr_dict = dict(zip(attr_names,attr_values))

            for attr in attr_names:
                print ' [%s] ' % attr,
                output[var] = attr_dict[attr]

    return output
#______________________________________________________________________________

def write_attributes_recursive(settings_dict, hdf5_group, level=0):
    """
    Auxiliary function that loops through each settings field inside dictionary
    (root_level) and writes an attribute in HDF5 group object.
    _________________________
    Inputs:
        settings_dict
            Dictionary with settings to loop through
        hdf5_group
            HDF5 group object where the attributes need to be written to

    Returns:
        None

    Versioning:
        20140409 CHF - Created
        20140701 CHF - Changed function to make it completely generic. It now recursively 
                        walks down through all sub-branches of a dictionary and creates
                        groups and attributes accordingly.
                        Changed name from attr_loop to write_attributes_recursive

    _________________________
    """

    level += 1
    # If it's a dictionary
    if isinstance(settings_dict,dict):
        # Loop for each key
        for fname in settings_dict.keys():
            # If this key holds a value that is a dictionary
            # (i.e. nested dictionary)
            if isinstance(settings_dict[fname],dict):
                print '\n' + '.'*level + '%s:' % fname,
                # Make a subgroup
                subgroup = hdf5_group.create_group(fname)
                # Perform a recursive call to this function
                write_attributes_recursive(settings_dict[fname], subgroup,level=level)
            # If it's not a nested dictionary
            else:
                print ' [%s] ' % fname,
                attr = settings_dict[fname]
                # Write the attribute
                hdf5_group.attrs[fname] = attr

#______________________________________________________________________________

def load_datasets_recursive(initem, output=None, inkeys=None):
    """
    This function recursively traverses an hdf5 file's group structure
    and returns the contents as a dictionary.
    _________________________
    Inputs:
        initem
            The .hdf5 group, ex. f['events']
        output
            The output dict. Best leave blank unless expert.
        inkeys
            The input key list. Best leave blank unless expert.
    Returns:
        output
            The output Python dictionary.

    Versioning:
        20140516 JRV - Created

    To do:
    _________________________
    """
    if output is None:
        output = dict()
    if inkeys is None:
        inkeys = list()
    keys = list(inkeys)

    for item in initem.items():
        print item
        if isinstance(item[1], h5py.Group):
            keys.append(item[0])
            load_datasets_recursive(item[1], output, keys)
            keys.pop()
        elif isinstance(item[1], h5py.Dataset):
            keys.append(item[0])
            nested_set(output, keys, item[1][...])
            keys.pop()

#______________________________________________________________________________

def nested_set(dic, keys, value):
    """
    This function sets a dict value using an arbitrary number of keys.
    _________________________
    Inputs:
        dic
            Dictionary to be modified
        keys
            List of keys defining nested dict path to set value
        value
            Value to set at position defined by keys

    Versioning:
        20140516 JRV - Created

    To do:
    _________________________
    """
    for key in keys[:-1]:
        dic = dic.setdefault(key, {})
    dic[keys[-1]] = value
    return dic

#______________________________________________________________________________



