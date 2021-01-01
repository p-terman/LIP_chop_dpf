"""
LUXLUGQueriesPyMod.py

This file contains all LUG queries for the DP framework.
The queries themselves were originally writen by CHF,
but then compiled here.

        user=xml_credentials['credentials']['user'], \
20121218 JRV - Created
20130204 CHF - Reorganized
20130703 JRV - Added LUXGetElectronLifetime, this is a good template for calibration IQs
               Added LUXGetXYZCorrection, this gets all XYZ correction IQs
20130729 JRV - Added LUXGetXYCorIQ, this gets the iq_cor IQ for the Mercury edge fix
               TODO: Replace individual IQ queries with generic queries where possible.
                     example: LUXGetLatestIQ(filename_prefix, iq_type, algorithm, version)
20140625 CHF - Added 'group by filename_prefix' clause in all IQ queries to ensure duplicate entries
                     are collapsed into one
20140717 CHF - Added clause in LUXGetPMTHV to send email to designated DP person if PMT HV is not set in LUG for dataset
20140805 CHF - Changed LUXGetRequiredIQs. It used to use "if X in "iq_type_unique_list[ii]", 
                which is not unique for pmt_gains and veto_pmt_gains. Using == exact matching of iq_name explcitly
               Added elif iq_name == "pmt_veto_gains" case for veto PMT gains
               Created LUXGetVetoPMTGainCalibration function to return latest veto PMT gain calibration
20140822 CHF - Added function LUXCheckEventBuilderSettings. This checks whether there is an Extended Event Builder eb number
                to reuse. If not, it will create a new entry in the LUG.
20140822 CHF - Minor fix in checking the insert success in LUXCheckEventBuilderSettings
20140905 AL  - Added explicit use of VUV PMT gains (IQ 8506) for simulated datasets after the 5th of September 2014 (inclusive)
"""

from LUXMachineLookupPyMod import LUXMachineLookup
import LUXDatabaseSQLPyMod
from filename_prefix2unix_timePyMod import filename_prefix2unix_time
from str2numPyMod import str2num
from ReportLUXcodePathPyMod import ReportLUXcodePath, ReportDataProcessingPath
import xml2dict
import itertools
from CompareRecursiveXmlDictPyMod import CompareRecursiveXmlDict
from ReportSubversionRevisionPyMod import ReportSubversionRevision
from bisect import bisect
from time import time

def LUXPrepareDPSettings(filename_prefix, eb, data_processing_xml_string):
    """
    20121210 CHF - Created
    """

    cp = None

    # LUXGetRequiredIQs
    iq_list_string = LUXGetRequiredIQs(filename_prefix, data_processing_xml_string)
    #iq_list_string = ','.join(map(str, iq_list))
    # iq_list_string = '12,13,24'

    # *** Preparation
    # Check if gs are the same
    gs, gs_new_entry_flag = LUXCheckGlobalSettings(data_processing_xml_string)

    if gs_new_entry_flag == 1:
        print '*** LUXCheckGlobalSettings inserted the new entry gs = %d' % gs

    # Get machine ID
    machine_id = LUXMachineLookup()

    # *** Insert Complete Process Record - DONE
    cp = LUXInsertCompleteProcessRecord(filename_prefix, eb, gs, machine_id, iq_list_string)

    return cp

def LUXGetDPSettings(cp):
    """
    This function will return two xml strings:
    1) Data Processing xml
    2) IQs xml

    20121207 - CHF - Created
    """

    # Define outputs
    data_processing_xml_string = None
    iq_xml_string = None

    # Importing LUG query credentials from XML file.
    #xml_credentials = xml2dict.xml2dict('/matlab/LUXcode/Trunk/DataProcessing/DataProcessingFramework/Utilities/lug_query_credentials_readonly.xml')
    dp_path = ReportDataProcessingPath()
    xml_credentials = xml2dict.xml2dict(dp_path + 'DataProcessingFramework/Utilities/lug_query_credentials_readonly.xml')

    # Define db object with database credentials from xml file. This is read-only!
    db = LUXDatabaseSQLPyMod.LUXDatabaseSQL(host=xml_credentials['credentials']['host'], \
        user=xml_credentials['credentials']['user'], \
        password=xml_credentials['credentials']['password'], \
        database=xml_credentials['credentials']['database'], \
        port=str2num(xml_credentials['credentials']['port']))

    # Connect to db
    db.connect()

    # *** Get the Complete Process Settings
    # This simply gets all the values that match record cp
    # Result values:
    # [0][0] cp
    # [0][1] last_updated_date
    # [0][2] filename_prefix
    # [0][3] eb
    # [0][4] gs
    # [0][5] machine_id
    # [0][6] iq_id_list
    query_str_cp = 'select * from control.lug_complete_process_record where cp = %d;' % cp
    result_cp = db.query_db(query_str_cp)

    # Make sure cp return matches and there is only one entry
    if not result_cp:
        print '***ERROR: Query did not return any result for cp = %d' % cp
        return data_processing_xml_string, iq_xml_string
    elif result_cp[0][0] != cp:
        print '***ERROR: Record not found for cp = %d!' % cp
        return data_processing_xml_string, iq_xml_string
    elif len(result_cp) > 1:
        print '***ERROR: Too many matches for cp = %d!' % cp
        return data_processing_xml_string, iq_xml_string

    # *** Get the Global Settings, which we will use as Data Processing Settings
    # This simply gets all the values that match record gs
    # Result values:
    # [0][0] gs
    # [0][1] last_updated_date
    # [0][2] process_settings_xml
    # [0][3] comments
    query_str_gs = 'select * from control.lug_global_settings_record where gs = %d;' % result_cp[0][4]
    result_gs = db.query_db(query_str_gs)

    gs = result_gs[0][0]

    # Make sure gs return matches and there is only one entry
    if not result_gs:
        print '***ERROR: Query did not return any result for gs = %d' % gs
        return data_processing_xml_string, iq_xml_string
    elif result_gs[0][0] != gs:
        print '***ERROR: Record not found for gs = %d!' % gs
        return data_processing_xml_string, iq_xml_string
    elif len(result_gs) > 1:
        print '***ERROR: Too many matches for gs = %d!' % gs
        return data_processing_xml_string, iq_xml_string

    data_processing_xml_string = result_gs[0][2]

    # *** Get the IQs xml strings
    # Loop over list, then concatenate list into one string

    # Turn iq list string into actual list
    iq_str_raw = result_cp[0][6]
    iq_list_str = iq_str_raw.split(',')
    # In the case that no iq's are requested, return an empty iq xml.
    if(iq_list_str == [""]):
      return data_processing_xml_string, "<iq></iq>"
    iq_list = map(int, iq_list_str)
    iq_raw_xml_string_list = list()

    # Make IQ query
    # Result values:
    # [0][0] iq
    # [0][1] filename_prefix
    # [0][2] values_xml
    # [0][3] comments
    for ii in range(len(iq_list)):
        query_str_iq = 'select iq,filename_prefix,values_xml,comments from control.lug_iqs where iq = %d;' % iq_list[ii]
        result_iq = db.query_db(query_str_iq)

        if not result_iq:
            print '***ERROR: Query did not return any result for iq = %d' % iq_list[ii]
            return data_processing_xml_string, iq_xml_string
        elif result_iq[0][0] != iq_list[ii]:
            print '***ERROR: Record not found for iq = %d!' % iq_list[ii]
            return data_processing_xml_string, iq_xml_string
        elif len(result_iq) > 1:
            print '***ERROR: Too many matches for iq = %d!' % iq_list[ii]
            return data_processing_xml_string, iq_xml_string
        iq_raw_xml_string_list.append(result_iq[0][2])

    # Concatenate all the returned IQs into a single string
    iq_xml_string = ''.join(iq_raw_xml_string_list)

    # Return value, terminado
    return data_processing_xml_string, iq_xml_string


def LUXGetRequiredIQs(filename_prefix, data_processing_xml_string):
    """
    docstring for LUXGetRequiredIQs

    20121214 CHF - Created
    20130131 CHF - Minor fix, killed first entry in iq_type_unique_list since it's an empty by design
    20140805 CHF - Changed logic for if statement that select appropriate IQ. 
                    It used to use "if X in "iq_type_unique_list[ii]", which is not unique for pmt_gains and veto_pmt_gains
                    Logic was changed to use the actual name portion with exact == matching
    """

    iq_list_string = []

    # This list must match the case-by-case
    #known_iqs_list = ['lrfs', 'pmt_gains', 'pmt_afterpulsing', 'energy_calibration', 'electron_lifetime', 'detector_tilt']

    # Convert string into dictionary
    data_processing_xml_dict = xml2dict.xml2dict(data_processing_xml_string)

    # If there is only one module specified, xml2dict doesn't return a list of
    # module settings. Let's toss it in to a list of just itself.
    if(type(data_processing_xml_dict['data_processing_settings']['module']) != list):
              data_processing_xml_dict['data_processing_settings']['module'] = [data_processing_xml_dict['data_processing_settings']['module']]
    # Get list of required IQs
    iq_long_string = ''
    N = len(data_processing_xml_dict['data_processing_settings']['module'])
    for ii in range(N):
        iq_long_string += data_processing_xml_dict['data_processing_settings']['module'][ii]['required_iqs'] + ','
    iq_type_unique_list = list(set(iq_long_string.strip(",").split(",")))

    Q = len(iq_type_unique_list)
    print iq_type_unique_list
    iq_list = list()

    # *** This is where each IQ will be queried with a custom function
    for ii in range(Q):
    
        algorithm = ''
        version = ''
    
        split_iq = iq_type_unique_list[ii].split(':')
        iq_name = split_iq[0]

        if len(split_iq) == 2:
            algorithm = split_iq[1]
        elif len(split_iq) == 3:
            algorithm = split_iq[1]
            version = split_iq[2]
    
        if iq_type_unique_list[ii] != '':

            # Introduce direct IQ specification using a # escape. Allows for 
            # exact specification for which IQ one wants to use from within
            # the xml file.
            if iq_type_unique_list[ii].startswith("#"):
              try:
                iq_list.append(iq_type_unique_list[ii][1:])
              except:
                print " *** Warning: You've provided an improper IQ number ",
                print 'using #. The IQ "%s" will not be used.' % (iq_type_unique_list[ii])

            elif iq_name == 'pmt_gains':
                
                # Check if this is for a simulation
                if filename_prefix[0:5] == 'luxsm':
                    # TODO: This needs to be fixed to not have IQ numbers hardcoded for luxsm datasets- 20130926 JRV
                    filename_date_str = filename_prefix[6:14]
                    datetime_int = int(filename_date_str)
                    if datetime_int < 20130913: # Sept 13, 2013
                      # IQ = 47 is an IQ that has all PMTs at 16 mVns
                      iq_list.append('47')
                    elif datetime_int < 2040905: #Sept 05, 2014
                      # IQ = 517 is the golden insitu gains IQ.
                      iq_list.append('517')
		    else:
                      # IQ = 8506 is the VUV gains IQ.
                      iq_list.append('8506')

                    print " *** Warning: Running luxsm dataset, using hardcoded IQ number %d for PMT gains " % (int(iq_list[-1]))
                else:
                    if algorithm == "golden_insitu":
                        iq_pmt_gains = LUXGetInSituGains(filename_prefix, algorithm, version)
                    else:
                        # Get HV for that time
                        PMT_HV_list = LUXGetPMTHV(filename_prefix)
                        # Find best match for gains at that time
                        iq_pmt_gains = LUXGetPMTGainCalibration(filename_prefix, PMT_HV_list, algorithm, version)
                    # Add this iq number to the list
                    iq_list.append(iq_pmt_gains)

            elif iq_name == 'lrfs':                
                # Get iq number for LRF nearest to this time
                iq_lrfs = LUXGetLRFs(filename_prefix,algorithm,version)
                # Add this iq number to the list
                iq_list.append(iq_lrfs)

            elif iq_name == 'pmt_veto_gains':
                # pmt_veto_gains was chosen over 'veto_pmt_gains' since if statement logic is flawed. 
                # It used to use "if X in "iq_type_unique_list[ii]", which is not unique for pmt_gains and veto_pmt_gains
                # Logic was changed to use the actual name portion with exact == matching
                
                # Ideally, get the HV setttings for the veto PMTs for this dataset
                # Since veto PMT HV is not in SC, we can't really use this.
                # veto_PMT_HV_list = LUXGetVetoPMTHV(filename_prefix)
                veto_PMT_HV_list = ''

                # Use those HV settings for veto PMTs to find gain IQ that best matches these HV values
                # For now, this returns the latest calibration entry. When veto PMT HV is in SC, this should make a decision
                # of which is best IQ to get based on closest HV match
                iq_pmt_veto_gains = LUXGetVetoPMTGainCalibration(filename_prefix,veto_PMT_HV_list,algorithm,version)

                iq_list.append(iq_pmt_veto_gains)

            elif iq_name == 'pmt_afterpulsing':
                pass

            elif iq_name == 'energy_calibration':
                pass

            elif iq_name == 'xy_rec_cor':
                iq_list.append(LUXGetXYCorIQ(filename_prefix, algorithm, version))

            elif iq_name == 'electron_lifetime':
                iq_electron_lifetime = LUXGetElectronLifetime(filename_prefix, algorithm, version)
                iq_list += iq_electron_lifetime

            elif iq_name == 'z_dep_s1_correction':
                iq_list.append(LUXGetXYZCorrection(filename_prefix, 'z_dep_s1_correction', algorithm, version))

            elif iq_name == 's2_xy_correction':
                iq_list.append(LUXGetXYZCorrection(filename_prefix, 's2_xy_correction', algorithm, version))

            elif iq_name == 's1_xy_correction':
                iq_list.append(LUXGetXYZCorrection(filename_prefix, 's1_xy_correction', algorithm, version))

            elif iq_name == 's1_xyz_correction':
                iq_list.append(LUXGetXYZCorrection(filename_prefix, 's1_xyz_correction', algorithm, version))

            elif iq_name == 'detector_tilt':
                pass

            elif iq_name == 'single_electron':
                iq_1es2_size = LUXGetES2Size(filename_prefix, algorithm,version)
                if iq_1es2_size:
                    iq_list.append(iq_1es2_size)

            elif iq_name == 'per_pod_spurious_area':
                iq_per_pod_spurious_area = LUXGetSpuriousArea(algorithm,version)
                if iq_per_pod_spurious_area:
                    iq_list.append(iq_per_pod_spurious_area)

            else:
                print ' *** Warning: The requested IQ "%s" inside the data processing settings XML file is not supported!' % iq_type_unique_list[ii]


    iq_list_string = ', '.join(map(str, iq_list))

    return iq_list_string



def LUXGetPMTHV(filename_prefix):
    """
    docstring for LUXGetPMTHV

    20130308 CHF - Created
    20140717 CHF - Using 'if results[0][0]' instead of just 'if results' since this only returns True when accessing the actual entry
                   Adding 

    """

    PMT_HV_list = list()
    dp_path = ReportDataProcessingPath()
    xml_credentials = xml2dict.xml2dict(dp_path + 'DataProcessingFramework/Utilities/lug_query_credentials_readonly.xml')

    # Define db object with database credentials from xml file. This is read-only!
    db = LUXDatabaseSQLPyMod.LUXDatabaseSQL(host=xml_credentials['credentials']['host'], \
        user=xml_credentials['credentials']['user'], \
        password=xml_credentials['credentials']['password'], \
        database=xml_credentials['credentials']['database'], \
        port=str2num(xml_credentials['credentials']['port']))

    # Connect to db
    db.connect()

    query_HV = 'select pmt_hv_list_V from control.lug_acquisitions where filename_prefix = "%s" limit 1' % filename_prefix
    results = db.query_db(query_HV)

    if results[0][0]:
        PMT_HV_list = [float(x) for x in results[0][0].split(',')]
    else:
        # Adding this so the designated DP user gets an email if a dataset doesn't have the PMT HV set in the LUG.
        from SendEmail_PyMod import SendDPAlert
        dp_alert_user = 'alex.lindote@gmail.com'
        SendDPAlert(dp_alert_user,'The dataset %s failed because PMT HV was not manually set in LUG.' % filename_prefix)
    
    return PMT_HV_list


def LUXGetWienerPMTHV(filename_prefix):
    """
    docstring for LUXGetLeCroyPMTHV

    20121214 CHF - Created
    20130220 CHF - Added clause where, if PMT HV crate is off so there were no SC records during dataset for PMT HV,
                   the returned values are 0 V for each PMT
    20130308 CHF - Renamed from LUXGetPMTHV to LUXGetLeCroyPMTHV. This function should not be used with the DP framework for now. 
                   The HV will now be stored with the LUG acquisition entry (in ordered 1:122 fashion!), so please use LUXGetPMTHV instead.
    20130607 CHF - Adapted for use with Wiener PMT HV crate
    """

    PMT_HV_list = list()
    # Importing LUG query credentials from XML file. THIS NEEDS TO BE DYNAMICALLY CHANGED BASED ON PATH. Hardcoded for now
    #xml_credentials = xml2dict.xml2dict('/matlab/LUXcode/Trunk/DataProcessing/DataProcessingFramework/Utilities/lug_query_credentials_readonly.xml')
    dp_path = ReportDataProcessingPath()
    xml_credentials = xml2dict.xml2dict(dp_path + 'DataProcessingFramework/Utilities/lug_query_credentials_readonly.xml')

    # Define db object with database credentials from xml file. This is read-only!
    db = LUXDatabaseSQLPyMod.LUXDatabaseSQL(host=xml_credentials['credentials']['host'], \
        user=xml_credentials['credentials']['user'], \
        password=xml_credentials['credentials']['password'], \
        database=xml_credentials['credentials']['database'], \
        port=str2num(xml_credentials['credentials']['port']))

    # Connect to db
    db.connect()

    unix_time = filename_prefix2unix_time(filename_prefix)

    # *** Get the PMT HV information straight out of HV boards
    time_offset_mins = 2

    HV_slots = range(2)
    unix_time_lower = unix_time
    unix_time_upper = unix_time + (time_offset_mins * 60)

    raw_HV_results = list()
    for ii in HV_slots:
        query_HV = 'select Voltage from control.sc_sens_MPOD_HV_U%02d where time > %d and time < %d order by time desc limit 1' % (ii, unix_time_lower, unix_time_upper)
        temp = db.query_db(query_HV)
        if not temp:
            raw_HV_results.append([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        else:
            raw_HV_results.append(temp[0][0].split(','))

    chain_HV_results = itertools.chain(*raw_HV_results)
    PMT_HV_list_straight = list(chain_HV_results)

    # *** Get the indexing of HV boards to PMT channels

    # Loop per PMT to get closest analog chain entry BEFORE dataset was taken
    HV_inds = list()
    for ii in range(122):
        query_analog = ('select run_number, hv_supply_board, hv_supply_channel from lug_analog_chain \
        where pmt_channel = %d and entry_date <= %d order by entry_date desc limit 1;') % (ii + 1, unix_time)
        raw_analog = db.query_db(query_analog)      
        if not raw_analog:
            HV_inds.append(ii)
        else:
            HV_inds.append(raw_analog[0][1] * 16 + (raw_analog[0][2] + 1))

    # Reorganize PMT_HV_list_straight
    PMT_HV_list = list()
    for ii in range(122):
        this_ind = int(HV_inds[ii])
        PMT_HV_list.append(float(PMT_HV_list_straight[this_ind - 1]))

    if not PMT_HV_list:
        print('*** ERROR: Could not find PMT HV entry for this dataset')

    PMT_HV_list_int = [-1*int(x) for x in PMT_HV_list]

    return PMT_HV_list_int


def LUXGetLeCroyPMTHV(filename_prefix):
    """
    docstring for LUXGetLeCroyPMTHV

    20121214 CHF - Created
    20130220 CHF - Added clause where, if PMT HV crate is off so there were no SC records during dataset for PMT HV,
                   the returned values are 0 V for each PMT
    20130308 CHF - Renamed from LUXGetPMTHV to LUXGetLeCroyPMTHV. This function should not be used with the DP framework for now. 
                   The HV will now be stored with the LUG acquisition entry (in ordered 1:122 fashion!), so please use LUXGetPMTHV instead.
    """

    PMT_HV_list = list()
    # Importing LUG query credentials from XML file. THIS NEEDS TO BE DYNAMICALLY CHANGED BASED ON PATH. Hardcoded for now
    #xml_credentials = xml2dict.xml2dict('/matlab/LUXcode/Trunk/DataProcessing/DataProcessingFramework/Utilities/lug_query_credentials_readonly.xml')
    dp_path = ReportDataProcessingPath()
    xml_credentials = xml2dict.xml2dict(dp_path + 'DataProcessingFramework/Utilities/lug_query_credentials_readonly.xml')

    # Define db object with database credentials from xml file. This is read-only!
    db = LUXDatabaseSQLPyMod.LUXDatabaseSQL(host=xml_credentials['credentials']['host'], \
        user=xml_credentials['credentials']['user'], \
        password=xml_credentials['credentials']['password'], \
        database=xml_crexyntials['credentials']['database'], \
        port=str2num(xml_credentials['credentials']['port']))

    # Connect to db
    db.connect()

    unix_time = filename_prefix2unix_time(filename_prefix)
    
    # *** Get the PMT HV information straight out of HV boards
    time_offset_mins = 2

    HV_slots = range(12)
    unix_time_lower = unix_time
    unix_time_upper = unix_time + (time_offset_mins * 60)

    raw_HV_results = list()
    for ii in HV_slots:
        query_HV = 'select Voltage from control.sc_sens_PMT_HV_S%02d where time > %d and time < %d order by time desc limit 1' % (ii, unix_time_lower, unix_time_upper)
        temp = db.query_db(query_HV)
        if not temp:
            raw_HV_results.append([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        else:
            raw_HV_results.append(temp[0][0].split(','))

    chain_HV_results = itertools.chain(*raw_HV_results)
    PMT_HV_list_straight = list(chain_HV_results)

    # *** Get the indexing of HV boards to PMT channels

    # Loop per PMT to get closest analog chain entry BEFORE dataset was taken
    HV_inds = list()
    for ii in range(122):
        query_analog = ('select run_number, hv_supply_board, hv_supply_channel from lug_analog_chain \
        where pmt_channel = %d and entry_date <= %d order by entry_date desc limit 1;') % (ii + 1, unix_time)
        raw_analog = db.query_db(query_analog)      
        if not raw_analog:
            HV_inds.append(ii)
        else:
            HV_inds.append(raw_analog[0][1] * 12 + (raw_analog[0][2] + 1))

    # Reorganize PMT_HV_list_straight
    PMT_HV_list = list()
    for ii in range(122):
        this_ind = int(HV_inds[ii])
        PMT_HV_list.append(float(PMT_HV_list_straight[this_ind - 1]))

    if not PMT_HV_list:
        print('*** ERROR: Could not find PMT HV entry for this dataset')

    return PMT_HV_list

def LUXGetPMTGainCalibration(filename_prefix, PMT_HV_list, algorithm=None, version=None):
    """
    To do: Currently, there is no interpolation and submission back to the LUG (new entry).
    This must be done in the future.

    20121214 CHF - Created
    20130120 CHF - Now ignores iqs with strikeme = 1
    20130204 CHF - Added check to ignore crap entries with < 122 PMTs
                   Added check to see if any calibrations were found within tolerance_days. If not, end with empty return.
    20130813 JRV - Now compares to current time, not dataset time. This ensures a the golden gain IQ is used for all processes.
    20140625 CHF - Added 'group by filename_prefix' to queries
    """
    # Initialize values to return at the end
    iq_pmt_gains = None
    PMT_gains_mVns = None
    PMT_error_gains_mVns = None

    sec_per_day = 86400
    tolerance_days = 999999  # how many days back to look for a calibration
    print 'Looking for gain IQs as far back as %d days...' % tolerance_days
    tolerance_seconds = tolerance_days * sec_per_day

    #unix_time = filename_prefix2unix_time(filename_prefix)
    unix_time = int(time())  # compare to current time, not dataset time

    # Importing LUG query credentials from XML file.
    dp_path = ReportDataProcessingPath()
    xml_credentials_insert = xml2dict.xml2dict(dp_path + 'DataProcessingFramework/Utilities/lug_query_credentials_insert.xml')

    # Define db object with database credentials from xml file. This is read-only!
    db = LUXDatabaseSQLPyMod.LUXDatabaseSQL(host=xml_credentials_insert['credentials']['host'], \
        user=xml_credentials_insert['credentials']['user'], \
        password=xml_credentials_insert['credentials']['password'], \
        database=xml_credentials_insert['credentials']['database'], \
        port=str2num(xml_credentials_insert['credentials']['port']))
    # Connect to db
    db.connect()

    # Make query. Only select entries with the iq's filename_prefix (in unix time) within tolerance_seconds of this dataset
    #
    # iqs_raw[0][0] time
    # iqs_raw[0][1] iq
    # iqs_raw[0][2] values_xml
    #
    if algorithm and version:
        query_gains_all = ('select unix_timestamp(str_to_date(filename_prefix,"lux10_%%Y%%m%%dT%%H%%i")) as time, iq, values_xml from control.lug_iqs\
        where iq_type = "pmt_gains" and algorithm = "%s" and algorithm_version = %s and strikeme = 0 group by filename_prefix having abs(time - %d) < %d order by time desc;') % (algorithm, version, unix_time, tolerance_seconds)
    elif algorithm:
        query_gains_all = ('select unix_timestamp(str_to_date(filename_prefix,"lux10_%%Y%%m%%dT%%H%%i")) as time, iq, values_xml from control.lug_iqs\
        where iq_type = "pmt_gains" and algorithm = "%s" and strikeme = 0 group by filename_prefix having abs(time - %d) < %d order by time desc;') % (algorithm, unix_time, tolerance_seconds)
    else:
        query_gains_all = ('select unix_timestamp(str_to_date(filename_prefix,"lux10_%%Y%%m%%dT%%H%%i")) as time, iq, values_xml from control.lug_iqs\
        where iq_type = "pmt_gains" and strikeme = 0 group by filename_prefix having abs(time - %d) < %d order by time desc;') % (unix_time, tolerance_seconds)
    
    
    iqs_raw = db.query_db(query_gains_all)

    # Check if it found any entries. If not, return empty
    if not iqs_raw:
        err_str = '*** ERROR: Did not find any PMT gain calibrations in the last %d days.' % tolerance_days
        raise NameError(err_str)

    N = len(iqs_raw)
    HV_residuals = list()
    PMT_gains_mVns_all = [[] for x in xrange(N)]
    PMT_error_gains_mVns_all = [[] for x in xrange(N)]

    # Loop over all entries found, get the best match by computing HV residuals
    for ii in range(N):
        HV_dict = xml2dict.xml2dict(iqs_raw[ii][2])
        # If this was for only 20 PMTs, or some Run01 bullcrap, skip it and assign very large residual
        if len(HV_dict['iq']['fit']['channel']) != 122:
            HV_residuals.append(9999999999)
        else:
            this_residual = 0
            for pp in range(122):
                PMT_gains_mVns_all[ii].append(float(HV_dict['iq']['fit']['channel'][pp]['mVns_per_phe']))
                this_HV = round(float(HV_dict['iq']['fit']['channel'][pp]['PMT_bias_V']))
                # If this PMT was on... (otherwise don't take a penalty!)
                if abs(PMT_HV_list[pp]) > 200:
                    this_residual += (PMT_HV_list[pp] - this_HV) ** 2
            HV_residuals.append(this_residual ** 0.5)

    y = min(HV_residuals)
    # this picks the first instance. Since the MySQL query ordered them in time desc,
    # even if residuals are the same, the newest one will be picked by default
    ind = HV_residuals.index(y)

    # Get the best match
    PMT_gains_mVns = PMT_gains_mVns_all[ind]
    #    PMT_error_gains_mVns = PMT_error_gains_mVns_all[ind]
    # this is just a number!
    iq_pmt_gains = int(iqs_raw[ind][1])

    # Just return the iq number
    return iq_pmt_gains

def LUXGetVetoPMTGainCalibration(filename_prefix, PMT_HV_list, algorithm=None, version=None):
    """
    To do: One veto PMT HV is in SC, find the IQ with closest HV match

    20140805 CHF - Created based on LUXGetPMTGainCalibration

    """
    # Initialize values to return at the end
    iq_pmt_gains = None
    PMT_gains_mVns = None
    PMT_error_gains_mVns = None

    sec_per_day = 86400
    tolerance_days = 999999  # how many days back to look for a calibration
    print 'Looking for gain IQs as far back as %d days...' % tolerance_days
    tolerance_seconds = tolerance_days * sec_per_day

    #unix_time = filename_prefix2unix_time(filename_prefix)
    unix_time = int(time())  # compare to current time, not dataset time

    # Importing LUG query credentials from XML file.
    dp_path = ReportDataProcessingPath()
    xml_credentials_insert = xml2dict.xml2dict(dp_path + 'DataProcessingFramework/Utilities/lug_query_credentials_insert.xml')

    # Define db object with database credentials from xml file. This is read-only!
    db = LUXDatabaseSQLPyMod.LUXDatabaseSQL(host=xml_credentials_insert['credentials']['host'], \
        user=xml_credentials_insert['credentials']['user'], \
        password=xml_credentials_insert['credentials']['password'], \
        database=xml_credentials_insert['credentials']['database'], \
        port=str2num(xml_credentials_insert['credentials']['port']))
    # Connect to db
    db.connect()

    # Make query. Only select entries with the iq's filename_prefix (in unix time) within tolerance_seconds of this dataset
    #
    # iqs_raw[0][0] time
    # iqs_raw[0][1] iq
    # iqs_raw[0][2] values_xml
    #
    if algorithm and version:
        query_gains_all = ('select unix_timestamp(str_to_date(filename_prefix,"lux10_%%Y%%m%%dT%%H%%i")) as time, iq, values_xml from control.lug_iqs \
        where iq_type = "pmt_veto_gains" and algorithm = "%s" and algorithm_version = %s and strikeme = 0 group by filename_prefix having abs(time - %d) < %d order by time desc;') % (algorithm, version, unix_time, tolerance_seconds)
    elif algorithm:
        query_gains_all = ('select unix_timestamp(str_to_date(filename_prefix,"lux10_%%Y%%m%%dT%%H%%i")) as time, iq, values_xml from control.lug_iqs \
        where iq_type = "pmt_veto_gains" and algorithm = "%s" and strikeme = 0 group by filename_prefix having abs(time - %d) < %d order by time desc;') % (algorithm, unix_time, tolerance_seconds)
    else:
        query_gains_all = ('select unix_timestamp(str_to_date(filename_prefix,"lux10_%%Y%%m%%dT%%H%%i")) as time, iq, values_xml from control.lug_iqs \
        where iq_type = "pmt_veto_gains" and strikeme = 0 group by filename_prefix having abs(time - %d) < %d order by time desc;') % (unix_time, tolerance_seconds)

    iqs_raw = db.query_db(query_gains_all)

    # Check if it found any entries. If not, return empty
    if not iqs_raw:
        err_str = '*** ERROR: Did not find any PMT gain calibrations in the last %d days.' % tolerance_days
        raise NameError(err_str)

    # *** IMPORTANT *** For now, return the latest calibration in time
    iq_pmt_gains = int(iqs_raw[0][1])

    # Just return the iq number
    return iq_pmt_gains

def LUXGetInSituGains(filename_prefix, algorithm, version):
    """ Documentation for LUXGetMeasuredGains.

    This function will query for the "measured_pmt_gains" IQ. This IQ is an
    averaged PMT gain measured with the in situ gain measurement software.

    WHO MADE THIS? EDIT HISTORY? (CHF)

    20140625 CHF - Added 'group by filename_prefix' to queries

    """

    # Initialize values to return at the end
    iq_pmt_gains = None
    PMT_gains_mVns = None
    PMT_error_gains_mVns = None

    sec_per_day = 86400
    tolerance_days = 9999  # how many days back to look for a calibration
    tolerance_seconds = tolerance_days * sec_per_day
    print 'Looking for in situ gain IQs as far back as %d days...' % tolerance_days

    unix_time = filename_prefix2unix_time(filename_prefix)

    # Importing LUG query credentials from XML file.
    dp_path = ReportDataProcessingPath()
    xml_credentials_insert = xml2dict.xml2dict(dp_path + 'DataProcessingFramework/Utilities/lug_query_credentials_insert.xml')

    # Define db object with database credentials from xml file. This is read-only!
    db = LUXDatabaseSQLPyMod.LUXDatabaseSQL(host=xml_credentials_insert['credentials']['host'], \
        user=xml_credentials_insert['credentials']['user'], \
        password=xml_credentials_insert['credentials']['password'], \
        database=xml_credentials_insert['credentials']['database'], \
        port=str2num(xml_credentials_insert['credentials']['port']))
    # Connect to db
    db.connect()

    # Make query. Only select entries with the iq's filename_prefix (in unix time) within tolerance_seconds of this dataset
    #
    # iqs_raw[0][0] time
    # iqs_raw[0][1] iq
    # iqs_raw[0][2] values_xml
    #
    if algorithm and version:
        query_gains_all = ('select unix_timestamp(str_to_date(filename_prefix,"lux10_%%Y%%m%%dT%%H%%i")) as time, iq, values_xml from control.lug_iqs\
        where iq_type = "pmt_gains" and algorithm = "%s" and algorithm_version = %s and strikeme = 0 group by filename_prefix having abs(time - %d) < %d order by time desc;') % (algorithm, version, unix_time, tolerance_seconds)
    elif algorithm:
        query_gains_all = ('select unix_timestamp(str_to_date(filename_prefix,"lux10_%%Y%%m%%dT%%H%%i")) as time, iq, values_xml from control.lug_iqs\
        where iq_type = "pmt_gains" and algorithm = "%s" and strikeme = 0 group by filename_prefix having abs(time - %d) < %d order by time desc;') % (algorithm, unix_time, tolerance_seconds)
    else:
        print "*** Warning: There were insufficient parameters specified to attempt importing in situ gains."
        print "*** Warning: This will undoubtedly cause errors shortly..."
    
    
    iqs_raw = db.query_db(query_gains_all)

    # Check if it found any entries. If not, return empty
    if not iqs_raw:
        err_str = '*** ERROR: Did not find any PMT gain calibrations in the last %d days.' % tolerance_days
        raise NameError(err_str)

    # Look for the pmt_gains:golden_insitu IQ that is closest in time to the dataset you are processing.
    nearest_iq = {"index":-1, "delta_t":1e50} # iq_index of closest matching IQ and the time difference between them.
    for result in iqs_raw:
        time, iq_index, iq = result
        if abs(time - unix_time) < nearest_iq["delta_t"]:
            nearest_iq["index"] = iq_index
            nearest_iq["delta_t"] = abs(time - unix_time)
    # Just return the iq number
    return nearest_iq["index"]



def LUXGetLRFs(filename_prefix,algorithm,version):
    """
    LUXGetLRFs will return the closest LRF iq value in time to filename_prefix
    (within X number of tolerance days, currently set to many, many years)
    It will ignore any iq entries with strikeme = 1

    20130120 CHF - Created
    20130204 CHF - Added a check to see if any calibrations were found within tolerance_days. If not, end with empty return
    20130814 JRV - Now uses current time for comparison, not dataset time; also orders by iq not time
    20140625 CHF - Added 'group by filename_prefix' to queries

    """
    # Initialize values to return at the end
    iq_lrfs = None

    sec_per_day = 86400
    tolerance_days = 9999999  # how many days back to look for a calibration
    print 'Looking for LRF IQs as far back as %d days...' % tolerance_days
    tolerance_seconds = tolerance_days * sec_per_day

    #unix_time = filename_prefix2unix_time(filename_prefix)
    unix_time = int(time())  # compare to current time, not dataset time

    # Importing LUG query credentials from XML file.
    dp_path = ReportDataProcessingPath()
    xml_credentials_insert = xml2dict.xml2dict(dp_path + 'DataProcessingFramework/Utilities/lug_query_credentials_insert.xml')

    # Define db object with database credentials from xml file. This is read-only!
    db = LUXDatabaseSQLPyMod.LUXDatabaseSQL(host=xml_credentials_insert['credentials']['host'], \
        user=xml_credentials_insert['credentials']['user'], \
        password=xml_credentials_insert['credentials']['password'], \
        database=xml_credentials_insert['credentials']['database'], \
        port=str2num(xml_credentials_insert['credentials']['port']))
    # Connect to db
    db.connect()

    # Make query. Only select entries with the iq's filename_prefix (in unix time) within tolerance_seconds of this dataset
    
    if algorithm and version:
        query_gains_all = ('select unix_timestamp(str_to_date(filename_prefix,"lux10_%%Y%%m%%dT%%H%%i")) as time, iq, values_xml from control.lug_iqs\
        where iq_type = "lrfs" and algorithm = "%s" and algorithm_version = %s and strikeme = 0 group by filename_prefix having abs(time - %d) < %d order by iq desc;') % (algorithm, version, unix_time, tolerance_seconds)
    elif algorithm:
        query_gains_all = ('select unix_timestamp(str_to_date(filename_prefix,"lux10_%%Y%%m%%dT%%H%%i")) as time, iq, values_xml from control.lug_iqs\
        where iq_type = "lrfs" and algorithm = "%s" and strikeme = 0 group by filename_prefix having abs(time - %d) < %d order by iq desc;') % (algorithm, unix_time, tolerance_seconds)
    else:
        query_gains_all = ('select unix_timestamp(str_to_date(filename_prefix,"lux10_%%Y%%m%%dT%%H%%i")) as time, iq, values_xml from control.lug_iqs\
        where iq_type = "lrfs" and strikeme = 0 group by filename_prefix having abs(time - %d) < %d order by iq desc;') % (unix_time, tolerance_seconds)
    
    iqs_raw = db.query_db(query_gains_all)
    # Check if it found any entries. If not, return empty
    if not iqs_raw:
        err_str = '*** ERROR: Did not find any LRF IQs in the last %d days.' % tolerance_days
        raise NameError(err_str)

    # Pick the latest one, which is entry 0 since we order it by descending iq
    iq_lrfs = int(iqs_raw[0][1])

    # Just return the iq number
    return iq_lrfs


def LUXGetXYCorIQ(filename_prefix,algorithm,version):
    """
    LUXGetXYCorIQ will return the closest xy_cor iq value in time to filename_prefix
    (within X number of tolerance days, currently set to many, many years)
    It will ignore any iq entries with strikeme = 1

    20130729 JRV - Created
    20140625 CHF - Added 'group by filename_prefix' to queries

    """
    # Initialize values to return at the end
    iq = None

    sec_per_day = 86400
    tolerance_days = 9999999  # how many days back to look for a calibration
    print 'Looking for xy_rec_cor IQs as far back as %d days...' % tolerance_days
    tolerance_seconds = tolerance_days * sec_per_day

    unix_time = filename_prefix2unix_time(filename_prefix)

    # Importing LUG query credentials from XML file.
    dp_path = ReportDataProcessingPath()
    xml_credentials_insert = xml2dict.xml2dict(dp_path + 'DataProcessingFramework/Utilities/lug_query_credentials_insert.xml')

    # Define db object with database credentials from xml file. This is read-only!
    db = LUXDatabaseSQLPyMod.LUXDatabaseSQL(host=xml_credentials_insert['credentials']['host'], \
        user=xml_credentials_insert['credentials']['user'], \
        password=xml_credentials_insert['credentials']['password'], \
        database=xml_credentials_insert['credentials']['database'], \
        port=str2num(xml_credentials_insert['credentials']['port']))
    # Connect to db
    db.connect()

    # Make query. Only select entries with the iq's filename_prefix (in unix time) within tolerance_seconds of this dataset
    
    if algorithm and version:
        query_gains_all = ('select unix_timestamp(str_to_date(filename_prefix,"lux10_%%Y%%m%%dT%%H%%i")) as time, iq, values_xml from control.lug_iqs\
        where iq_type = "xy_rec_cor" and algorithm = "%s" and algorithm_version = %s and strikeme = 0 group by filename_prefix having abs(time - %d) < %d order by time desc;') % (algorithm, version, unix_time, tolerance_seconds)
    elif algorithm:
        query_gains_all = ('select unix_timestamp(str_to_date(filename_prefix,"lux10_%%Y%%m%%dT%%H%%i")) as time, iq, values_xml from control.lug_iqs\
        where iq_type = "xy_rec_cor" and algorithm = "%s" and strikeme = 0 group by filename_prefix having abs(time - %d) < %d order by time desc;') % (algorithm, unix_time, tolerance_seconds)
    else:
        query_gains_all = ('select unix_timestamp(str_to_date(filename_prefix,"lux10_%%Y%%m%%dT%%H%%i")) as time, iq, values_xml from control.lug_iqs\
        where iq_type = "xy_rec_cor" and strikeme = 0 group by filename_prefix having abs(time - %d) < %d order by time desc;') % (unix_time, tolerance_seconds)
    
    iqs_raw = db.query_db(query_gains_all)
    # Check if it found any entries. If not, return empty
    if not iqs_raw:
        err_str = '*** ERROR: Did not find any iq_rec_cor IQs in the last %d days.' % tolerance_days
        raise NameError(err_str)

    # Pick the latest one, which is entry 0 since we order it by descending time
    iq = int(iqs_raw[0][1])

    # Just return the iq number
    return iq

def LUXGetES2Size(filename_prefix,algorithm,version):
    """
    LUXGetES2Size will return the 1eS2 size iq value for a given filename_prefix
    It will ignore any iq entries with strikeme = 1, and limits to 1

    20130907 PHP - Created
    20140625 CHF - Didn't change anything, but noticed queries are flawed. Author should fix.
    """
    # Initialize values to return at the end
    iq = None

    print 'Looking for single_electron IQs for %s...' % filename_prefix

    unix_time = filename_prefix2unix_time(filename_prefix)

    # Importing LUG query credentials from XML file.
    dp_path = ReportDataProcessingPath()
    xml_credentials_insert = xml2dict.xml2dict(dp_path + 'DataProcessingFramework/Utilities/lug_query_credentials_insert.xml')

    # Define db object with database credentials from xml file. This is read-only!
    db = LUXDatabaseSQLPyMod.LUXDatabaseSQL(host=xml_credentials_insert['credentials']['host'], \
        user=xml_credentials_insert['credentials']['user'], \
        password=xml_credentials_insert['credentials']['password'], \
        database=xml_credentials_insert['credentials']['database'], \
        port=str2num(xml_credentials_insert['credentials']['port']))
    # Connect to db
    db.connect()

    # Make query. Only select entries with the iq's filename_prefix (in unix time) within tolerance_seconds of this dataset
    
    if algorithm and version:
        query_1es2_all = ('select unix_timestamp(str_to_date(filename_prefix,"lux10_%%Y%%m%%dT%%H%%i")) as time, iq, values_xml from control.lug_iqs\
        where filename_prefix = "%s" and iq_type = "single_electron" and algorithm = "%s" and algorithm_version = %s and strikeme = 0 limit 1;') % (filename_prefix, algorithm, version)
    elif algorithm:
        query_1es2_all = ('select unix_timestamp(str_to_date(filename_prefix,"lux10_%%Y%%m%%dT%%H%%i")) as time, iq, values_xml from control.lug_iqs\
        where filename_prefix = "%s" and iq_type = "single_electron" and algorithm = "%s" and strikeme = 0 order by algorithm_version desc limit 1;') % (filename_prefix, algorithm)
    else:
        query_1es2_all = ('select unix_timestamp(str_to_date(filename_prefix,"lux10_%%Y%%m%%dT%%H%%i")) as time, iq, values_xml from control.lug_iqs\
        where filename_prefix = "%s" and iq_type = "single_electron" and strikeme = 0 order by iq desc limit 1;') % (filename_prefix)
        
    iqs_raw = db.query_db(query_1es2_all)
    
    # Check if it found any entries. If not, return empty
    if not iqs_raw:
        err_str = '*** ERROR: Did not find any 1eS2 size IQs for dataset %s' % filename_prefix
        print err_str
        iq = None
    else:
        # Pick the latest one, which is entry 0 since we order it by descending time
        iq = int(iqs_raw[0][1])
        # Just return the iq number

    return iq


def LUXGetXYZCorrection(filename_prefix, iq_type, algorithm, version):
    """
    LUXGetXYZCorrection will return the closest XYZ iq value in time to filename_prefix
    (within X number of tolerance days, currently set to many, many years)
    It will ignore any iq entries with strikeme = 1

    This function is used to grab the following IQs using the iq_type input:
    z_dep_s1_correction
    s2_xy_correction
    s1_xy_correction
    s1_xyz_correction

    20130703 - JRV - Created
    20140625 - CHF - Added 'group by filename_prefix' to queries
    """
    # Initialize values to return at the end
    iq_xyz_correction = None

    sec_per_day = 86400
    tolerance_days = 9999999  # how many days back to look for a calibration
    print 'Looking for %s IQs as far back as %d days...' % (iq_type, tolerance_days)
    tolerance_seconds = tolerance_days * sec_per_day

    unix_time = filename_prefix2unix_time(filename_prefix)

    # Importing LUG query credentials from XML file.
    dp_path = ReportDataProcessingPath()
    xml_credentials_insert = xml2dict.xml2dict(dp_path + 'DataProcessingFramework/Utilities/lug_query_credentials_insert.xml')

    # Define db object with database credentials from xml file. This is read-only!
    db = LUXDatabaseSQLPyMod.LUXDatabaseSQL(host=xml_credentials_insert['credentials']['host'], \
        user=xml_credentials_insert['credentials']['user'], \
        password=xml_credentials_insert['credentials']['password'], \
        database=xml_credentials_insert['credentials']['database'], \
        port=str2num(xml_credentials_insert['credentials']['port']))
    # Connect to db
    db.connect()

    # Make query. Only select entries with the iq's filename_prefix (in unix time) within tolerance_seconds of this dataset
    
    if algorithm and version:
        query_gains_all = ('select unix_timestamp(str_to_date(filename_prefix,"lux10_%%Y%%m%%dT%%H%%i")) as time, iq, values_xml from control.lug_iqs\
        where iq_type = "%s" and algorithm = "%s" and algorithm_version = %s and strikeme = 0 group by filename_prefix having abs(time - %d) < %d order by time desc;') % (iq_type, algorithm, version, unix_time, tolerance_seconds)
    elif algorithm:
        query_gains_all = ('select unix_timestamp(str_to_date(filename_prefix,"lux10_%%Y%%m%%dT%%H%%i")) as time, iq, values_xml from control.lug_iqs\
        where iq_type = "%s" and algorithm = "%s" and strikeme = 0 group by filename_prefix having abs(time - %d) < %d order by time desc;') % (iq_type, algorithm, unix_time, tolerance_seconds)
    else:
        query_gains_all = ('select unix_timestamp(str_to_date(filename_prefix,"lux10_%%Y%%m%%dT%%H%%i")) as time, iq, values_xml from control.lug_iqs\
        where iq_type = "%s" and strikeme = 0 group by filename_prefix having abs(time - %d) < %d order by time desc;') % (iq_type, unix_time, tolerance_seconds)
    
    iqs_raw = db.query_db(query_gains_all)
    # Check if it found any entries. If not, return empty
    if not iqs_raw:
        err_str = '*** ERROR: Did not find any %s IQs in the last %d days.' % (iq_type, tolerance_days)
        raise NameError(err_str)

    date_list = [iq[0] for iq in iqs_raw]
    closest_iq = date_list.index(min(date_list, key=lambda x:abs(x - unix_time)))

    iq_xyz_correction = iqs_raw[closest_iq][1]

    # Just return the iq number
    return iq_xyz_correction


def LUXGetElectronLifetime(filename_prefix,algorithm,version):
    """
    LUXGetElectronLifetime will return the closest e-lifetime iq value in time to filename_prefix
    (within X number of tolerance days, currently set to many, many years)
    It will ignore any iq entries with strikeme = 1

    20130703 - JRV - Created from CHFs lrfs iq function
    20140625 - CHF - Added 'group by filename_prefix' to queries
    """
    # Initialize values to return at the end
    iq_electron_lifetime = []

    sec_per_day = 86400
    tolerance_days = 9999999  # how many days back to look for a calibration
    print 'Looking for electron_lifetime IQs as far back as %d days...' % tolerance_days
    tolerance_seconds = tolerance_days * sec_per_day

    unix_time = filename_prefix2unix_time(filename_prefix)

    # Importing LUG query credentials from XML file.
    dp_path = ReportDataProcessingPath()
    xml_credentials_insert = xml2dict.xml2dict(dp_path + 'DataProcessingFramework/Utilities/lug_query_credentials_insert.xml')

    # Define db object with database credentials from xml file. This is read-only!
    db = LUXDatabaseSQLPyMod.LUXDatabaseSQL(host=xml_credentials_insert['credentials']['host'], \
        user=xml_credentials_insert['credentials']['user'], \
        password=xml_credentials_insert['credentials']['password'], \
        database=xml_credentials_insert['credentials']['database'], \
        port=str2num(xml_credentials_insert['credentials']['port']))
    # Connect to db
    db.connect()

    # Make query. Only select entries with the iq's filename_prefix (in unix time) within tolerance_seconds of this dataset
    
    if algorithm and version:
        query_all = ('select unix_timestamp(str_to_date(filename_prefix,"lux10_%%Y%%m%%dT%%H%%i")) as time, iq, values_xml from control.lug_iqs\
        where iq_type = "electron_lifetime" and algorithm = "%s" and algorithm_version = %s and strikeme = 0 group by filename_prefix having abs(time - %d) < %d order by time asc;') % (algorithm, version, unix_time, tolerance_seconds)
    elif algorithm:
        query_all = ('select unix_timestamp(str_to_date(filename_prefix,"lux10_%%Y%%m%%dT%%H%%i")) as time, iq, values_xml from control.lug_iqs\
        where iq_type = "electron_lifetime" and algorithm = "%s" and strikeme = 0 group by filename_prefix having abs(time - %d) < %d order by time asc;') % (algorithm, unix_time, tolerance_seconds)
    else:
        query_all = ('select unix_timestamp(str_to_date(filename_prefix,"lux10_%%Y%%m%%dT%%H%%i")) as time, iq, values_xml from control.lug_iqs\
        where iq_type = "electron_lifetime" and strikeme = 0 group by filename_prefix having abs(time - %d) < %d order by time asc;') % (unix_time, tolerance_seconds)
    
    iqs_raw = db.query_db(query_all)
    # Check if it found any entries. If not, return empty
    if not iqs_raw:
        err_str = '*** ERROR: Did not find any electron_lifetime IQs in the last %d days.' % tolerance_days
        raise NameError(err_str)

    # Pick the latest one, which is entry 0 since we order it by descending time
    #iq_lrfs = int(iqs_raw[0][1])

    # get list of iq dates and find iqs before and after target date
    date_list = [iq[0] for iq in iqs_raw]
    ind_grt_target = bisect(date_list, unix_time)
    ind_less_target = ind_grt_target - 1

    if len(date_list) > ind_grt_target:
        iq_electron_lifetime.append(iqs_raw[ind_grt_target][1])  # if cal exists after target date, append it
    iq_electron_lifetime.append(iqs_raw[ind_less_target][1])  # append cal exists before target date


    # Just return the iq number
    return iq_electron_lifetime

def LUXGetSpuriousArea(algorithm, version):
    """
    LUXGetSpuriousArea will return the iq of specified algorithm and version
    It will ignore any iq entries with strikeme = 1

    20140616 - AC
    20140625 - CHF - Added 'group by filename_prefix' to queries
    """ 
    # Initialize values to return at the end
    iq_spurious_area = []
    
    sec_per_day = 86400
    tolerance_days = 9999999  # how many days back to look for a calibration
    print 'Looking for per_pod_spurious_area IQs as far back as %d days...' % tolerance_days
    tolerance_seconds = tolerance_days * sec_per_day

    #unix_time = filename_prefix2unix_time(filename_prefix)
    unix_time = int(time())  # compare to current time, not dataset time

    # Importing LUG query credentials from XML file.
    dp_path = ReportDataProcessingPath()
    xml_credentials_insert = xml2dict.xml2dict(dp_path + 'DataProcessingFramework/Utilities/lug_query_credentials_insert.xml')

    # Define db object with database credentials from xml file. This is read-only!
    db = LUXDatabaseSQLPyMod.LUXDatabaseSQL(host=xml_credentials_insert['credentials']['host'], \
        user=xml_credentials_insert['credentials']['user'], \
        password=xml_credentials_insert['credentials']['password'], \
        database=xml_credentials_insert['credentials']['database'], \
        port=str2num(xml_credentials_insert['credentials']['port']))
    # Connect to db
    db.connect()

    # Make query. Only select entries with the iq's filename_prefix (in unix time) within tolerance_seconds of this dataset
    
    if algorithm and version:
        query_all = ('select unix_timestamp(str_to_date(filename_prefix,"lux10_%%Y%%m%%dT%%H%%i")) as time, iq, values_xml from control.lug_iqs\
        where iq_type = "per_pod_spurious_area" and algorithm = "%s" and algorithm_version = %s and strikeme = 0 group by filename_prefix having abs(time - %d) < %d order by iq desc;') % (algorithm, version, unix_time, tolerance_seconds)
    elif algorithm:
        query_all = ('select unix_timestamp(str_to_date(filename_prefix,"lux10_%%Y%%m%%dT%%H%%i")) as time, iq, values_xml from control.lug_iqs\
        where iq_type = "per_pod_spurious_area" and algorithm = "%s" and strikeme = 0 group by filename_prefix having abs(time - %d) < %d order by iq desc;') % (algorithm, unix_time, tolerance_seconds)
    else:
        query_all = ('select unix_timestamp(str_to_date(filename_prefix,"lux10_%%Y%%m%%dT%%H%%i")) as time, iq, values_xml from control.lug_iqs\
        where iq_type = "per_pod_spurious_area" and strikeme = 0 group by filename_prefix having abs(time - %d) < %d order by iq desc;') % (unix_time, tolerance_seconds)
    
    iqs_raw = db.query_db(query_all)
    # Check if it found any entries. If not, return empty
    if not iqs_raw:
        err_str = '*** ERROR: Did not find any per_pod_spurious_area IQs in the last %d days.' % tolerance_days
        raise NameError(err_str)

    # Pick the latest one, which is entry 0 since we order it by descending iq
    iq_spurious_area = int(iqs_raw[0][1])    

    return iq_spurious_area
    

def LUXCheckGlobalSettings(data_processing_xml_string):
    """
    LUXCheckGlobalSettings will go through each entry in the LUG global settings table
    and find the latest identical match to data_processing_xml settings.
    If it finds a match, it will return the gs value (and the gs_new_entry_flag = 0).
    If it does not find a match, it will make a new entry and return the new gs value
    (with gs_new_entry_flag = 1).

    The logic for finding the identical match is in CompareRecursiveXmlDictPyMod

    20121207 CHF - Created
    20130120 CHF - Added functionality for reusing gs settings
    20130131 CHF - Added comparison checks for svn_location, svn_last_changed_rev
    20130305 JRV - ???
    """

    # Initialize
    gs = None
    gs_new_entry_flag = None
    data_processing_xml_dict = xml2dict.xml2dict(data_processing_xml_string)
    recycle_gs_flag = 0
    svn_location = None
    svn_last_changed_rev = None

    # Importing LUG query credentials from XML file.
    dp_path = ReportDataProcessingPath()
    xml_credentials_insert = xml2dict.xml2dict(dp_path + 'DataProcessingFramework/Utilities/lug_query_credentials_insert.xml')

    # Figure out what SVN folder we are using (Trunk or Stable_Releases)
    svn_rev, svn_last_changed_rev, svn_url = ReportSubversionRevision(dp_path)

    if svn_url.find('Trunk') > -1:
        svn_location = 'Trunk'
    elif svn_url.find('Stable_Releases') > -1:
        svn_location = 'Stable_Releases'
    else:
        print('*** ERROR: Could not determine if framework is running on Trunk or Stable_Releases!')
        return gs, gs_new_entry_flag

    print 'Running DP framework from LUXcode/%s, last changed revision %d' % (svn_location, svn_last_changed_rev)

    # Define db object with database credentials from xml file. This is read-only!
    db = LUXDatabaseSQLPyMod.LUXDatabaseSQL(host=xml_credentials_insert['credentials']['host'], \
        user=xml_credentials_insert['credentials']['user'], \
        password=xml_credentials_insert['credentials']['password'], \
        database=xml_credentials_insert['credentials']['database'], \
        port=str2num(xml_credentials_insert['credentials']['port']))

    # Connect to db
    db.connect()

    # *** Get ALL the global settings entries
    # [0][0] gs
    # [0][1] process_settings_xml
    # [0][2] svn_location
    # [0][3] svn_last_changed_rev
    # Entries returned are in descending gs order, so [0][0] should be latest gs entry
    # Only considering entries that match the svn location and last changed revision
    query_gs = "select gs,process_settings_xml,svn_location,svn_last_changed_rev from control.lug_global_settings_record where svn_location = '%s' and svn_last_changed_rev = %d order by gs desc;" % (svn_location,svn_last_changed_rev)
    results_gs = db.query_db(query_gs)

    # For every entry in the global settings table, check if it's identical or not
    # Since entries are in descending gs order, when you find the first identical entry (the entries might not all be unique...)
    # then just pick that one. This guarantees you use the latest identical entry.
    for qq in range(len(results_gs)):
        gs_xml_dict = xml2dict.xml2dict(results_gs[qq][1])
        different_flag = CompareRecursiveXmlDict(data_processing_xml_dict,gs_xml_dict)
        if different_flag == 0:
            recycle_gs_flag = 1;
            gs = int(results_gs[qq][0])
            print '*** Found a match! gs = %d' % gs
            break

    # *** Make a decision: does the record exist already, or should you make a new one?
    if recycle_gs_flag == 1:
        # Done, reuse these settings
        gs_new_entry_flag = 0
        return gs, gs_new_entry_flag
    else:
        # *** INSERT NEW ENTRY
        # Query to get latest entry
        previous_gs = db.query_db('select gs from control.lug_global_settings_record order by gs desc limit 1;')

        # Insert this entry
        insert_str_begin = 'insert into control.lug_global_settings_record (last_updated_date,process_settings_xml,svn_location,svn_last_changed_rev) '
        insert_str_end = 'values(unix_timestamp(current_timestamp),"%s","%s",%d)' % (data_processing_xml_string,svn_location,svn_last_changed_rev)
        insert_str = insert_str_begin + insert_str_end

        result = db.query_db(insert_str)  # 'result' is not very useful, returns empty either case

        # Query to get latest entry
        this_gs = db.query_db('select gs from control.lug_global_settings_record order by gs desc limit 1;')

        # Do a basic check to see if the record was inserted
        # This checks that after insertion, the primary key cp increased by at least 1
        if int(this_gs[0][0]) - int(previous_gs[0][0]) <= 0:
            print '***ERROR: Query was not inserted successfully'
            return gs, gs_new_entry_flag

        # Ask back for the newest record that has the settings we just inserted, see if that matches this_cp[0][0]
        check_str = 'select gs from control.lug_global_settings_record where process_settings_xml = "%s" \
                and svn_location = "%s" and svn_last_changed_rev = %d order by gs desc limit 1;' % (data_processing_xml_string,svn_location,svn_last_changed_rev)

        gs_check = db.query_db(check_str)

        # Check that the newest record we just inserted matches the one returned
        if int(gs_check[0][0]) != int(this_gs[0][0]):
            print '***ERROR: A record was inserted, but there was a problem with the record values (value mismatch between insert and check query)'
            return gs, gs_new_entry_flag

        # If you got this far, all is well
        gs_new_entry_flag = 1
        gs = int(this_gs[0][0])

    return gs, gs_new_entry_flag


def LUXInsertCompleteProcessRecord(filename_prefix, eb, gs, machine_id, iq_list_string):
    """
    This function inserts the Complete Process Record settings into the LUG.
    It does a couple of sanity checks to ensure the record was inserted successfully.

    20121207 CHF - Created
    """

    cp = None

    # Importing LUG query credentials from XML file.
    dp_path = ReportDataProcessingPath()
    xml_credentials_insert = xml2dict.xml2dict(dp_path + 'DataProcessingFramework/Utilities/lug_query_credentials_insert.xml')

    # Define db object with database credentials from xml file. This is read-only!
    db = LUXDatabaseSQLPyMod.LUXDatabaseSQL(host=xml_credentials_insert['credentials']['host'], \
        user=xml_credentials_insert['credentials']['user'], \
        password=xml_credentials_insert['credentials']['password'], \
        database=xml_credentials_insert['credentials']['database'], \
        port=str2num(xml_credentials_insert['credentials']['port']))

    # Connect to db
    db.connect()

    # Query to get latest entry
    previous_cp = db.query_db('select cp from control.lug_complete_process_record order by cp desc limit 1;')

    # Insert this entry
    insert_str_begin = 'insert into control.lug_complete_process_record (last_updated_date,filename_prefix,eb,gs,machine_id,iq_list) '
    insert_str_end = 'values(unix_timestamp(current_timestamp), "%s", %d, %d, %d, "%s")' % (filename_prefix, eb, gs, machine_id, iq_list_string)
    insert_str = insert_str_begin + insert_str_end

    result = db.query_db(insert_str)  # 'result' is not very useful, returns empty either case

    # Query to get latest entry
    this_cp = db.query_db('select cp from control.lug_complete_process_record order by cp desc limit 1;')

    # Do a basic check to see if the record was inserted
    # This checks that after insertion, the primary key cp increased by at least 1
    if int(this_cp[0][0]) - int(previous_cp[0][0]) <= 0:
        print '***ERROR: Query was not inserted successfully'
        return cp

    # Ask back for the newest record that has the settings we just inserted, see if that matches this_cp[0][0]
    check_str_begin = 'select cp from control.lug_complete_process_record where '
    check_str_end = 'filename_prefix="%s" and eb=%d and gs=%d and machine_id=%d and iq_list="%s" order by cp desc limit 1;' % (filename_prefix, eb, gs, machine_id, iq_list_string)
    check_str = check_str_begin + check_str_end

    cp_check = db.query_db(check_str)

    # Check that the newest record we just inserted matches the one returned
    if int(cp_check[0][0]) != int(this_cp[0][0]):
        print '***ERROR: A record was inserted, but there was a problem with the record values (value mismatch between insert and check query)'
        return cp

    cp = int(this_cp[0][0])

    return cp


def LUXCheckEventBuilderSettings(eb_xml_string):
    """
    LUXCheckEventBuilderSettings will go through each entry in the LUG global settings table
    and find the latest identical match to data_processing_xml settings.
    If it finds a match, it will return the gs value (and the eb_new_entry_flag = 0).
    If it does not find a match, it will make a new entry and return the new eb value
    (with eb_new_entry_flag = 1).

    The logic for finding the identical match is in CompareRecursiveXmlDictPyMod

    20140821 CHF - Created
    """

    # Initialize
    eb = None
    eb_new_entry_flag = None
    eb_xml_dict = xml2dict.xml2dict(eb_xml_string)
    recycle_eb_flag = 0

    # Get some data out of the eb_xml_dict
    pretrigger_samples = int(eb_xml_dict['data_processing_settings']['module']['parameters']['trigger']['pretrigger'])
    posttrigger_samples = int(eb_xml_dict['data_processing_settings']['module']['parameters']['trigger']['posttrigger'])
    trigger_method = eb_xml_dict['data_processing_settings']['module']['parameters']['trigger']['method']
    salt_evt_flag = int(eb_xml_dict['data_processing_settings']['module']['parameters']['salt_evt'])

    livetime_veto_flag = eb_xml_dict['data_processing_settings']['module']['parameters']['livetime_veto']['method']
    if (livetime_veto_flag == 'NULL') or (not livetime_veto_flag) or (livetime_veto_flag == 'None'):
        livetime_veto_flag = 0
    else:
        livetime_veto_flag = 1

    preselection_flag = eb_xml_dict['data_processing_settings']['module']['parameters']['preselection']['method']
    if (preselection_flag == 'NULL') or (not preselection_flag) or (preselection_flag == 'None'):
        preselection_flag = 0
    else:
        preselection_flag = 1

    salt_shake_id = eb_xml_dict['data_processing_settings']['module']['parameters']['salt_shake_id']
    if (salt_shake_id == 'NULL') or (not salt_shake_id) or (salt_shake_id == 'None'):
        salt_shake_id = 0
    else:
        salt_shake_id = int(salt_shake_id)


    # Importing LUG query credentials from XML file.
    dp_path = ReportDataProcessingPath()
    xml_credentials_insert = xml2dict.xml2dict(dp_path + 'DataProcessingFramework/Utilities/lug_query_credentials_insert.xml')

    # Get EEB version
    svn_rev, svn_last_changed_rev, svn_url = ReportSubversionRevision(dp_path + '/CppModules/ExtendedEventBuilder/EventBuilder.cpp')

    # Define db object with database credentials from xml file. This is read-only!
    db = LUXDatabaseSQLPyMod.LUXDatabaseSQL(host=xml_credentials_insert['credentials']['host'], \
        user=xml_credentials_insert['credentials']['user'], \
        password=xml_credentials_insert['credentials']['password'], \
        database=xml_credentials_insert['credentials']['database'], \
        port=str2num(xml_credentials_insert['credentials']['port']))

    # Connect to db
    db.connect()

    # *** Get ALL the global settings entries
    # [0][0] eb
    # [0][1] event_builder_version
    # [0][2] settings_xml
    # Entries returned are in descending gs order, so [0][0] should be latest gs entry
    # Only considering entries that match the pretrigger and postrigger settings
    query_eb = 'select eb,event_builder_version,settings_xml from control.lug_event_builder_record where \
        svn_last_changed_rev = %d \
        and pretrigger_samples = %d and posttrigger_samples = %d \
        and salt_evt_flag = %d and salt_shake_id = %d \
        and trigger_method = "%s" and livetime_veto_flag = %d \
        and preselection_flag = %d \
        order by eb desc;' % (svn_last_changed_rev,pretrigger_samples,posttrigger_samples,salt_evt_flag,salt_shake_id,trigger_method,livetime_veto_flag,preselection_flag)
    results_eb = db.query_db(query_eb)

    # For every entry in the eb settings table, check if it's identical or not
    # Since entries are in descending eb order, when you find the first identical entry (the entries might not all be unique...)
    # then just pick that one. This guarantees you use the latest identical entry.
    if results_eb:
        for qq in range(len(query_eb)):
            checking_eb_xml_dict = xml2dict.xml2dict(results_eb[qq][2])
            different_flag = CompareRecursiveXmlDict(checking_eb_xml_dict,eb_xml_dict)
            if different_flag == 0:
                recycle_eb_flag = 1;
                eb = int(results_eb[qq][0])
                print '*** Found a match! eb = %d' % eb
                break

    # *** Make a decision: does the record exist already, or should you make a new one?
    if recycle_eb_flag:
        # Done, reuse these settings
        eb_new_entry_flag = 0
        return eb, eb_new_entry_flag
    else:
        # *** INSERT NEW ENTRY
        # Query to get latest entry
        previous_eb = db.query_db('select eb from control.lug_event_builder_record order by eb desc limit 1;')

        # Insert this entry
        insert_str_begin = 'insert into control.lug_event_builder_record (last_updated_date,svn_last_changed_rev,pretrigger_samples,posttrigger_samples,\
            salt_evt_flag,salt_shake_id,trigger_method,livetime_veto_flag,preselection_flag,settings_xml) '

        insert_str_end = 'values(unix_timestamp(current_timestamp),%d,%d,%d,%d,%d,"%s",%d,%d,"%s");' % \
            (svn_last_changed_rev,pretrigger_samples,posttrigger_samples,\
            salt_evt_flag,salt_shake_id,trigger_method,livetime_veto_flag,preselection_flag,eb_xml_string)

        insert_str = insert_str_begin + insert_str_end

        result = db.query_db(insert_str)  # 'result' is not very useful, returns empty either case

        # Query to get latest entry
        this_eb = db.query_db('select eb from control.lug_event_builder_record order by eb desc limit 1;')

        # Do a basic check to see if the record was inserted
        # This checks that after insertion, the primary key cp increased by at least 1
        if int(this_eb[0][0]) - int(previous_eb[0][0]) <= 0:
            print '***ERROR: Query was not inserted successfully'
            return eb, eb_new_entry_flag

        # Ask back for the newest record that has the settings we just inserted, see if that matches this_eb[0][0]
        check_str = 'select eb from control.lug_event_builder_record where settings_xml = "%s" order by eb desc;' % (eb_xml_string)

        eb_check = db.query_db(check_str)

        # Check that the newest record we just inserted matches the one returned
        if int(eb_check[0][0]) != int(this_eb[0][0]):
            print '***ERROR: A record was inserted, but there was a problem with the record values (value mismatch between insert and check query)'
            return eb, eb_new_entry_flag

        # If you got this far, all is well
        eb_new_entry_flag = 1
        eb = int(this_eb[0][0])

    return eb, eb_new_entry_flag

