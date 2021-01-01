"""
LUXUpdateLUGStatusPyMod.py

This Python Module updates the LUG CP record table with the \
number of files processed, as well as the running status.

20131122 CHF - Created
20131125 CHF - LUXUpdateLUGStatusPyMod is now imported

"""

import LUXDatabaseSQLPyMod
from str2numPyMod import str2num
import xml2dict
from ReportLUXcodePathPyMod import ReportDataProcessingPath
# On Oscar, need to run before:
# execfile('/users/jverbus/LUXcode/Trunk/DataProcessing/DataProcessingFramework/Utilities/LUXSetPythonPaths.py') 
# for these modules to be found

def LUXUpdateLUGStatus(cp, number_rq_files_finished, total_number_evt_files, status_flag):
    """
    Inputs:
        cp                       is the Complete Process Record identified
        number_rq_files_finished is an integer, number of files that have been processed
        total_number_evt_files   is an integer, total number of event files in dataset
        status_flag              is a int flag that means:
                                    status_flag = 1     Running
                                    status_flag = 2     Finished
    
    Outputs:
        status_flag              if it's 1, the LUG submission was successful!
    
    The code will insert a string such as 10/200 into the LUG, where in the example 10 is the number of
    rq files finished, and 200 is the total number of evt files. It will also insert the status_flag that will
    appear with the proper code ("Running" or "Finished") in the LUG display.
    
    This code needs to be run from a whitelisted computer or cluster (e.g. Oscar CCV) or from a computer on Sanford VPN
    in order for the query to work.
    
    20131122 CHF - Created
    """
    
    success_flag = None
    
    dp_path = ReportDataProcessingPath()
    xml_credentials = xml2dict.xml2dict(dp_path + 'DataProcessingFramework/Utilities/lug_query_credentials_insert.xml')
    
    # Define db object with database credentials from xml file.
    db = LUXDatabaseSQLPyMod.LUXDatabaseSQL(host=xml_credentials['credentials']['host'], \
        user=xml_credentials['credentials']['user'], \
        password=xml_credentials['credentials']['password'], \
        database=xml_credentials['credentials']['database'], \
        port=str2num(xml_credentials['credentials']['port']))
    
    # Connect to db
    db.connect()
    
    rq_files_finished_text = '%d/%d' % (number_rq_files_finished,total_number_evt_files)

    # Make the update
    query_str_update = "update control.lug_complete_process_record set files_processed_status \
    = '%s', status_flag = %d where cp = %d" % (rq_files_finished_text,status_flag,cp)
    result_update = db.query_db(query_str_update)

    # Now query values back to see if they got updated. If they did, exit with a success_flag = 1
    query_str_update_check = 'select files_processed_status, status_flag from control.lug_complete_process_record where cp = %d' % cp
    result_check = db.query_db(query_str_update_check)

    returned_number_rq_files_finished = result_check[0][0]
    returned_total_number_evt_files = result_check[0][1]
    
    success_flag = returned_number_rq_files_finished == rq_files_finished_text and int(returned_total_number_evt_files) == status_flag
    
    return success_flag
    
def LUXGetEB(cp):
    """
    Inputs:
        cp      is the Complete Process Record identified

    Outputs:
        eb      is the Event Builder Record identified

    The code return the Event Builder eb record for a particular Complete Process cp run.

    20131203 CHF - Created
    """

    eb = None

    dp_path = ReportDataProcessingPath()
    xml_credentials = xml2dict.xml2dict(dp_path + 'DataProcessingFramework/Utilities/lug_query_credentials_readonly.xml')

    # Define db object with database credentials from xml file.
    db = LUXDatabaseSQLPyMod.LUXDatabaseSQL(host=xml_credentials['credentials']['host'], \
        user=xml_credentials['credentials']['user'], \
        password=xml_credentials['credentials']['password'], \
        database=xml_credentials['credentials']['database'], \
        port=str2num(xml_credentials['credentials']['port']))

    # Connect to db
    db.connect()

    # Make the update
    query_str_update = "select eb from control.lug_complete_process_record where cp = %d" % cp
    result_update = db.query_db(query_str_update)
    
    if result_update:
        eb = int(result_update[0][0])
        
    return eb
        
        