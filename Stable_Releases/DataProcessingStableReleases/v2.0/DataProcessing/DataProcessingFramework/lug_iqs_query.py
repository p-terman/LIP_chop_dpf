#!/usr/bin/python

"""
Just read iqs from lug for cp dataset and make file with them
190321 PAT
"""

import os
import sys 

home_dir = os.path.expanduser("~")
execfile(home_dir + '/LUXCode/Trunk/DataProcessing/DataProcessingFramework/Utilities/LUXSetPythonPaths.py')\

# Using pymysql instead of MySQLdb
try:
    import pymysql
    pymysql.install_as_MySQLdb()
except ImportError:
    pass

#from LUXMachineLookupPyMod import LUXMachineLookup
import LUXDatabaseSQLPyMod
from filename_prefix2unix_timePyMod import filename_prefix2unix_time
from str2numPyMod import str2num
from ReportLUXcodePathPyMod import ReportLUXcodePath, ReportDataProcessingPath
import xml2dict
import itertools
from CompareRecursiveXmlDictPyMod import CompareRecursiveXmlDict
from ReportSubversionRevisionPyMod import ReportSubversionRevision
from bisect import bisect

cp = sys.argv[1]
print home_dir
# Define outputs
iq_xml_string = None

# Importing LUG query credentials from XML file.
#xml_credentials = xml2dict.xml2dict('/matlab/LUXCode/Trunk/DataProcessing/DataProcessingFramework/Utilities/lug_query_credentials_readonly.xml')
dp_path = ReportDataProcessingPath()
xml_credentials = xml2dict.xml2dict(home_dir + '/LUXCode/Trunk/DataProcessing/DataProcessingFramework/Utilities/lug_query_credentials_readonly.xml')

# Define db object with database credentials from xml file. This is read-only!
db = LUXDatabaseSQLPyMod.LUXDatabaseSQL(host=xml_credentials['credentials']['host'], \
    user=xml_credentials['credentials']['user'], \
    password=xml_credentials['credentials']['password'], \
    database=xml_credentials['credentials']['database'], \
    port=str2num(xml_credentials['credentials']['port']))
#print cp
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
query_str_cp = cp #'select * from control.lug_complete_process_record where cp =' + cp +';'
result_cp = db.query_db(query_str_cp)

# Make sure cp return matches and there is only one entry
if not result_cp:
    print '***ERROR: Query did not return any result for cp = %d' % cp
#    return iq_xml_string
elif result_cp[0][0] != cp:
    print '***ERROR: Record not found for cp = %d!' % cp
#    return iq_xml_string
elif len(result_cp) > 1:
    print '***ERROR: Too many matches for cp = %d!' % cp
#    return iq_xml_string


# *** Get the IQs xml strings
# Loop over list, then concatenate list into one string

# Turn iq list string into actual list
iq_str_raw = result_cp[0][6]
iq_list_str = iq_str_raw.split(',')
# In the case that no iq's are requested, return an empty iq xml.
#if(iq_list_str == [""]):
#  return data_processing_xml_string, "<iq></iq>"
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
        #return data_processing_xml_string, iq_xml_string
    elif result_iq[0][0] != iq_list[ii]:
        print '***ERROR: Record not found for iq = %d!' % iq_list[ii]
       # return data_processing_xml_string, iq_xml_string
    elif len(result_iq) > 1:
        print '***ERROR: Too many matches for iq = %d!' % iq_list[ii]
      #  return data_processing_xml_string, iq_xml_string
    iq_raw_xml_string_list.append(result_iq[0][2])

# Concatenate all the returned IQs into a single string
iq_xml_string = ''.join(iq_raw_xml_string_list)

iq_xml_path = home_dir + '/lug_iqs.xml' # smart_rq_path.get_final_rq_path() + '.lug_iqs_' + rand_str_generator() + '.xml'

f = open(iq_xml_path, 'w')
f.write(iq_xml_string)
f.close()
os.chmod(iq_xml_path, 0o775)
