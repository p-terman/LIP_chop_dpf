#PAT wrote based upon AL's various LUG programs 130320

# Set Python Paths
import sys
import os
home_dir = os.path.expanduser("~")
#from subprocess import Popen
#Popen(["source /users/pterman/pauls_world/bin/activate"])
execfile(home_dir + '/LUXCode/Trunk/DataProcessing/DataProcessingFramework/Utilities/LUXSetPythonPaths.py')

# Using pymysql instead of MySQLdb
try:
    import pymysql
    pymysql.install_as_MySQLdb()
except ImportError:
    pass

# Useful libraries and functions...
import LUXDatabaseSQLPyMod
import xml2dict
from str2numPyMod import str2num
import re

# Load read-only credentials to access the LUG
xml_credentials = xml2dict.xml2dict(home_dir + '/LUXCode/Trunk/DataProcessing/DataProcessingFramework/Utilities/lug_query_credentials_readonly.xml')

# Define db object with database credentials from xml file. This is read-only!
db = LUXDatabaseSQLPyMod.LUXDatabaseSQL(host=xml_credentials['credentials']['host'], \
                                            user=xml_credentials['credentials']['user'], \
                                            password=xml_credentials['credentials']['password'], \
                                            database=xml_credentials['credentials']['database'], \
                                            port=str2num(xml_credentials['credentials']['port']))

# Connect to db
db.connect()

# Query the LUG 
cp = sys.argv[1]

query_str_cp = 'select * from control.lug_complete_process_record where cp = %s;' % cp
result_cp = db.query_db(query_str_cp)

iq_str_raw = result_cp[0][6]
iq_list_str = iq_str_raw.split(',')

iq_list = map(int, iq_list_str)
iq_raw_xml_string_list = list()


# *** Get the IQs xml strings
# Loop over list, then concatenate list into one string

# Make IQ query
# Result values:
# [0][0] iq
# [0][1] filename_prefix
# [0][2] values_xml
# [0][3] comments
for ii in range(len(iq_list)):
    query_str_iq = 'select iq,filename_prefix,values_xml,comments from control.lug_iqs where iq = %d;' % iq_list[ii]
    result_iq = db.query_db(query_str_iq)
    iq_raw_xml_string_list.append(result_iq[0][2])

# Concatenate all the returned IQs into a single string
iq_xml_string = ''.join(iq_raw_xml_string_list)

iq_xml_path = home_dir + '/lug_iqs.xml' # smart_rq_path.get_final_rq_path() + '.lug_iqs_' + rand_str_generator() + '.xml'

f = open(iq_xml_path, 'w')
f.write(iq_xml_string)
f.close()
os.chmod(iq_xml_path, 0o775)
