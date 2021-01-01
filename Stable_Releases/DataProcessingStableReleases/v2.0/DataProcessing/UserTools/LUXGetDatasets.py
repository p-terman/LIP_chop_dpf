"""

Download all mat files from all datasets in an acquisition collection, from cp.

change variable destination_path to put your mat files. 
change name of acquisition collection in initial query_str.

2013-05-10 JJC - created.
"""

import LUXDatabaseSQLPyMod
import xml2dict
from str2numPyMod import str2num
import os
xml_credentials = xml2dict.xml2dict('lug_query_credentials_readonly.xml')


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

query_str = 'select filename_prefix from lug_acquisitions where collection = "BG S1+S2 - Prelim WIMP Run 3 Low Thresh";'
result = db.query_db(query_str)

destination_path = '/Volumes/Chapgate/'

for ii in range(len(result)):
    #print result[ii][0]
    query_str = 'select cp from lug_complete_process_record where filename_prefix = "%s";' % result[ii][0]
    result2 = db.query_db(query_str)
    if len(result2) > 0:
        #print result2[0][0]
        filename_prefix_cp = '%s_cp' % result[ii][0]
        filename_prefix_cp = filename_prefix_cp + '%05d' % result2[0][0]
        print filename_prefix_cp
		
        rsync_string = 'rsync -pruv lux_mirror_user@128.148.26.214:/Volumes/Unicorn/data/rq/%s/matfiles/ ' % filename_prefix_cp + '%s/' % destination_path + '%s/' % filename_prefix_cp
        print rsync_string
        try:
            os.system(rsync_string)
        except Exception:
            print 'did not find matfiles...'
            pass
        else:
            print 'downloading...'
