"""
filename_prefix2unix_timePyMod.py

If using in a MySQL query, just use:
   select unix_timestamp(str_to_date('lux10_20091023T1016',"lux10_%Y%m%dT%H%i")) as time
instead of this function! (use double %% to escape the actual character %)

20121214 CHF - Created
20130119 CHF - Realized that, with an internet connection, a MySQL query can do the conversion for you... :/

"""

from datetime import datetime
import calendar
from GetHourOffsetDSTPyMod import GetHourOffsetDST

def filename_prefix2unix_time(filename_prefix):
    """docstring for filename_prefix2unix_time"""
    
    unix_time = None
    
    sanford_tz_offset_hrs = -7
    fmt = '%Y%m%dT%H%M'
    date = datetime.strptime(filename_prefix[6:], fmt)
    local_unix_time = calendar.timegm(date.utctimetuple())
    
    dst_hrs = GetHourOffsetDST(date.year,date.month,date.day,date.hour)
    total_offset_hrs = sanford_tz_offset_hrs + dst_hrs
    print 'Time offset: %d hrs' % (total_offset_hrs)
    
    unix_time = local_unix_time - total_offset_hrs*3600
    
    return unix_time
    