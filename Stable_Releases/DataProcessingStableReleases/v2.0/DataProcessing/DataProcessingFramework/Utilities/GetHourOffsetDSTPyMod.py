"""
GetHourOffsetDSTPyMod.py



20121214 CHF - Created

"""

import calendar

def GetHourOffsetDST(year,month,day,hour):
    """docstring for GetHourOffsetDST"""
    
    dst = None
    
    dow = calendar.weekday(year,month,day) + 1 # Monday = 1, Sunday = 7
    previous_sunday = day - dow
    
    # *** To DST, or not to DST...
    
    # definitely not DST
    if month < 3 or month > 11:
        dst = 0
    # definitely DST
    elif month > 3 and month < 11:
        dst = 1
    # case for March: DST is on second Sunday of March after 2 am
    elif month == 3:
        # If the previous Sunday was 8th or more, this is the third week of March,
        # which means that DST already started
        if previous_sunday >= 8:
            dst = 1
        # If this IS the second Sunday of the month AND it's after 2 am
        elif dow == 7 and day >= 8 and hour > 2:
            dst = 1
        # If none of the above, then it's before DST transition
        else:
            dst = 0
    # case for November: DST ends first Sunday of November
    elif month == 11:
        # If the previous Sunday was the 1st or after, this is the second week of November,
        # which means that DST already finished
        if previous_sunday >= 1:
            dst = 0
        # IF this IS the first Sunday of the month and it's after 2 am
        elif dow == 7 and day <= 7 and hour > 2:
            dst = 0
        else:
            dst = 1
    
    return dst
    