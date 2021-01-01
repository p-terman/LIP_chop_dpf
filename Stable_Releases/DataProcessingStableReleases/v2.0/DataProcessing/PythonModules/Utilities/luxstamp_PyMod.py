# -*- coding: utf-8 -*-
"""
Converts filename_prefix and event_timestamp_samples to luxstamp
luxstamp defined in DAQ samples (10 ns) since Jan 1 2010 at 00:00, 0th sample.

Tue 20140514 CHF - Created

"""

import numpy as np
from datetime import datetime

def filename2luxstamp(filename_prefix,event_timestamp_samples=0):
    epoch = 'lux10_20100101T0000'
    fmat = 'lux10_%Y%m%dT%H%M'
    delta_min = datetime.strptime(filename_prefix, fmat) - datetime.strptime(epoch, fmat)
    luxstamp_samples = np.uint64(delta_min.total_seconds())*1e8
    return luxstamp_samples+np.uint64(event_timestamp_samples)