#!/usr/bin/python
"""
LUXSetPythonPaths.py

Use this script to set up your Python paths before running
th data processing framework code.

You will need to modify this for your DataProcessing path

2012-11-13 - JRV - Created
2013-01-24 - JRV - The default assumes your LUXcode is in your home directory.
"""

import sys
import os
from os.path import expanduser

home_path = expanduser("~")
dp_path = home_path + '/LUXcode/Stable_Releases/DataProcessingStableReleases/v2.0/DataProcessing/' 

# all paths need for python framework need to be added here
relative_paths = ['DataProcessingFramework',\
        '/',\
        'PythonModules/',
	'PythonModules/Utilities/']

# put all full paths into list
paths = [dp_path + s for s in relative_paths]

# add all paths to python path
for p in paths:
    for root, dirs, files in os.walk(p):
        if not os.path.basename(root).startswith('.'):
            if root not in sys.path:
                sys.path.append(root)
