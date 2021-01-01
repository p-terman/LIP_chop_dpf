"""
ReportLUXcodePathPyMod.py

The functions in this module will report the path to LUXcode
and the path to LUXcode/Trunk/DataProcessing

2012-12-17 - JRV - Created
"""

import os


def ReportLUXcodePath():
    """ returns the absolute path to LUXcode """

    LUXcode_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..')) + '/'

    return LUXcode_path


def ReportDataProcessingPath():
    """ returns the absolute path to LUXcode/Trunk/DataProcessing/ """

    data_processing_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')) + '/'

    return data_processing_path
