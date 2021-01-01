"""
ReportSubversionRevisionPyMod.py

The function in this module will report the subversion revision number for the
path passed in as an argument.

2013-01-29 - JRV - Created
"""

from subprocess import Popen, PIPE
import os
import re


def ReportSubversionRevision(path):
    """ This function reports the SVN revision number.

    Inputs:
        path - Path to svn controlled directory or file (string)

    Outputs:
        curr_rev -              current SVN revsion number (int)
        last_changed rev -      SVN revision number of last change (int)

    2013-01-29 - JRV - Created
    2013-03-05 - JRV - Now returns URL string as well
    """

    if not os.path.isfile(path) and not os.path.isdir(path):
        print "ERROR: Path doesn't exist"
        return None

    p = Popen(['svn', 'info', path], stdout=PIPE)
    p.wait()
    out = p.communicate()[0]

    curr_rev = int(re.findall('\sRevision: (\d+)\s', out)[0])
    last_changed_rev = int(re.findall('\sLast Changed Rev: (\d+)\s', out)[0])
    url = re.findall('\sURL: (.+)\s', out)[0]

    return curr_rev, last_changed_rev, url
