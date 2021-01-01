"""
str2num.py

This function converts a string into a int or float.

2012-11-13 - JRV - Created
"""

import exceptions


def str2num(s):
    try:
        return int(s)
    except exceptions.ValueError:
        try:
            return float(s)
        except exceptions.ValueError:
            return None
