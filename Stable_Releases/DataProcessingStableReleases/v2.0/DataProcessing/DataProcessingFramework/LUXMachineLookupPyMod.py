"""
LUXMachineLookupPyMod.py

Machine IDs (use the same as in LUXMachineLookup.m):
010 - James's laptop
009 - James's desktop
310 - on-site cluster

20121210 - CHF - Created
20121217 - JRV - Changed default to 0 and added my personal machines
                 * We should put the machine/ID accociations in xml *
                 Actually, they are in a table in the LUG. We should do a query (CHF)
"""

from socket import gethostname
from getpass import getuser

def LUXMachineLookup():
    """docstring for LUXMachineLookup"""

    # Need to check hostname and do a long list of case statements similar to MATLAB's LUXMachineLookup.m
    # Default for now is zero (unknown)
    machine_id = 0

    hostname = gethostname()
    username = getuser()

    if hostname == 'iolanthe' or \
            hostname == 'iolanthe.local':
        machine_id = 10  # James's Laptop
    elif hostname == 'pinafore':
        machine_id = 9  # James's Desktop

    if hostname == 'loc.sanfordlab.org':
        machine_id = 310 # on-site cluster - will need to add ORs for other nodes

    if hostname == 'luxanalysis.physics.yale.edu':
        machine_id = 401 # nightly build check machine at Yale
    
    if hostname == 'login001' or hostname == 'login002' or hostname.startswith('node') or hostname.startswith('smp'):
        if username == 'jverbus' or username == 'cah3':
            machine_id = 120
        else:
            machine_id = 121

    return machine_id
