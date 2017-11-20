#!/usr/bin/env python
#-*- coding: utf-8 -*-



######################################################
# Read run parameters from mprms.py,
# and store them all in StoredParams.py  
######################################################
def ReadMPRMS():
    return








######################################################
# Get starting pool and library
#
#  Library comes from, in order of precedence:
#      1. Restart files ('itX.lib.gz') (for restarts only)
#      2. Read from mprms.seedfile
#      3. Read from mprms.seedlib
#      4. Preset in previous options are not met
#         (benzene + cyclohexane)
#
#  Pool comes from:
#      1. Restart file ('pool.lib.gz') (for restarts only)
#      2. Read from mprms.poolfile
#      3. Read from mprms.poollib
#      4. Copy of library
######################################################
def StartLibrary(restart):
    return 

