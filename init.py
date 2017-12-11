#!/usr/bin/env python
#-*- coding: utf-8 -*-
import os
from rdkit import Chem 


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
    
    #case 1: restart file ('itX.lib.gz')

    #case 2: read from mprms.seedfile

    #case 3: read from mprms.seedlib
    
    #case 4: predefined as (benzene + cyclohexane)
    else:
        print 'Seeding library with presets'
        lib = []
        lib.append(Chem.MolFromSmiles('C1CCCCC1'))
        lib.append(Chem.MolFromSmiles('C1=CC=CC=C1'))

    return 


def StartPool(restart):
    
    #case 1: restart file ('pool.lib.gz')

    #case 2: read from mprms.poolfile

    #case 3: read from mprms.poollib

    #case 4: copy from library

    return
