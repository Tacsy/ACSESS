#!/usr/bin/env python
#-*- coding: utf-8 -*-
import os
from rdkit import Chem 
import mprms


######################################################
# Read run parameters from mprms.py,
# and store them all in StoredParams.py 
######################################################
def ReadMPRMS():
    import filters
    filters.FilterInit()
    return








######################################################
# Get starting library and pool
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
def StartLibAndPool(restart):
    startiter = 0

    if hasattr(mprms, 'nSeed'):
        nSeed = mprms.nSeed
    else:
        nSeed = -1
    
    #####################
    #   Start library   #
    #####################

    #case 1: restart file ('itX.lib.gz')
    if restart:
        pass
    #case 2: read from mprms.seedfile
    elif hasattr(mprms, 'seedFile'):
        seedFile = mprms.seedFile
        supplier = Chem.SmilesMolSupplier(seedFile)

        lib = [mol for mol in supplier]
        print "Seeding with " + str(len(lib)) + " molecules from " + seedFile

    #case 3: read from mprms.seedlib
    elif hasattr(mprms, 'seedLib'):
        lib = mprms.seedLib
        if nSeed > 0:
            lib = lib[:nSeed]
            lib = [Chem.MolFromSmiles(mol) for mol in lib]
        if len(lib) == 0:
            raise ValueError, "Seedlib is empty."
        print "Seeding with mprms.seedlib"

    #case 4: predefined as (benzene + cyclohexane)
    else:
        print 'Seeding library with presets'
        lib = []
        lib.append(Chem.MolFromSmiles('C1CCCCC1'))
        lib.append(Chem.MolFromSmiles('C1=CC=CC=C1'))
    
    #####################
    #     Start pool    #
    #####################
    
    #case 1: restart file ('pool.lib.gz')
    '''
    may work together with starting library
    '''
    #case 2: read from mprms.poolfile
    if hasattr(mprms, 'poolFile'):
        pool = [mol for mol in lib]
        for poolFile in mprms.poolFile:
            supplier = Chem.SmilesMolSupplier(poolFile)

            newpool = [mol for mol in supplier]
            pool += newpool

            print 'Initializing pool from ' + poolFile + ' of ' +\
                   str(len(newpool)) + ' molecules'

    #case 3: read from mprms.poollib
    elif hasattr(mprms, 'poolLib'):
        #poolLib is a list of SMILES
        pool = mprms.poolLib
        pool = [Chem.MolFromSmiles(mol) for mol in lib]
        pool = pool + lib
        if len(pool) == 0:
            raise ValueError, "Poollib is empty."
        print "Initializing pool from mprms.poolLib"
    
    #case 4: copy from library
    else:
        print 'Initializing pool from library'
        pool = [mol for mol in lib]

    return startiter, lib, pool
