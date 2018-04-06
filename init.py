#!/usr/bin/env python
#-*- coding: utf-8 -*-
import os, sys
sys.path.append('.')
sys.path.append('./Filters/')
from rdkit import Chem 
from rdkithelpers import *
import mprms
import importlib
import random


######################################################
# Read run parameters from mprms.py,
# and store them all in StoredParams.py 
######################################################
def ReadMPRMS():
    primitiveTypes=(str, float, bool, int, list, tuple, type(None))
    normalvar = lambda var:type(var) in primitiveTypes
    notbuiltin= lambda var:not var.startswith('_')
    def goodvar(var, mod):return notbuiltin(var) and normalvar(getattr(mod, var))

    def SetModVal(name,ModAssign):
        if hasattr(mprms,name):
            if not restart: print>> myopts,name+'='+str(getattr(mprms,name))
            setattr(ModAssign,name,
                    getattr(mprms,name) )
    def SetMPRMDefault(name,value):
        if not hasattr(mprms,name):
            setattr(mprms,name,value)
        else:
            if not restart: print >> myopts,name+'='+str(getattr(mprms,name))

    if hasattr(mprms, 'rseed'):
        random.seed(mprms.rseed)
    else:
        seed = random.randint(1,1000)
        print "random seed:", seed
        random.seed(seed)
    if not hasattr(mprms, 'optimize'):
        mprms.optimize=False

    mprms.MxAtm=50
    mprms.EdgeRatio=0.1
    mprms.EdgeLen=10

    # Decide which modules have to be imported:
    _modules = ['mutate', 
                'filters', 'Filters.DefaultFilters',
                'drivers',
                'rdkithelpers']
    if hasattr(mprms, 'metric'):
        _modules.append('distance')
    else:
        _modules.append('similarity')
    if mprms.optimize:
        _modules.append('objective')

    # Import the modules / Set the global variables / Initiate modules
    for module in _modules:
        # get normal global variables of the modules
        _Mod = importlib.import_module(module)
        modvars = [ var for var in dir(_Mod) if goodvar(var, _Mod)]
        for modvar in modvars:
            # check if mprms has same attr:
            if hasattr(mprms, modvar):
                #than change it:
                print "set attr: {}.{}".format(module, modvar)
                setattr(_Mod, modvar, getattr(mprms, modvar))

        # initialize module if it has a Init function
        if hasattr(_Mod, 'Init') and callable(getattr(_Mod, 'Init')):
            getattr(_Mod, 'Init')()
    
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
        for i in reversed(xrange(mprms.nGen+1)):
            filename= "it{}.smi".format(i)
            if os.path.isfile(filename):
                startiter=i+1
                break
        else: raise SystemExit('RESTART but no iteration files found')
        supplier = Chem.SmilesMolSupplier(filename, sanitize=False)
        lib = [mol for mol in supplier]
        print 'RESTARTING calculation from iteration', startiter-1,
        print 'with {} molecules'.format(len(lib))
    #case 2: read from mprms.seedfile
    elif hasattr(mprms, 'seedFile'):
        seedFile = mprms.seedFile
        supplier = Chem.SmilesMolSupplier(seedFile, sanitize=False)

        lib = [mol for mol in supplier]
        print "Seeding with " + str(len(lib)) + " molecules from " + seedFile

    #case 3: read from mprms.seedlib
    elif hasattr(mprms, 'seedLib'):
        lib = mprms.seedLib
        if nSeed > 0:
            lib = lib[:nSeed]
            lib = [Chem.MolFromSmiles(mol, sanitize=False) for mol in lib]
        if len(lib) == 0:
            raise ValueError, "Seedlib is empty."
        print "Seeding with mprms.seedlib"

    #case 4: predefined as (benzene + cyclohexane)
    else:
        print 'Seeding library with presets'
        lib = []
        lib.append(Chem.MolFromSmiles('C1CCCCC1', sanitize=False))
        lib.append(Chem.MolFromSmiles('C1=CC=CC=C1', sanitize=False))
    
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
        supplier = Chem.SmilesMolSupplier(mprms.poolFile, sanitize=False)
        newpool = [mol for mol in supplier]
        pool += newpool
        print 'Initializing pool from ' + mprms.poolFile + ' of ' +\
                   str(len(newpool)) + ' molecules'

    #case 3: read from mprms.poollib
    elif hasattr(mprms, 'poolLib'):
        #poolLib is a list of SMILES
        pool = mprms.poolLib
        pool = [Chem.MolFromSmiles(mol, False) for mol in lib]
        pool = pool + lib
        if len(pool) == 0:
            raise ValueError, "Poollib is empty."
        print "Initializing pool from mprms.poolLib"
    
    #case 4: copy from library
    else:
        print 'Initializing pool from library'
        pool = [mol for mol in lib]

    # set isosmi
    setisosmi = lambda mol:mol.SetProp('isosmi', Chem.MolToSmiles(mol, True))
    map(setisosmi, pool)
    map(setisosmi, lib)
    ######### sanitize lib & pool
    ll1,lp1= (len(lib), len(pool))
    lib    = filter(Sane, lib)
    pool   = filter(Sane, pool)
    ll2,lp2= (len(lib), len(pool))
    if not ll1==ll2:print "removed {:d} unsane molecules from lib  | libsize : {:d}".format(ll1-ll2, ll2)
    if not lp1==lp2:print "removed {:d} unsane molecules from pool | poolsize: {:d}".format(lp1-lp2, lp2)

    return startiter, lib, pool
