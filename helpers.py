#!/usr/bin/env python
#-*- coding: utf-8 -*-

import os
import sys
from rdkit import Chem
from rdkit.Chem import AllChem

'''
this is the helper module that contains helping functions that relates to
both chemical and non-chemical perspectives
'''

########################################
#       NON-CHEMICAL RELATED
########################################

# set a unique scratch directory (specifically for ET-MEI cluster use)
_scratch_loc = '/scr/' + os.environ['USER'] + '/'
_scratch_dir = ''
def SetScratchDir(location=None):
    global _scratch_dir, _scratch_dir

    if location is not None:
        _scratch_dir = location
        return _scratch_dir

    #scratch already set
    if _scratch_dir != '':
        return _scratch_dir

    #parallel case
    try:
        from mpi4py import MPI
        myRank = MPI.COMM_WORLD.Get_rank()
    except ImportError:
        myRank = 0

    try:
        jobID = os.environ['JOB_ID']
    except KeyError:
        jobID = 'off_queue_job'

    _scratch_dir = _scratch_loc + str(jobID) + '_' + str(myRank) + '/'
    if not os.path.isdir(_scratch_dir):
        os.makedirs(_scratch_dir)

    print 'Scratch directory set to: ' + _scratch_dir
    return _scratch_dir

_warnScratch=False
def GetScratchDir():
    global _warnScratch
    
    if _scratch_dir == '' and not _warnScratch:
        print 'WARNING: scratch dir not set, current directory will be used.'
        _warnScratch = True

    return _scratch_dir


########################################
#       CHEMICAL RELATED
########################################

def Compute3DCoords(mol):
    #mol: rdkit RWMol or Mol
    molAddH = Chem.AddHs(mol)
    #calculate mol 3D coordinate
    AllChem.EmbedMolecule(molAddH, AllChem.ETKDG())

    molStr = Chem.MolToMolBlock(molAddH)
    #parse mol string file into list
    molCoords = []
    
    molStr2List = molStr.split('\n')
    #get number of atoms
    num = molStr2List[3].split()[0]

    for i in xrange(4,4+num):
        coords = molStr2List[i].split()[0:4]
        molCoords.append(coords)

    return molCoords
