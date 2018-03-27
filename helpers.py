#!/usr/bin/env python
#-*- coding: utf-8 -*-

import os
import sys
import cPickle as pickle
import gzip
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

'''
this is the helper module that contains helping functions that relates to
both chemical and non-chemical perspectives
'''
########################################
#       PRE-DEFINED VARIABLES
########################################

atomNum = {'Ac': 89, 'Ag': 47, 'Al': 13, 'Am': 95, 'Ar': 18, 'As': 33, 'At': 85, 'Au': 79,
        'B': 5, 'Ba': 56, 'Be': 4, 'Bh': 107, 'Bi': 83, 'Bk': 97, 'Br': 35, 'C': 6,
        'Ca': 20, 'Cd': 48, 'Ce': 58, 'Cf': 98, 'Cl': 17, 'Cm': 96, 'Cn': 112, 'Co': 27,
        'Cr': 24, 'Cs': 55, 'Cu': 29, 'Db': 105, 'Ds': 110, 'Dy': 66, 'Er': 68, 'Es': 99,
        'Eu': 63, 'F': 9, 'Fe': 26, 'Fm': 100, 'Fr': 87, 'Ga': 31, 'Gd': 64, 'Ge': 32,
        'H': 1, 'He': 2, 'Hf': 72, 'Hg': 80, 'Ho': 67, 'Hs': 108, 'I': 53, 'In': 49, 'Ir': 77,
        'K': 19, 'Kr': 36, 'La': 57, 'Li': 3, 'Lr': 103, 'Lu': 71, 'Md': 101, 'Mg': 12,
        'Mn': 25, 'Mo': 42, 'Mt': 109, 'N': 7, 'Na': 11, 'Nb': 41, 'Nd': 60, 'Ne': 10,
        'Ni': 28, 'No': 102, 'Np': 93, 'O': 8, 'Os': 76, 'P': 15, 'Pa': 91, 'Pb': 82,
        'Pd': 46, 'Pm': 61, 'Po': 84, 'Pr': 59, 'Pt': 78, 'Pu': 94, 'Ra': 88, 'Rb': 37,
        'Re': 75, 'Rf': 104, 'Rg': 111, 'Rh': 45, 'Rn': 86, 'Ru': 44, 'S': 16, 'Sb': 51,
        'Sc': 21, 'Se': 34, 'Sg': 106, 'Si': 14, 'Sm': 62, 'Sn': 50, 'Sr': 38, 'Ta': 73,
        'Tb': 65, 'Tc': 43, 'Te': 52, 'Th': 90, 'Ti': 22, 'Tl': 81, 'Tm': 69, 'U': 92,
        'Uuh': 116, 'Uuo': 118, 'Uup': 115, 'Uuq': 114, 'Uus': 117, 'Uut': 113, 'V': 23,
        'W': 74, 'Xe': 54, 'Y': 39, 'Yb': 70, 'Zn': 30, 'Zr': 40}
elements = {}
for element, num in atomNum.iteritems(): 
    elements[num] = element

########################################
#       NON-CHEMICAL RELATED
########################################

# read a zip file directly into memory, then iterate through it.
# much faster than using gzip.open, which is broken
def FastGZ(fname):
    import cStringIO
    io_method = cStringIO.StringIO
    import subprocess
    p = subprocess.Popen(["zcat", fname], stdout = subprocess.PIPE)
    fh = io_method(p.communicate()[0])
    assert p.returncode == 0
    
    return fh

def Depickle(fname, gz=False):
    if fname.split('.')[-1] == 'gz':
        gz = True
    if gz:
        myfile = FastGZ(fname)
    else:
        myfile = open(fname, 'r')
    tmp = pickle.load(myfile)
    myfile.close()

    return tmp

def Enpickle(obj, fname, protocal=2, gz=False):
    if fname.split('.')[-1] == 'gz':
        gz = True
    if gz:
        myfile = gzip.open(fname, 'wb')
    else:
        myfile = open(fname, 'wb')
    pickle.dump(obj, myfile, protocal)
    myfile.close()

# get SMILES name of specific atom
def GetSmilesName(atom):
    global elements
    symbol = elements[atom.GetAtomicNum()]
    if atom.GetIsAromatic():
        symbol = symbol.lower()
    else:
        symbol = symbol.upper()
    
    return symbol

# extract number in a string
def ExtractNum(string):
    outstr = ''
    for char in string:
        if char.isdigit():
            outstr += charg
    
    return int(outstr)

# get the file format
def GetFileFormat(filename):
    fields = filename.split('.')
    if fields[-1].lower() == 'gz':
        return fields[-2].lower()
    else:
        return fields[-1].lower()

def GetBaseName(filename, stripstr=True):
    if stripstr == True:
        return filename.split('/')[-1].split('.')[0]
    elif stripstr:
        return filename.split('/')[-1].rstrip(stripstr)
    else:
        return name.split('/')[-1]

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

def Normalize(vec):
    n = np.sqrt(np.dot(vec, vec))
    return vec/n

########################################
#       CHEMICAL RELATED
########################################

# calculate 3D coordinate of a molecule using RDKit build-in functions
def Compute3DCoords(mol):
    #mol: rdkit RWMol or Mol
    molAddH = Chem.AddHs(mol)
    #calculate mol 3D coordinate
    AllChem.EmbedMolecule(molAddH, AllChem.ETKDG())

    molStr = Chem.MolToMolBlock(molAddH)
    print molStr

    #parse mol string file into list
    molCoords = []
    
    molStr2List = molStr.split('\n')

    print "molStr2List:", molStr2List
    #get number of atoms
    num = int(molStr2List[3].split()[0])

    for i in xrange(4,4+num):
        coords = molStr2List[i].split()[0:4]
        molCoords.append(coords)

    return molCoords

# get the Murckoscaffold of a molecule
def GetMurckoScaffold(mol):
    #mol: rdkit RWMol or Mol
    from rdkit.Chem.Scaffolds import MurckoScaffold

    scaffold = MurckoScaffold.MakeScaffoldGeneric(mol)

    #return scaffold rdkit.mol object
    return scaffold
