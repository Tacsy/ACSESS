#!/usr/bin/env python
#-*- coding: utf-8 -*-
import numpy as np

from rdkit import Chem
from rdkit.Chem import AllChem

import parallel
import mongoserver
import output


'''
this module include various previous modules that corresponds to chemical space
representation and selection, including:
    CellDiversity.py
    ChemGPS.py
    Distance.py
    Coords.py
    PCAmodule.py
    Kohonen.py
    BuildRandomKoho.py
    MakeSubsetDistance.py
'''

############################################################
#           Functions from Coords.py
############################################################

# Initialize coords
def CoordInit():
    global Coords

    # set coordination system
    if metric == 'mqn':
        from molproperty import CalcMQNs as Coords
    elif metric == 'autocorr':
        from molproperty import MoreauBrotoACVec as Coords
    elif metric == 'none':
        print 'No coordinate system set!'
    else:
        raise KeyError('Unknown metric specified in parameter file: ' + metric)

# Get coordinates
def SetCoords(mol):
    # set the coordinate of the molecule and return as a np.array object
    # if not, calculate first and then assign it as a property of the molecule
    if mol.HasProp('coords'):
        coord = mol.GetProp('coords')
    else:
        coord = Coords(mol)
        mol.SetProp('coords', coord)
    
    return np.array(coord)

# Drive MPI coordinate calculation
def ScatterCoords(mols):
    
    if not parallel.mpi:
        output.StartTimer('COORDS')
        for mol in mols:
            if not mol.HasProp('coords'):
                SetCoords(mol)
        output.EndTimer('COORDS')

        return None

    needCalc = [mol for mol in mols if not mol.HasProp('coords')]
    if len(needCalc) == 0:
        return None

    # Check if values are already computed
    if mprms.UseMongo:
        import mongoserver as mongo

        needSMI = [Chem.MolToSmiles(mol) for mol in needCalc]
        myMongo = mongo.LookupDB(metric, needSMI)
        for v, m in zip(myMongo, needCalc):
            if v is not None:
                m.SetProp('coords', v)

        oldn = len(needCalc)
        needCalc = [mol for mol in mols if not mol.HasProp('coords')]
        print 'Used %d memorized coordinate values' %(oldn - len(needCalc))

        if len(needCalc) == 0:
            return None

    output.StartTimer('COORDS')
    print 'Scattering chemical space coordinate calculation ...', len(needCalc)
    parallel.MyTask.SetFunction(MPICoordCalc)
    coords = parallel.MyTask.RunMPI(needCalc)

    for coord, mol in zip(coords, needCalc):
        mol.SetProp('coords', coord)

    # Update values in mongo database
    if mprms.UseMongo:
        needSMI = [Chem.MolToSmiles(mol) for mol in needCalc]
        mongo.UpdateDB(metric, {s: m.GetProp('coords') for s, m in zip(needSMI, needCalc)})
        print '%d new memorized coordinate values' %len(needCalc)

    output.EndTimer('COORDS')

# Function for scattering coordinate calculations
@parallel.MPIScatter
def MPICoordCalc(smistring):
    mol = Chem.MolFromSmiles(smistring)

    return SetCoords(mol)
