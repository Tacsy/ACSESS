#!/usr/bin/env python
''' This module acts as an interface between ACSESS and CINDES
    The expected input:
    - list of Chem.Mol objects
    - CINDES.INDES.procedures.run object with calc. specifications

    The result is that all Chem.Mol objects have an objective attribute

    by J.L. Teunissen
'''
import os,sys
import numpy as np

from rdkit import Chem

from helpers import Compute3DCoords, xyzfromrdmol
#s3d.UseBalloon=True
#s3d.S3DInit()

from CINDES.utils.molecule import SmiMolecule
table=None

def Init():
    from CINDES.INDES.inputreader import readfile
    from CINDES.INDES.procedures import BaseRun
    # here maybe a Run object can be created?
    param = readfile(open('INPUT','r'))
    run=BaseRun(**param)
    print run
    global table
    from CINDES.INDES.procedures import set_table
    table=set_table(run)
    return run

def calculate(rdmols, run, QH2=False, gen=0):
    # -2 prepare optionally for logging to file
    global table
    print "len(table)".format(len(table)),
    if len(table)<10: print table

    # -1 imports
    from CINDES.INDES.calculator import do_calcs, set_target_properties
    from CINDES.INDES.predictions import predictor
    from CINDES.INDES.construction import check_in_table

    # 1. RDKit.Chem.Mol(rdmol) objects have to have a molecular specification
    # CINDES.utils.molecule.SmiMolecule objects
    mols = [ SmiMolecule(Chem.MolToSmiles(rdmol, True)) for rdmol in rdmols ]
    for mol, rdmol in zip(mols, rdmols):
        #from Canonical import Finalize
        #Finalize(rdmol)
        mol.xyz= getXYZ(rdmol)
        mol.rdmol=rdmol
        #print "xyz:", mol.xyz
        mol.xyz2 = mol.xyz

    # check if already in database
    mols_todo, mols_nodo = check_in_table(mols, table, run.props, check_ignored=True)

    # Optional: perform prescreaning in a predictions.
    if run.predictions:
        from CINDES.INDES.procedures import get_property_table
        property_table = get_property_table(table, run)
        mols_nocal, mols_tocal, made_pred = predictor(
            run,
            property_table,
            mols_todo,
            mols_nodo,
            count=gen,
            )
    else:mols_nocal, mols_tocal = (mols_nodo, mols_todo)

    #2. do the calculations
    do_calcs(mols_tocal, run)
    mols = mols_tocal+mols_nocal

    #3. calculate final objectives
    set_target_properties( mols, run)

    # Optional. log results to screen
    if True:
        from CINDES.INDES.loggings import log_screen
        log_screen(mols)

    #4. set the results as attributes from the rdmols
    for mol in mols:
        mol.rdmol.SetDoubleProp('Objective', float(mol.Pvalue))

    from CINDES.INDES.loggings import log_table
    print "logging table..."
    log_table(mols, table)

    return

def getXYZ(rdmol):
    xyz = xyzfromrdmol(rdmol)
    return xyz


def xyzfromstring(string):
    # split to get the first line which has the number of atoms
    splitted=string.split('\n',2)
    nxyz=int(splitted[0])
    # take only as many atoms as are in one configuration:
    xyz = '\n'.join(splitted[2].split('\n',nxyz)[:-1])
    # add an empty line
    xyz += '\n\n'
    return xyz

if __name__=="__main__":
    import sys
    class Unbuffered(object):
        def __init__(self,stream):
            self.stream = stream
        def write(self,data):
            self.stream.write(data)
            self.stream.flush()
        def __getattr__(self,attr):
            return getattr(self.stream, attr)
    sys.stdout = Unbuffered(sys.stdout)

    import Filters as fl
    fl.FLInit()

    if sys.argv[1][-4:]=='.smi':
        print "testset SMILES file detected\n\tperforming test calculation"
        run=getrun()
        with open(sys.argv[1]) as f:
            smiles=[ line.strip('\n') for line in f.readlines() ]
            print smiles
        rdmols=[]
        for smi in smiles:
            mol=Chem.MolFromSmiles(smi, True)
            rdmols.append(mol)
        calculate(rdmols, run, log=True)
    print "done"

