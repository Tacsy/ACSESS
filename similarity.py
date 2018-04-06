#!/usr/bin/env python
#-*- coding: utf-8 -*-
import numpy as np
import random

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs

fpcodes = {}

mcodes = {}

mfuncs = {}

#########################
# Module initialization #
#########################

fp = 'path'
measure = 'tanimoto'
mcode = None
fpcode = None
mfunc = None

def Init():

    global mcode, fpcode, mfunc
    
    #check requested fingerprint existance
    if not fpcode.has_key(fp):
        print 'Unrecognized fingerprint (mprms.fp): ' + fp
        print 'Implemented fingerprints: ', fpcodes.keys()
        raise KeyError('Unknown fingerprint')

    #check requested similarity measure existance
    if not mcodes.has_key(measure):
        print 'Unrecognized similarity measure (mprms.measure): ' + measure
        print 'Implemented similarity measures: ', mcodes.keys()
        raise KeyError('Unknown similarty measure')

    #set fingerprint type and similarity measure
    mcode = mcodes[measure]
    fpcode = fpcodes[fp]
    mfunc = mfuncs[measure]


def FPDB(mols,fpcode):
    '''
import random
Generate fingerprint database based on given fpcode
    '''
    return [fpcode(mol) for mol in mols]


def GenFPSimMatrix(mols,mfunc,fpcode):
    '''
    Compute similarity matrix using given fingerprint and similarity measure
    
    Note: lingo will produce fatal errors for molecules of less than 3 atoms.
    We either need special cases to handle this or just not use lingo.
    '''
    #get the fingerprint database
    fps = FPDB(mols,fpcode)

    n = len(mols)
    sim = np.zeros((n,n))
    #diagonal terms are set to be -1
    for i in range(n):
        sim[i,i] = -1
    #off-diagonal terms
    for i in range(1,n):
        for j in range(i):
            sim_ij = mfunc[fps[i],fps[j]]
            sim[i,j] = sim_ij
            sim[j,i] = sim[i,j]

    return sim


def FPMaximin(mols,nMol,sim=None,mfunc=None,fpcode=None):
    '''
    Maximin maximum diversity selection using fingerprint similarity as a
    measure
    It's currently coded to work with similarity rather than dissimilarity
    also only for things where we don't have descriptors, just pairwise distance
    
    If no similarity matrix is specified, it will calculate simialarity matrix
    first (i.e., the explicit similarity matrix is not stored)
    '''
    if sim is None:
        #no explicit similarity matrix, use "direct" version
        #different from previous version of Maximin in OPENEYE-ACSESS
        if mfunc is None or fpcode is None:
            print 'Need to assign fingerprint and similarity measure'
            raise KeyError('Unassigned parameters')
        else:
            sim = GenFPSimMatrix(mols,mfunc,fpcode)
    
    if nMol > len(mols):
        print 'ERROR: requested subset is larger than parent library.'
        raise Exception
    selected = [False] * len(mols)
    newmols = []
    
    for i in range(nMol):
        if i == 0:
            minsimmol = random.randrange(0,len(mols),1)
            minsim = 0.0
        else:
            minsim = 2.0
            for i in range(len(mols)):
                if selected[i]:
                    continue
                maxsim = -1.0
                for j in newmols:
                    if sim[i,j] > maxsim:
                        maxsim = sim[i,j]
                if maxsim < minsim:
                    minsim = maxsim
                    minsimmol = i
        newmols.append(minsimmol)
        selected[minsimmol] = True
    for i in range(len(newmols)):
        newmols[i] = mols[newmols[i]]

    return newmols

def NNSimilarity(mols,sim=None,mfunc=None,fpcode=None):
    '''
    Calculate average nearest neighbor similarity
    (a diverse library should MINIMIZE this function)
    '''
    div = 0.0
    if sim is None:
        #no explicit similarity matrix, use "direct" version
        #different from previous version of NNSimilarity in OPENEYE-ACSESS
        if mfunc is None or fpcode is None:
            print 'Need to assign fingerprint and similarity measure'
            raise KeyError('Unassigned parameters')
        else:
            sim = GenFPSimMatrix(mols,mfunc,fpcode)
        for imol in range(len(mols)):
            div += sim[imol].max()
    else:
        #with storaged similarity matrix version
        for imol in range(len(mols)):
            div += sim[imol].max()
    
    return div/len(mols)
