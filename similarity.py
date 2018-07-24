#!/usr/bin/env python
#-*- coding: utf-8 -*-
import numpy as np
import random

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols

from rdkit.Chem import MACCSkeys

# Global variables
fp = 'rdkit'
measure = 'tanimoto'
morganradius = 2
dmetric = None
fpcode= None

# old OE codes:
#fpcodes={'lingo':OEFPType_Lingo,
#         'maccs':OEFPType_MACCS166,
#         'path':OEFPType_Path}

fpcodes= {'rdkit': FingerprintMols.GetRDKFingerprint,
          'maccs': MACCSkeys.GenMACCSKeys,
          'morgan': lambda mol:AllChem.GetMorganFingerprintAsBitVect(mol, morganradius),
          'atompairs': lambda mol:Chem.GetAtomPairFingerPrint(mol, 2)
         }

# tversky(F1, F2)= F1@F2 / ( a*sum(F1) + b*sum(F2) - (1-a-b)* F1@F2 )
dmetrics = {'tanimoto': DataStructs.TanimotoSimilarity,  # a, b = 1, 1
           'dice': DataStructs.DiceSimilarity, # a, b = 0.5, 0.5
           'cosine': DataStructs.CosineSimilarity,
           'tversky': lambda m1, m2: DataStructs.TverskySimilarity(m1, m2, 0.5, 0.5),
           'sokal': DataStructs.SokalSimilarity }

#########################
# Module initialization #
#########################


def Init():

    global dmetric, fpcode

    #check requested fingerprint existance
    if not fpcodes.has_key(fp):
        print 'Unrecognized fingerprint (mprms.fp): ' + fp
        print 'Implemented fingerprints: ', fpcodes.keys()
        raise KeyError('Unknown fingerprint')

    #check requested similarity measure existance
    if not dmetrics.has_key(measure):
        print 'Unrecognized similarity measure (mprms.measure): ' + measure
        print 'Implemented similarity measures: ', dmetrics.keys()
        raise KeyError('Unknown similarty measure')

    #set fingerprint type and similarity measure
    dmetric = dmetrics[measure]
    fpcode = fpcodes[fp]

    print "## Similarity scheme will be used:"
    print "## - Fingerprint: {}".format(fpcode)
    print "## - SimilarityM: {}".format(measure)

    return 


def FPDB(mols, fpcode):
    '''
import random
Generate fingerprint database based on given fpcode
    '''
    return [fpcode(mol) for mol in mols]


def GenFPSimMatrix(mols):
    '''
    Compute similarity matrix using given fingerprint and similarity measure
    
    Note: lingo will produce fatal errors for molecules of less than 3 atoms.
    We either need special cases to handle this or just not use lingo.
    '''
    #get the fingerprint database
    fps = FPDB(mols, fpcode)

    n = len(mols)
    sim = np.zeros((n, n))
    #diagonal terms are set to be -1
    for i in range(n):
        sim[i, i] = -1
    #off-diagonal terms
    for i in range(1, n):
        for j in range(i):
            sim_ij = dmetric(fps[i], fps[j])
            sim[i, j] = sim_ij
            sim[j, i] = sim[i, j]
    return sim


def FPMaximin(mols, nMol, sim=None):
    '''
    Maximin maximum diversity selection using fingerprint similarity as a
    measure
    It's currently coded to work with similarity rather than dissimilarity
    also only for things where we don't have descriptors, just pairwise distance
    
    If no similarity matrix is specified, it will calculate simialarity matrix
    first (i.e., the explicit similarity matrix is not stored)
    '''
    if sim is None:
        #different from previous version of Maximin in OPENEYE-ACSESS
        if dmetric is None or fpcode is None:
            print 'Need to assign fingerprint and similarity measure'
            raise KeyError('Unassigned parameters')
        else:
            sim = GenFPSimMatrix(mols)

    if nMol > len(mols):
        print 'ERROR: requested subset is larger than parent library.'
        raise Exception

    selected = [False] * len(mols)
    newmols = []
    for i in range(nMol):
        if i == 0:
            # initially pick a random molecule
            minsimmol = random.randrange(0, len(mols), 1)
            minsim = 0.0
        else:
            minsim = 2.0
            for i in range(len(mols)):
                if selected[i]:
                    continue
                maxsim = -1.0
                for j in newmols:
                    if sim[i, j] > maxsim:
                        maxsim = sim[i, j]
                if maxsim < minsim:
                    minsim = maxsim
                    minsimmol = i
        newmols.append(minsimmol)
        selected[minsimmol] = True
    for i in range(len(newmols)):
        newmols[i] = mols[newmols[i]]

    return newmols


def NNSimilarity(mols, sim=None, average=False):
    '''
    New version:
    NNSimilarity = sum of maximum similarities (excluding identity)

    Older version:
    Calculate average nearest neighbor similarity
    (a diverse library should MINIMIZE this function)
    '''
    if sim is None:
        if dmetric is None or fpcode is None:
            print 'Need to assign fingerprint and similarity measure'
            raise KeyError('Unassigned parameters')
        else:
            sim = GenFPSimMatrix(mols)

    if average:
        avdiv = np.average(sim[np.tril_indices(len(sim), -1)])
        return avdiv
    else:
        div = 0.0
        for imol in range(len(mols)):
            div += sim[imol].max()
        return div / len(mols)
