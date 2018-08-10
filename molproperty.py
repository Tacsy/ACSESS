#!/usr/bin/env python
#-*- coding: utf-8 -*-
import numpy as np

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.six.moves import cPickle
from rdkit.six import iteritems
import mprms

import math
from collections import defaultdict

import os

'''
this module include various previous modules that corresponding to molecular
property calculation, including:
    Properties.py
    SAScore.py
    SMCM.py
    MQN.py
'''
maxBonds = 7

def AutoCorr2D(mol):
    from rdkit.Chem import rdMolDescriptors
    vector = rdMolDescriptors.CalcAUTOCORR2D(mol)
    return vector

def AutoCorrMordred(mol):
    from mordred import Calculator, Autocorrelation
    calc = Calculator()
    # ATS, ATSC, AATS, AATSC ?
    descriptor = Autocorrelation.ATS
    calc.register(descriptor)
    res = calc(mol)
    res = res.fill_missing()
    # Z: atomic num, pe=pauling electronegativity, p=polarizability, x=unweighted(identity), v=vdw-volume
    # dv= nValence
    props= ['Z', 'pe', 'p', 'dv' ]
    keys = [ 'ATS{d}{p}'.format(d=d, p=p) for d in range(maxBonds+1) for p in props ]
    #print "keys:", keys
    res = { k:v for k, v in res.asdict().iteritems() if k in keys }
    for key in keys:
        if not key in res:
            print key
    #print "res:", res
    vector = [ value for (key, value) in sorted(res.items())]
    #print "len(vector):", len(vector)
    return vector

def MoreauBrotoPyBioMed(mol):
    #from PyBioMed import Pymolecule
    #_mol = Pymolecule.PyMolecule()
    #smi = Chem.MolToSmiles(mol)
    #_mol.ReadMolFromSmile('C(c1ccccc1)C(C)O')
    #vector = _mol.GetMoreauBroto().values()
    from PyBioMed.PyMolecule.moreaubroto import GetMoreauBrotoAuto
    vector = GetMoreauBrotoAuto(mol)
    # vector is a dict so return it as a list sorted based on the dict keys
    return [ value for (key, value) in sorted(vector.items())]


############################################################
# Assign polarizabilities to every atom.
# Units are Angstrom^3/2. We assume all valences are filled.
#
# This will work without explicit hydrogens, but it will
# obviously not assign properties to implicit hydrogens.
#
# Values are from Table I of
# Miller and Savchik, JACS 101(24) 7206-7213, 1979.
# dx.doi.org/10.1021/ja00518a014
############################################################
def AssignAtomicPolarizability(mol):

    for atom in mol.GetAtoms():
        #get degree for specific molecule
        nBonds = atom.GetTotalDegree()

        #Hydrogen (H)
        if atom.GetAtomicNum() == 1:
            atom.SetDoubleProp('polarizability', 0.314)

        #Carbon (C)
        elif atom.GetAtomicNum() == 6:
            if nBonds == 4:
                atom.SetDoubleProp('polarizability', 1.294)
            elif nBonds == 2:
                atom.SetDoubleProp('polarizability', 1.393)
            elif nBonds == 3:
                if atom.GetNumExplicitHs() + atom.GetNumImplicitHs() > 0:
                    atom.SetDoubleProp('polarizability', 1.428)
                else:
                    '''
                    in this part I employ a different logic than that in
                    previous version of ACSESS and may be consult with Aaron for
                    further details
                    '''
                    cross = True
                    for nbor in atom.GetNeighbors():
                        if not nbor.GetIsAromatic():
                            cross = False
                            break
                    if cross:
                        atom.SetDoubleProp('polarizability', 1.800)
                    else:
                        atom.SetDoubleProp('polarizability', 1.428)

        #Nitrogen (N)
        elif atom.GetAtomicNum() == 7:
            if atom.GetIsAromatic():
                if nBonds == 2:
                    atom.SetDoubleProp('polarizability', 1.262)
                else:
                    atom.SetDoubleProp('polarizability', 1.220)
            else:
                if nBonds == 1:
                    atom.SetDoubleProp('polarizability', 1.304)
                else:
                    atom.SetDoubleProp('polarizability', 1.435)

        #Oxygen (O)
        elif atom.GetAtomicNum() == 8:
            if atom.GetIsAromatic():
                atom.SetDoubleProp('polarizability', 1.099)
            else:
                if nBonds == 1:
                    atom.SetDoubleProp('polarizability', 1.216)
                else:
                    atom.SetDoubleProp('polarizability', 1.290)

        #Sulfur (S)
        elif atom.GetAtomicNum() == 16:
            if atom.IsAromatic():
                atom.SetDoubleProp('polarizability', 2.982)
            elif nBonds == 2:
                atom.SetDoubleProp('polarizability', 3.496)
            else:
                atom.SetDoubleProp('polarizability', 3.967)

        #Halogens
        elif atom.GetAtomicNum() == 9:
            atom.SetDoubleProp('polarizability', 1.046)
        elif atom.GetAtomicNum() == 15:
            atom.SetDoubleProp('polarizability', 3.000)
        elif atom.GetAtomicNum() == 17:
            atom.SetDoubleProp('polarizability', 3.130)
        elif atom.GetAtomicNum() == 35:
            atom.SetDoubleProp('polarizability', 5.577)
        elif atom.GetAtomicNum() == 53:
            atom.SetDoubleProp('polarizability', 8.820)

#Iridium (I)
#This param value was obtained by fitting the above known tau values
#In general polarizability increases with atomic number so we used
#linear fit to get the value
#This is a crudest approx so could be wrong!
        elif atom.GetAtomicNum() == 77:
            atom.SetData('polarizability', 12.77)

        else:
            raise KeyError('No polarizabilities for atomic number' +
                           str(atom.GetAtomicNum()))


############################################################
# Topological steric effect index (TSEI)
# Cao and Liu, J Chem Inf Comput Sci 44(2), 678-687, 2004.
# dx.doi.org/10.1021/ci034266b


# Covalent radii (in angstroms)
# These are the 2008 values from CSD, except for carbon which is from Cao
#
# This returns numbers a little different than those returned by Marvin.
# If this is to be used for serious business, we should figure out why.
# First check would probably be covalent radius values
# Second would be path termination length
############################################################
def AssignTSEI(mol):
    # all mol parsed into this function should be hydrogen added

    radius = {
        1: 0.315,
        5: 0.84,
        6: 0.772,
        7: 0.711,
        8: 0.662,
        9: 0.57,
        14: 1.11,
        15: 1.06,
        16: 1.05,
        17: 1.024,
        35: 1.14,
        53: 1.33,
        77: 1.42
    }
    nAtoms = mol.GetNumAtoms()
    '''
    this is a revision part of code that deviate from previous code in ACSESS,
    we should consult with Aaron for more detail
    '''

    #Need to get shortest path between each pair of atoms
    #this could be faster by following paths outward from each atom
    #if atom is hydrogen, the TSEI is set as 0.0 directly
    for i in xrange(nAtoms):
        iAtom = mol.GetAtomWithIdx(i)
        if iAtom.GetAtomicNum() == 1:
            iAtom.SetDoubleProp('TSEI', 0.0)
            continue
        for j in xrange(i + 1, nAtoms):
            jAtom = mol.GetAtomWithIdx(j)
            if jAtom.GetAtomicNum() == 1:
                continue
            else:
                #calculate shoretest path
                path = Chem.GetShortestPath(mol, i, j)
                length = 0.0
                for k in path:
                    kAtom = mol.GetAtomWithIdx(k)
                    if k == 0 or k == len(path) - 1:
                        length += radius[kAtom.GetAtomicNum()]
                    else:
                        length += 2.0 * radius[kAtom.GetAtomicNum()]

                #factor of 2 cancels length vs. radius
                iRad = radius[iAtom.GetAtomicNum()]
                jRad = radius[jAtom.GetAtomicNum()]
                if iAtom.HasProp('TSEI'):
                    iAtom.SetDoubleProp('TSEI',
                                        iAtom.GetDoubleProp('TSEI') +
                                        (2.0 * jRad / length)**3)
                else:
                    iAtom.SetDoubleProp('TSEI', (2.0 * jRad / length)**3)
                if jAtom.HasProp('TSEI'):
                    jAtom.SetDoubleProp('TSEI',
                                        jAtom.GetDoubleProp('TSEI') +
                                        (2.0 * iRad / length)**3)
                else:
                    jAtom.SetDoubleProp('TSEI', (2.0 * iRad / length)**3)


############################################################
# Calculate autocorrelation vector for atomic properties
# Method first described in Moreau and Broto. Nouv. J. Chim.1980, 4, 757-764.
# A more recent reference is dx.doi.org/10.1002/poc.610061008
############################################################
def AutoCorrelation(mol, props):

    atoms = mol.GetAtoms()

    #Initializes the autocorrelation vector for each property
    #To access the component for property x at bond separation y, use
    #ACVector[x][y]
    ACVector = [([0.0] * (maxBonds + 1)) for prop in props]

    #Retrive all properties for each atom
    #atomProps[i][j] gives property j for itom i
    atomProps = np.array(
        [[atom.GetDoubleProp(prop) for prop in props] for atom in atoms])

    #construct autocorrelation vector
    nAtoms = mol.GetNumAtoms()
    for i in xrange(nAtoms):
        iAtom = mol.GetAtomWithIdx(i)
        for j in xrange(i, nAtoms):
            jAtom = mol.GetAtomWithIdx(j)

            #get path length between this pair of atoms
            if i != j:
                pathlen = len(Chem.GetShortestPath(mol, i, j)) - 1
                '''
                this part is little bit confusing, we should consult with Aaron
                for his advice on this
                '''
                if pathlen > maxBonds:
                    pathlen = 0
            else:
                pathlen = 0

            #add component to autocorrelation vectors
            for iprop, prop in enumerate(props):
                ACVector[iprop][pathlen] += (
                    atomProps[i][iprop] * atomProps[j][iprop])

    return ACVector

def MoreauBrotoACVector(mol0):

    mol = Chem.AddHs(mol0)
    AssignTSEI(mol)
    AssignAtomicPolarizability(mol)
    mol.ComputeGasteigerCharges()
    for atom in mol.GetAtoms():
        atom.SetDoubleProp('anum', atom.GetAtomicNum())
        atom.SetDoubleProp('one',1)

    props = ['TSEI', 'anum', 'one', '_GasteigerCharge', 'polarizability']
    result = np.array( AutoCorrelation(mol, props) )

    #Normalize by typical carbon values
    result[0,:]=result[0,:]/(2.3*2.3)
    result[1,:]=result[1,:]/36.0
    result[3,:]=result[3,:]/.005
    result[4,:]=result[4,:]/(1.2*1.2)

    result.shape=( np.product(result.shape),)
    return result

############################################################
# Calculation of Molecular Quantum Numbers as given by Nguyen
# in ChemMedChem 4: 1803-5(2009)
# it is also a built-in function RDkit
############################################################


def CalcMQNs(mol):
    molAddH = Chem.AddHs(mol)
    MQNs = AllChem.MQNs_(molAddH)

    return MQNs


############################################################
# Calculation of synthetic accessibility score as described in :
# Journal of Cheminformatic 1:8 (2009) with Several small modifications
# Peter ertl & Greg landrum, Sept 2013, also in RDKit/Contrib/SA_Score
# need to retain original copyright notice
############################################################

_fscores = None


def ReadFragScores(name='fpscores'):
    import gzip
    global _fscores
    #generate the full path filename
    if name == "fpscores":
        name = os.path.join(os.path.dirname(__file__), name)

    _fscores = cPickle.load(gzip.open('%s.pkl.gz' % name))
    outDict = {}
    for i in _fscores:
        for j in range(1, len(i)):
            outDict[i[j]] = float(i[0])
    _fscores = outDict


def NumBridgeheadsAndSpiro(mol, ri=None):
    nSpiro = AllChem.CalcNumSpiroAtoms(mol)
    nBridgehead = AllChem.CalcNumBridgeheadAtoms(mol)

    return nBridgehead, nSpiro


def CalcSAScore(mol):
    if _fscores is None:
        ReadFragScores()

    #fragment score
    fp = AllChem.GetMorganFingerprint(
        mol, 2)  #<- 2 is the *radius* of the circular fingerprint
    fps = fp.GetNonzeroElements()
    score1 = 0.0
    nf = 0
    for bitId, v in iteritems(fps):
        nf += v
        sfp = bitId
        score1 += _fscores.get(sfp, -4) * v
    score1 /= nf

    #features score
    nAtoms = mol.GetNumAtoms()
    nChiralCenters = len(
        Chem.FindMolChiralCenters(mol, includeUnassigned=True))
    ri = mol.GetRingInfo()
    nBridgehead, nSpiro = NumBridgeheadsAndSpiro(mol, ri)
    nMacrocycles = 0
    for x in ri.AtomRings():
        if len(x) > 8:
            nMacrocycles += 1

    sizePenalty = nAtoms**1.005 - nAtoms
    stereoPenalty = math.log10(nChiralCenters + 1)
    spiroPenalty = math.log10(nSpiro + 1)
    bridgePenalty = math.log10(nBridgehead + 1)
    macrocyclePenalty = 0.0
    # -----------------------------
    # This differs from the paper, which defines:
    #   macrocyclePenalty = math.log10(nMacrocycles+1)
    # This form generates better results when 2 or more macrocycles are present
    if nMacrocycles > 0:
        macrocyclePenalty = math.log10(2)

    score2 = 0.0 - sizePenalty - stereoPenalty - spiroPenalty - bridgePenalty - macrocyclePenalty

    # correction for the fingerprint density
    # not in the original publication
    # to make highly symmetrical molecules easier to synthesize
    score3 = 0.0
    if nAtoms > len(fps):
        score3 = math.log(float(nAtoms) / len(fps)) * 0.5

    sascore = score1 + score2 + score3

    # need to transform "raw" value into scale between 1 and 10
    minv = -4.0
    maxv = 2.5
    sascore = 11.0 - (sascore - minv + 1) / (maxv - minv) * 9.0
    # smooth the 10-end
    if sascore > 8.0:
        sascore = 8.0 + math.log(sascore - 8.0)
    if sascore > 10.0:
        sascore = 10.0
    elif sascore < 1.0:
        sascore = 1.0

    return sascore


############################################################
# Calculation of molar LogP and molar refractivity as given by
# Crippen in J Chem Inf Comp Sci 1999 39(5) 868-873
# Also provided as built-in function in RDKit
############################################################


def CalcCrippenLogPAndMR(mol):
    LogP, MR = AllChem.CalcCrippenDescriptors(mol)

    return LogP, MR


############################################################
# Determine a molecule type according to ring analysis,
# belongs to the following categories: aromatic, heteroaromatic,
# fused-heterocycle, fused-cycle, polyheterocyclic, monoheterocyclic,
# polycyclic, monocyclic, heteroacyclic, acyclic
############################################################
'''
currently it is shown that this function is not called in any other modules, so
we will keep it on hold for further revision since it could be useful in some
cases.
'''


def ClassifyMolType(mol):
    return 0
