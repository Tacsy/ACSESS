#!/usr/bin/env python
#-*- coding: utf-8 -*-
from filters import NewFilter, NewPatternFilter
from GDBFilters import FixByRemovingHeteroatoms
from rdkit import Chem
from rdkithelpers import *
from copy import deepcopy
from itertools import combinations
ExtraFilters = dict()
extraSmarts  = []
extraSmartsAromatic = False



# -----------------
ExtraFilters['ExtraSmarts'] = NewFilter('ExtraSmarts')
def HasSmarts(mol):
    if not extraSmartsAromatic:
        try:
            Chem.Kekulize(mol, True)
        except ValueError as e:
            return 'unkekulizable'
    smarts = [ Chem.MolFromSmarts(p) for p in extraSmarts ]
    answer = any( mol.HasSubstructMatch(smart) for smart in smarts )
    Chem.SetAromaticity(mol)
    if answer: return False
    else: return 'no smarts match'
ExtraFilters['ExtraSmarts'].SetFilterRoutine(HasSmarts)

# -------------------
ExtraFilters['ringshare3'] = NewFilter('ringshare3')
def ringshare3(mol):
    rings = mol.GetRingInfo().AtomRings()
    if not rings: return False
    combis= combinations(rings, 2)
    lintersect = lambda combi: len(set(combi[0]) & set(combi[1]))
    intersects = map(lintersect, combis)
    if not intersects: return False
    if max(intersects)>2:
        return 'two rings share 3 atoms'
    else:
        return False
ExtraFilters['ringshare3'].SetFilterRoutine(ringshare3)

# -----------------
# EXTRA AROMATICITY FILTER
ExtraFilters['aromatic'] = NewFilter('aromatic')

def HasAromaticity(oldmol):
    mol = deepcopy(oldmol)
    #a1 = Chem.MolFromSmarts("C1=C-C=C-C=C1")
    a2 = Chem.MolFromSmarts("C1~*-cc-*~*1")
    a3 = Chem.MolFromSmarts("O=C1-*(@[R]):,=*(@[R])-C(-*:,=-*1)=O")
    a4 = Chem.MolFromSmarts("[R]=*1-*=*-*(=*)-*=,:*1")
    a5 = Chem.MolFromSmarts("*=,:*1:*:*:*(=,:*)*:*1")
    A = [a2, a3, a4, a5]
    if any(mol.HasSubstructMatch(a) for a in A):
        return False
    #if len(mol.GetAromaticAtoms()) or mol.HasSubstructMatch(a):
    #    return False
    return 'not aromatic'

ExtraFilters['aromatic'].SetFilterRoutine(HasAromaticity)

# -----------------
# Sphericity filters. Other Descriptor 3D filters could be added
# preferably to this particular function to avoid repeatedly
# calculating 3D coordinates
maxSphericity = 0.3
ExtraFilters['sphericity'] = NewFilter('sphericity')
def Sphericity(mol):
    from rdkit.Chem import AllChem
    from rdkit.Chem import Descriptors3D
    m = Chem.AddHs(mol)
    try: # this function can raise a ValueError for unkekulizable molecules
        AllChem.EmbedMolecule(m)
    except ValueError:
        return False
    try:
        AllChem.UFFOptimizeMolecule(m)
    except ValueError:
        print "ValueError with UFFOptimize, mol:", Chem.MolToSmiles(m), Chem.MolToSmiles(mol)
        return False
    sphericity = Descriptors3D.SpherocityIndex(m)
    if sphericity > maxSphericity:
        return 'maxSphericity {}'.format(sphericity)
    return False
ExtraFilters['sphericity'].SetFilterRoutine(Sphericity)

###########################################################

newfilt = NewPatternFilter('Polycyclic3')
newfilt.SetFilterPattern(Chem.MolFromSmarts('*@N(@*)@C@=*'))
ExtraFilters['Polycyclic3'] = newfilt

newfilt = NewPatternFilter('nnn')
newfilt.SetFilterPattern(Chem.MolFromSmarts("n~n~[#7]"))
#newfilt.SetFixRoutine(FixByRemovingHeteroatoms)
ExtraFilters['nnn'] = newfilt

############################################################
#       Functions from QiuFilter.py
############################################################

ringSizes = [3, 4, 8]
bondOrders = [3.0]
ExtraFilters['qiu1'] = NewFilter('checkringsize')
ExtraFilters['qiu2'] = NewFilter('checkbondorder')


# ensure molecule do not have ring with ring size N
def CheckRingSize(mol):
    global ringSizes
    # calculate smallest set of rings (SSR)
    ssr = Chem.GetSymmSSSR(mol)
    for i in range(len(ssr)):
        if len(list(ssr[i])) in ringSizes:
            return True
    return False


# ensure molecule do not have a bond with bondorder N
def CheckBondOrder(mol):
    # get mol and return bool for filter or not
    # bondOrder should be double, 1.0 for single, 1.5 for aromatic, 2.0 for
    # double, 3.0 for triple
    # loop over bonds in molecule and check bond order
    for bond in mol.GetBonds():
        if bond.GetBondTypeAsDouble() in bondOrders:
            return True
    return False


ExtraFilters['qiu1'].SetFilterRoutine(CheckRingSize)
ExtraFilters['qiu2'].SetFilterRoutine(CheckBondOrder)

