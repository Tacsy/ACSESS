#!/usr/bin/env python
#-*- coding: utf-8 -*-
from rdkit import Chem
from rdkit.Chem import AllChem


'''
This filter.py contains all the possible filters in the previous version of
ACSESS, in order to make the code in the cleaner way.

previous filters included:
    Filters.py
    GDBFilters.py
    DruglikeFilters.py
    OEFilters.py
    SMUFilters.py
    Filters_BC.py
    QiuFilter.py
'''
'''
Some clearification:
    all filter-related function are acting on single molecule and return a
    bool value to indicate whether a molecule should be filtered out (True)
    or not (False)
'''
############################################################
#       Functions from Filter.py
############################################################

# ensure molecule has specific pattern
def CheckSubstructure(mol, patterns):
    # get mol and return bool for filter or not
    # patterns: list of pattern for substructure search 
    dumpMol = False
    # for substructure search
    for pattern in xrange(len(patterns)):
        if mol.HasSubstructMatch(pattern):
            dumpMol = True
            break

    return dumpMol

############################################################
#       Functions from QiuFilter.py
############################################################

# ensure molecule do not have ring with ring size N
def CheckRingSize(mol, ringSize):
    # get mol and return bool for filter or not
        
    dumpMol = False
    # calculate smallest set of rings (SSR)
    ssr = Chem.GetSymmSSSR(mol)
    for i in range(len(ssr)):
        if len(list(ssr[i])) == ringSize:
            dumpMol = True

    return dumpMol

# ensure molecule do not have a bond with bondorder N
def CheckBondOrder(mol, bondOrder):
    # get mol and return bool for filter or not  
    # bondOrder should be double, 1.0 for single, 1.5 for aromatic, 2.0 for
    # double, 3.0 for triple

    dumpMol = False
    # loop over bonds in molecule and check bond order
    for bond in mol.GetBonds():
        if bond.GetBondTypeAsDouble() == bondOrder:
            dumpMol = True

    return dumpMol
