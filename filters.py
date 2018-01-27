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

############################################################
#       Functions from QiuFilter.py
############################################################

# ensure molecules do not have ring with ring size N
def CheckRingSize(molList, ringSize):
    # get mol list and return list of filterd molecules
    filtered = []
    
    for mol in molList:
        dumpMol = False
        # calculate smallest set of rings (SSR)
        ssr = Chem.GetSymmSSSR(mol)
        for i in range(len(ssr)):
            if len(list(ssr[i])) == ringSize:
                dumpMol = True

        if dumpMol:
            filtered.append(mol)
    
    return filtered

# ensure molecules do not have a bond with bondorder N
def CheckBondOrder(molList, bondOrder):
    # get mol list and return list of filtered molecules
    # bondOrder should be double, 1.0 for single, 1.5 for aromatic, 2.0 for
    # double, 3.0 for triple
    filtered = []

    for mol in molList:
        dumpMol = False
        # loop over bonds in molecule and check bond order
        for bond in mol.GetBonds():
            if bond.GetBondTypeAsDouble() == bondOrder:
                dumpMol = True

        if dumpMol:
            filtered.append(mol)

    return filtered
