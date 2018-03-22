#!/usr/bin/env python
#-*- coding: utf-8 -*-
from rdkit import Chem
from rdkit.Chem import AllChem

from exceptions import *


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

Note from Jos:
    I don't think it is possible to include and combine every filter into this module
    since this module would become unmanageable. We could make a filters subfolder to 
    order them. 
'''
'''
Some clearification:
    all filter-related function are acting on single molecule and return a
    bool value to indicate whether a molecule should be filtered out (True)
    or not (False)
'''

ActiveFilters = dict()
FilterFlavor = 'GDB'
MAXTRY=10

############################################################
#       Functions from Filter.py
############################################################

def FilterInit():
    global ActiveFilters, FilterFlavor

    if FilterFlavor=='GDB':
        from filters import GDBFilters as FilterFlavor
    else:
        raise TypeError('FilterFlavor not recognized')
    # so every filterflavor-module has a dictionary which is formatted as:
    # filterfunction might be a filterclass with a __call__ attribute though.
    # {'filtername1':filterfunction1, 'filtername2':filterfunction2, ... }
    ActiveFilters.update(FilterFlavor.AllFilters)

    return

def FixAndFilter(mol):

    # HERE AN EXTRY MAXTRY LAYER SHOULD BE ADDED AS IN THE OLD VERSION.
    # THIS IS JUST TO TEST
    changed=False
    failure=False
    for _ in xrange(MAXTRY):
        # in old version: here prepare for 2D filters

        # run through all filters:
        for filtername in ActiveFilters:
            failure = ActiveFilters[filtername](mol)
            if failure:
                try:
                    success= ActiveFilters[filtername].Fix(mol)
                except MutateFail:
                    success=False
                    changed=True
            if not success: return changed, failure
            else:
                changed=True
                #cn.Finalize(mol)
                break
        if not failure: break #i.e. all filters passed without problems
    return changed, failure




















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
