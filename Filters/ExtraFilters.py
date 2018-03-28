#!/usr/bin/env python
#-*- coding: utf-8 -*-
from filters import NewFilter, NewPatternFilter
from rdkit import Chem

ExtraFilters=dict()

# EXTRA AROMATICITY FILTER
ExtraFilters['aromatic']=NewFilter('aromatic')
def HasAromaticity(mol):
    a = Chem.MolFromSmarts("*1=**=*=**1")
    if len(mol.GetAromaticAtoms()) or mol.HasSubstructMatch(a):
        return False
    return True
ExtraFilters['aromatic'].SetFilterRoutine(HasAromaticity)
#########################


############################################################
#       Functions from QiuFilter.py
############################################################

ringSizes=[3,4,8]
bondOrders=[3.0]
ExtraFilters['qiu1']=NewFilter('checkringsize')
ExtraFilters['qiu2']=NewFilter('checkbondorder')
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
