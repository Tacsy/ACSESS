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


