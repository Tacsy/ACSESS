#!/usr/bin/env python
#-*- coding: utf-8 -*-
from filters import NewFilter, NewPatternFilter
from rdkit import Chem
from rdkithelpers import *

ExtraFilters=dict()

# EXTRA AROMATICITY FILTER
ExtraFilters['aromatic']=NewFilter('aromatic')
def HasAromaticity(mol):
    a = Chem.MolFromSmarts("C1=C-C=C-C=C1")
    b = Chem.MolFromSmarts("c1ccccc1")
    #if mol.HasSubstructMatch(a) or mol.HasSubstructMatch(b):
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

ExtraFilters['quinoid']=NewFilter('notquinoid')
def findQuinoid(mol, strict=False):
    Chem.Kekulize(mol)

    # quionoid matches
    # the 12 match urges for at least one extra conjugated bond that is not terminal
    ss12=Chem.MolFromSmarts('[O,o]:,=[#6]~[#6](:,=[O,o])~[#6]=,:*~*')
    ss14=Chem.MolFromSmarts('[O,o]:,=[#6]-,:[#6,#7]=,:[#6,#7]-,:[#6]:,=[O,o]')
    ss16=Chem.MolFromSmarts('[O,o]:,=[#6]-,:[#6]=,:[#6,#7]-,:[#6,#7]=,:[#6,#7]-,:[#6]:,=[O,o]')
    ss18=Chem.MolFromSmarts('[O,o]:,=[#6]-,:[#6]=,:[#6,#7]-,:[#6,#7]=,:[#6,#7]-,:[#6,#7]=,:[#6]-,:[#6]:,=[O,o]')
    sss=[    ((0, 3), ss12),
             ((0, 5), ss14),
             ((0, 7), ss16),
             ((0, 9), ss18),]

    qmatch=0
    ssss=[ss12, ss14, ss16, ss18]
    if True:
        matches=set()
        for qinds_p, ss in sss:
            for match in mol.GetSubstructMatches(ss):
                for idx, ma in zip(match, GetIAtoms(match, mol)):
                    if ma.GetAtomicNum()==8:
                        matches.add(idx)
        if len(matches)<2:
            print "no match:", Chem.MolToSmiles(mol)
            return True
        elif len(matches)>2:
            print "multiple matches:", Chem.MolToSmiles(mol)
            return True
        mol.SetProp('quinoid_indices', " ".join(map(str, list(matches))))
    for ss in ssss:
        if mol.HasSubstructMatch(ss): qmatch+=1
    if True: # molecule needs to have at least conjugated path length of 5 atoms or a conjugated 5/6 ring
        rmatch=0
        # some extra requirements:
        ring5a = Chem.MolFromSmarts('*1:,=*:,-*:,=*~*1') #i.e. 5 ring with two internal double bonds
        ring5b = Chem.MolFromSmarts('*1~*(=*)~*=*~*1') #i.e. 5 ring with double bonds. one internal one external
        ring6a = Chem.MolFromSmarts('*1~*=,:*-,:*=,:*~*1') #i.e. 6 ring with two internal double bonds in a conjugated pattern
        ring6b = Chem.MolFromSmarts('*1~*=,:*-,:*(:,=*)~*~*1')
        conjp = Chem.MolFromSmarts('*=*-*=*:,-*:,=*')
        extrar=[ring5a, ring5b, ring6a, ring6b, conjp]
        for r in extrar:
            if mol.HasSubstructMatch(r): rmatch+=1
        requirement = qmatch>0 and rmatch>0
    else:requirement = qmatch>0
    if True: #set maximum number of carbonyl groups
        nmax_co=4
        n_co=0
        carbonyl = Chem.MolFromSmarts('[#6]=,:[#8;X1]')
        for match in mol.GetSubstructMatches(carbonyl):
            n_co+=1
        requirement = requirement and n_co<=nmax_co

    # since when okay has to return False:
    return not requirement
ExtraFilters['quinoid'].SetFilterRoutine(findQuinoid)

