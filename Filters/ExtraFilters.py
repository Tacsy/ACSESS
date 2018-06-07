#!/usr/bin/env python
#-*- coding: utf-8 -*-
from filters import NewFilter, NewPatternFilter
from rdkit import Chem
from rdkithelpers import *
from copy import deepcopy
ExtraFilters = dict()

# constants changeable by user input
maxSphericity = 0.3

# EXTRA AROMATICITY FILTER
ExtraFilters['aromatic'] = NewFilter('aromatic')

def HasAromaticity(oldmol):
    mol = deepcopy(oldmol)
    try:
        Chem.SetAromaticity(mol)
    except Exception as e:
        print "in HasAromaticity:", Chem.MolToSmiles(mol), e
        return e
    a = Chem.MolFromSmarts("C1=C-C=C-C=C1")
    b = Chem.MolFromSmarts("c1ccccc1")
    #if mol.HasSubstructMatch(a) or mol.HasSubstructMatch(b):
    if len(mol.GetAromaticAtoms()) or mol.HasSubstructMatch(a):
        return False
    return 'not aromatic'

ExtraFilters['aromatic'].SetFilterRoutine(HasAromaticity)

# Sphericity filters. Other Descriptor 3D filters could be added
# preferably to this particular function to avoid repeatedly
# calculating 3D coordinates
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









##################################################
#  HERE Some specific filters from Jos Teunissen #
#  These will be removed in later versions       #
##################################################




ExtraFilters['quinoid'] = NewFilter('notquinoid')


def findQuinoid(mol, strict=False):

    # quionoid matches
    # the 12 match urges for at least one extra conjugated bond that is not terminal
    ss12 = Chem.MolFromSmarts('[O,o]:,=[#6]~[#6](:,=[O,o])~[#6]=,:*~*')
    ss14 = Chem.MolFromSmarts(
        '[O,o]:,=[#6]-,:[#6,#7]=,:[#6,#7]-,:[#6]:,=[O,o]')
    ss16 = Chem.MolFromSmarts(
        '[O,o]:,=[#6]-,:[#6]=,:[#6,#7]-,:[#6,#7]=,:[#6,#7]-,:[#6]:,=[O,o]')
    ss18 = Chem.MolFromSmarts(
        '[O,o]:,=[#6]-,:[#6]=,:[#6,#7]-,:[#6,#7]=,:[#6,#7]-,:[#6,#7]=,:[#6]-,:[#6]:,=[O,o]'
    )
    sss = [
        ((0, 3), ss12),
        ((0, 5), ss14),
        ((0, 7), ss16),
        ((0, 9), ss18),
    ]

    qmatch = 0
    ssss = [ss12, ss14, ss16, ss18]
    if True:
        matches = set()
        for qinds_p, ss in sss:
            for match in mol.GetSubstructMatches(ss):
                for ind in qinds_p:
                    matches.add(match[ind])
        if len(matches) < 2:
            #print "no match:", Chem.MolToSmiles(mol)
            return True
        elif len(matches) > 2:
            print "multiple matches:", Chem.MolToSmiles(mol)
            return True
        mol.SetProp('quinoid_indices', " ".join(map(str, list(matches))))
    for ss in ssss:
        if mol.HasSubstructMatch(ss): qmatch += 1
    if True:  # molecule needs to have at least conjugated path length of 5 atoms or a conjugated 5/6 ring
        rmatch = 0
        # some extra requirements:
        ring5a = Chem.MolFromSmarts(
            '*1:,=*:,-*:,=*~*1')  #i.e. 5 ring with two internal double bonds
        ring5b = Chem.MolFromSmarts(
            '*1~*(=*)~*=*~*1'
        )  #i.e. 5 ring with double bonds. one internal one external
        ring6a = Chem.MolFromSmarts(
            '*1~*=,:*-,:*=,:*~*1'
        )  #i.e. 6 ring with two internal double bonds in a conjugated pattern
        ring6b = Chem.MolFromSmarts('*1~*=,:*-,:*(:,=*)~*~*1')
        conjp = Chem.MolFromSmarts('*=*-*=*:,-*:,=*')
        extrar = [ring5a, ring5b, ring6a, ring6b, conjp]
        for r in extrar:
            if mol.HasSubstructMatch(r): rmatch += 1
        requirement = qmatch > 0 and rmatch > 0
    else:
        requirement = qmatch > 0
    if True:  #set maximum number of carbonyl groups
        nmax_co = 4
        n_co = 0
        carbonyl = Chem.MolFromSmarts('[#6]=,:[#8;X1]')
        for match in mol.GetSubstructMatches(carbonyl):
            n_co += 1
        requirement = requirement and n_co <= nmax_co

    # since when okay has to return False:
    if aromatic: Chem.SetAromaticity(mol)
    return not requirement


ExtraFilters['quinoid'].SetFilterRoutine(findQuinoid)

############################################
#           MY RADICAL FILTERS            #
############################################
DE = True
ExtraFilters['radical'] = NewFilter('not a radical')


def FindRadical(mol):
    smi_before = Chem.MolToSmiles(mol)
    #print "smi_before:", smi_before

    # 1. Test for radical centers.
    nradcenters = 0
    for atom in mol.GetAtoms():
        num = atom.GetAtomicNum()
        idx = atom.GetIdx()
        if not num == 1.0:  #i.e. we dont want to look for an explicit H
            totalbondorder = atom.GetTotalValence()
            atomcharge = atom.GetFormalCharge()

        # test if not carbene, i.e. carbon with only two bonds and a free electron pair
        if num==6 and totalbondorder==2: return 'carbene in molecule'

        # Test if single radical centers present: also Si/P/S included
        possible_radicals = [(6, 3), (7, 2), (8, 1), (14, 3), (15, 2), (16, 1)]
        for RadAtNum, RadTotalBondorder in possible_radicals:
            if num == RadAtNum and atomcharge == 0:
                #print RadAtNum, RadTotalBondorder, num, totalbondorder
                if totalbondorder == RadTotalBondorder:
                    nradcenters += 1
    #print "nradcenters:", nradcenters, 

    # 2. Test for specific stable radical patterns
    rphenoxyl1a = Chem.MolFromSmarts('a1aaaaa1~[C,c]~[O,o;v1]')
    rphenoxyl1b = Chem.MolFromSmarts('a1aaaaa1~[O,o;v1]')
    rphenoxyl2a = Chem.MolFromSmarts('C1=CC=CC=C1-[O;v1]')
    rphenoxyl2b = Chem.MolFromSmarts('C1=CC=CC=C1-C(=*)-[O;v1]')
    rphenoxyl22 = Chem.MolFromSmarts(
        '[CX3]1~[CX3]~[CX3]~[CX3]~[CX3]~[CX3]1~[CX3]~[O,o;v1]')
    rphenalenyl1 = Chem.MolFromSmarts('a1aaaaa1~[C,c;v3]')
    rphenalenyl2 = Chem.MolFromSmarts('C1=CC=CC=C1[C;v3]')
    rphenoxyl3 = Chem.MolFromSmarts(
        '[O,o;v1]1~[#6,#7]=,:[#6;H0,#7;H0]2~[#6,#7](~[#6,#7]~[#6,#7]~[#6,#7]~[#6,#7]2)~[#6,#7]~[#6,#7]1'
    )
    rphenalenyl3 = Chem.MolFromSmarts(
        '[C,c;v3]1~[#6,#7]=,:[#6;H0,#7;H0]2~[#6,#7](~[#6,#7]~[#6,#7]~[#6,#7]~[#6,#7]2)~[#6,#7]~[#6,#7]1'
    )

    sss = {
        'PHENOXYL1a': rphenoxyl1a,
        'PHENOXYL1b': rphenoxyl1b,
        'PHENOXYL2a': rphenoxyl2a,
        'PHENOXYL2b': rphenoxyl2b,
        'PHENOXYL22': rphenoxyl22,
        'PHENOXYL3': rphenoxyl3,
        'PHENALENYL1': rphenalenyl1,
        'PHENALENYL2': rphenalenyl2,
        'PHENALENYL3': rphenalenyl3
    }
    nmatch = 0
    for name, ss in sss.iteritems():
        for match in mol.GetSubstructMatches(ss):
            nmatch += 1

    # 3. Test if requirements are fullfilled.
    if nradcenters == 0:
        return "no radical center"
    elif nradcenters > 0 and nmatch == 0:
        return "not a valid radical"
    elif nradcenters == 0 and nmatch > 0:
        return "this shouldn't be possible"
    elif nradcenters > 1 and nmatch > 0:
        return "too many radical centers"
    else:
        return False


def FindRadical_old(mol):
    # test presence of max 1 normal radical center:
    rc = Chem.MolFromSmarts('[C;v3;+0]')
    ro = Chem.MolFromSmarts('[O;v1;+0]')
    rn = Chem.MolFromSmarts('[N;v2;+0]')
    sss = [rc, ro, rn]
    once = 0
    for ss in sss:
        oe.OEPrepareSearch(mol, ss)
        for match in ss.Match(mol):
            for ma in match.GetAtoms():
                radcenter = ma.target.GetIdx()
                mol.SetData('radcenter', radcenter)
            once += 1
    if once == 1:  #only one radical center
        return False
    elif once == 0:
        return True
    else:
        return True

ExtraFilters['radical'].SetFilterRoutine(FindRadical)
