#!/usr/bin/env python
#-*- coding: utf-8 -*-
debug = False

from filters import NewFilter, NewPatternFilter
from rdkit import Chem
from rdkithelpers import *
from molfails import *
import random

AllFilters = {}

############ FILTERS


def GeomFilter(mol):
    # NotImplemented
    return False


# BREDT VIOLATIONS
AllFilters["Bredt's rule"] = NewFilter("Bredt's rule")
minMacroCycBredt=5
#BredtViolation=Chem.MolFromSmarts('[R&x2]@;=,:[R&x3](@[R&x2])@[R&x2]')
#BredtViolation = NewPatternFilter('bredt violation')
#BredtViolation.SetFilterPattern(Chem.MolFromSmarts('[R]@;=,:[R&x3](@[R])@[R]'))
#BredtViolation.SetExceptions(
#    [Chem.MolFromSmarts('[R]@[R&x3](@[R])@[x3,x4]')])
#    Chem.MolFromSmarts('[R]@;=,:[R&x3]@[R&x3]')])


def OldBredtsRule(mol):
    for match in BredtViolation.FilterWithExceptions(mol):
        macrocycle = False
        atoms = filter(lambda x: x.GetIdx() in match, mol.GetAtoms())
        for atom in atoms:
            SRZ = GetSmallestRingSize(atom)
            if SRZ >= minMacroCycBredt:
                macrocycle = True
                break
        if not macrocycle: return 'Bredt violation'
    return False

def NewBredt(mol):
    doubletworings=Chem.MolFromSmarts('[R]@;=,:[R&x3](@[R])@[R]')
    if mol.HasSubstructMatch(doubletworings):
        # 1. test if match exceptions
        ex1 = Chem.MolFromSmarts('[R]@[R&x3](@[R])@[R&x3,R&x4]')
        ex2 = Chem.MolFromSmarts('[R]@[R&x3](@[R])@[R](@[R])@[R]')
        if any(mol.HasSubstructMatch(ex) for ex in [ex1, ex2]):
            return False
        # 2. test if macrocycle
        for match in mol.GetSubstructMatches(doubletworings):
            macrocycle = False
            atoms = filter(lambda atom:atom.GetIdx() in match, mol.GetAtoms())
            for atom in atoms:
                SRZ = GetSmallestRingSize(atom)
                if SRZ >= minMacroCycBredt:
                    macrocycle = True
                    break
            if not macrocycle: return 'Bredt violation'
    return False


AllFilters["Bredt's rule"].SetFilterRoutine(NewBredt)


def FixBredt(mol):
    changed = False
    #Remove double bonds at bridgeheads
    mtemp = list(BredtViolation.FilterWithExceptions(mol))
    matches = []
    for match in mtemp:
        matches.append(match)
        #for atom in match.GetAtoms():
        for atom in GetIAtoms(match, mol):
            if GetSmallestRingSize(atom) >= 8:
                matches.pop()
                break
    for match in matches:
        for bond in GetAtomIBonds(match, mol):
            bond.SetBondType(Chem.BondType.SINGLE)
        if aromatic: Chem.SetAromaticity(mol)
        changed = True
    return changed


#AllFilters["Bredt's rule"].SetFixRoutine(FixBredt)

####################################################################
# HETEROATOM RATIOS
NandOtoC = 1.0
#NtoC=0.571
NtoC = 0.6  #not what they say in the paper, but it appears to be true in practice
OtoC = 0.666
HalogentoC = 0.5
StoC = 0.333
nitros = Chem.MolFromSmarts('[N+]([O-])=O')
nitriles = Chem.MolFromSmarts('C#N')
sulfones = Chem.MolFromSmarts('O=[S&H0]=O')

AllFilters['atomcounts'] = NewFilter('atomcounts')


def AtomCountFilter(mol):
    nnitros = len(list(mol.GetSubstructMatches(nitros)))
    nnitriles = len(list(mol.GetSubstructMatches(nitriles)))
    nsulfones = len(list(mol.GetSubstructMatches(sulfones)))
    carbon = 1.0 * len(GetXAtoms(mol, 6)) + nsulfones
    nitrogen = 1.0 * len(GetXAtoms(mol, 7)) - nnitriles - nnitros
    oxygen = 1.0 * len(GetXAtoms(mol, 8)) - 2.0 * nsulfones + nnitros
    sulfur = 1.0 * len(GetXAtoms(mol, 16))
    halogen = 1.0 * sum([len(GetXAtoms(mol, num)) for num in [9, 17, 35, 53]])

    if carbon == 0:
        return 'No carbon atoms'
    if nitrogen / carbon > NtoC:
        return 'N/C ratio too high'
    if oxygen / carbon > OtoC:
        return 'O/C ratio too high'
    if (nitrogen + oxygen) / carbon > NandOtoC:
        return '(N+O)/C ratio too high'
    if sulfur / carbon > StoC:
        return 'S/C ratio too high'
    if halogen / carbon > HalogentoC:
        return 'Halogen/C ratio too high'
    return False


AllFilters['atomcounts'].SetFilterRoutine(AtomCountFilter)


def FixAtomQuantities(mol, filtertype=None):
    changed = False
    # Heteroatom ratios
    nnitros = len(list(mol.GetSubstructMatches(nitros)))
    nnitriles = len(list(mol.GetSubstructMatches(nitriles)))
    nsulfones = len(list(mol.GetSubstructMatches(sulfones)))

    nCarb = 1.0 * len(GetXAtoms(mol, 6)) + nsulfones
    if nCarb == 0:  #You really gotta have carbon, sure!
        newcarbons = GetAtoms(mol, notprop='group')
        if len(newcarbons) == 0: raise MutateFail()
        random.choice(newcarbons).SetAtomicNum(6)
        changed = True
        nCarb = 1.0
    nNit = 1.0 * len(GetXAtoms(mol, 7)) - nnitriles - nnitros
    nOxy = 1.0 * len(GetXAtoms(mol, 8)) - 2.0 * nsulfones + nnitros
    nSulf = 1.0 * len(GetXAtoms(mol, 16))
    nHalo = 1.0 * sum([len(GetXAtoms(mol, num)) for num in [9, 17, 35, 53]])

    # These are candidates for changing to satisfy the ratios
    carbon = GetXAtoms(mol, 6, 'group')
    nitrogen = GetXAtoms(mol, 7, 'group')
    oxygen = GetXAtoms(mol, 8, 'group')
    sulfur = GetXAtoms(mol, 16, 'group')
    halogens = [(GetXAtoms(mol, num, 'group')) for num in [9, 7, 35, 53]]
    halogen = [atom for hatoms in halogens for atom in hatoms]

    while nNit / nCarb > NtoC:
        if len(nitrogen) == 0: raise MutateFail()
        topop = random.randrange(len(nitrogen))
        carbon.append(nitrogen.pop(topop))
        carbon[-1].SetAtomicNum(6)
        nCarb += 1
        nNit -= 1
        changed = True

    while nOxy / nCarb > OtoC:
        if len(oxygen) == 0: raise MutateFail()
        topop = random.randrange(len(oxygen))
        carbon.append(oxygen.pop(topop))
        carbon[-1].SetAtomicNum(6)
        nCarb += 1
        nOxy -= 1
        changed = True

    while (nNit + nOxy) / nCarb > NandOtoC:
        if len(oxygen) + len(nitrogen) == 0: raise MutateFail()
        topop = random.randrange(len(oxygen) + len(nitrogen))
        if topop >= len(oxygen):
            nNit -= 1
            topop -= len(oxygen)
            carbon.append(nitrogen.pop(topop))
        else:
            nOxy -= 1
            carbon.append(oxygen.pop(topop))
        nCarb += 1
        carbon[-1].SetAtomicNum(6)
        changed = True

    while nSulf / nCarb > StoC:
        if len(sulfur) == 0: raise MutateFail()
        topop = random.randrange(len(sulfur))
        carbon.append(sulfur.pop(topop))
        carbon[-1].SetAtomicNum(6)
        nSulf -= 1
        nCarb += 1
        changed = True

    while nHalo / nCarb > HalogentoC:
        if len(halogen) == 0: raise MutateFail()
        topop = random.randrange(len(halogen))
        carbon.append(halogen.pop(topop))
        carbon[-1].SetAtomicNum(6)
        nHalo - +1
        nCarb += 1
        changed = True

    return changed


AllFilters['atomcounts'].SetFixRoutine(FixAtomQuantities)

#######################################
# Unsaturations in rings of less than 9 members
AllFilters['triple bond in ring'] = NewFilter('triple bond in ring')


def TripleBondInRing(mol):
    ntripleinlargering = 0
    for bond in mol.GetBonds():
        if bond.IsInRing() and bond.GetBondTypeAsDouble() > 2.0:
            # test if Ringsize smaller than 9:
            if any(bond.IsInRingSize(n) for n in range(2, 10)):
                #if oe.OEBondGetSmallestRingSize(bond)<9 and bond.GetOrder()==3:
                if aromatic: Chem.SetAromaticity(mol)
                return "triple bond in ring"
            else:
                ntripleinlargering += 1
                if ntripleinlargering > 1:
                    return "more than one triple bond in large ring"
                #print "triple bond in large ring:", Chem.MolToSmiles(mol)
    if aromatic: Chem.SetAromaticity(mol)
    return False


AllFilters['triple bond in ring'].SetFilterRoutine(TripleBondInRing)


def FixTripleBondInRing(mol):
    changed = False
    for bond in mol.GetBonds():
        if bond.IsInRing() and bond.GetBondTypeAsDouble() > 2.0:
            if any(bond.IsInRingSize(n) for n in range(2, 10)):
                bond.SetBondType(Chem.BondType.SINGLE)
                changed = True
    return changed


AllFilters['triple bond in ring'].SetFixRoutine(FixTripleBondInRing)

##################################################
# ROUTINES TO FIX GENERAL FILTER VIOLATIONS      #
##################################################


########################################
# Groups that can be fixed by ... removing heteroatoms
# Assumes that the first heteroatom found is
# the one that we want to remove
def FixByRemovingHeteroatoms(mol, filter):
    changed = False
    Chem.SetAromaticity(mol)
    matches = filter.FilterWithExceptions(mol)
    Chem.Kekulize(mol, True)
    while len(matches) > 0:
        match = matches[0]
        fixd = False
        for atom in GetIAtoms(match, mol):
            if (atom.GetAtomicNum() not in (1,6)) and \
                   CanRemoveAtom(atom):
                if atom.HasProp('grouprep'):
                    mutate.RemoveGroup(mol, atom.GetProp('group'))
                atom.SetAtomicNum(6)
                changed = True
                Chem.SetAromaticity(mol)
                matches = filter.FilterWithExceptions(mol)
                Chem.Kekulize(mol, True)
                fixd = True
                break
        if not fixd: raise MutateFail()
    return changed


########################################
# Groups that can be fixed by ... eliminating unsaturated bonds
def FixBySaturating(mol, filter):
    changed = False
    matches = filter.FilterWithExceptions(mol)
    while len(matches) > 0:
        match = matches[0]
        fixd = False
        for bond in GetAtomIBonds(match, mol):
            if not bond.HasProp('group') and bond.GetBondTypeAsDouble() > 1.0:
                changed = True
                fixd = True
                bond.SetBondType(Chem.BondType.SINGLE)
                matches = filter.FilterWithExceptions(mol)
                break  #just need to modify the first bond
        if not fixd: raise MutateFail()
    return changed


########################################
# Groups that can be fixed by ... deleting specific bonds
def RemoveBondNumber(iBond):
    def fixer(mol, filter):
        changed = False
        matches = filter.FilterWithExceptions(mol)
        while len(matches) > 0:
            bonds = list(GetAtomIBonds(matches[0], mol))
            if not bonds[iBond].HasProp('group'):
                #mol.RemoveBond( bonds[iBond].GetBeginAtomIdx(), bonds[iBond].GetEndAtomIdx() )
                break
                changed = True
            else:
                break
            matches = filter.FilterWithExceptions(mol)
        return changed

    return fixer


#################################################################
# PATTERN MATCHING FILTERS
#################################################################

newfilt = NewPatternFilter('allene')
newfilt.SetFilterPattern(Chem.MolFromSmarts("[!O]=*=*"))
AllFilters['allene'] = newfilt

newfilt = NewPatternFilter('acidtaut')
newfilt.SetFilterPattern(Chem.MolFromSmarts("C=CO(O)"))
AllFilters['acidtaut'] = newfilt

newfilt = NewPatternFilter('aminal')
newfilt.SetFilterPattern(Chem.MolFromSmarts("[N&X3][C&X4][N,O]"))
AllFilters['aminal'] = newfilt

newfilt = NewPatternFilter('C=N')
newfilt.SetFilterPattern(
    Chem.MolFromSmarts("[C&X3]=[N&X2]"))  # ca veut dire: sp2C=sp2N
AllFilters['C=N'] = newfilt
newfilt.SetExceptions([
    Chem.MolFromSmarts("[#7]C=N"),
    Chem.MolFromSmarts("[C&X3]=N[N,O,n,o]"),
    Chem.MolFromSmarts("cn")
])
#if allow_imines:
#    newfilt.SetExceptions([Chem.MolFromSmarts("[#7]C=N"),
#                           Chem.MolFromSmarts("[C&X3]=N[N,O,n,o]"),
#                           Chem.MolFromSmarts('[#6]=[#6]-[#6]=N'), # C=N conjugated with another C=C
#                           Chem.MolFromSmarts("cn")])
#else:
#    newfilt.SetExceptions([Chem.MolFromSmarts("[#7]C=N"),
#                           Chem.MolFromSmarts("[C&X3]=N[N,O,n,o]"),
#                           Chem.MolFromSmarts("cn")])

newfilt = NewPatternFilter('Decarboxy1')
newfilt.SetFilterPattern(Chem.MolFromSmarts("CC(=O)[C&X4]C([O&H1])=O"))
AllFilters['Decarboxy1'] = newfilt
newfilt = NewPatternFilter('Decarboxy2')
newfilt.SetFilterPattern(Chem.MolFromSmarts("O=[C&H1][C&X4]C([O&H1])=O"))
AllFilters['Decarboxy2'] = newfilt
#Need bridgehead detection for exceptions

newfilt = NewPatternFilter('Enamine')
newfilt.SetFilterPattern(Chem.MolFromSmarts("[C&X3]=C[N&X3]"))
AllFilters['Enamine'] = newfilt
newfilt.SetExceptions([
    Chem.MolFromSmarts("[C&X3]=CNC=[O,N]"),
    Chem.MolFromSmarts("[O,N]=CC=C[N&X3]"),
    Chem.MolFromSmarts("[c&X3]c[n,N]"),
    Chem.MolFromSmarts("[C,N,O]1C([C,O])=C([C,O])NC([C,O])=C1([C,O])"),  #0000
    Chem.MolFromSmarts(
        "[C,N,O]1[CH1]=C([C,O])NC([C,O])=C1([C,O])"),  #1000,0001
    Chem.MolFromSmarts("[C,N,O]1[CH1]=[CH1]NC([C,O])=C1([C,O])"),  #1100,0011
    Chem.MolFromSmarts("[C,N,O]1[CH1]=C([C,O])NC([C,O])=[CH1]1"),  #1001
    Chem.MolFromSmarts("[C,N,O]1[CH1]=[CH1]NC([C,O])=[CH1]1"),  #1101,1011
    Chem.MolFromSmarts("[C,N,O]1[CH1]=[CH1]NC[CH1]=[CH1]1"),  #1111
    Chem.MolFromSmarts("[C,N,O]1[CH1]=[CH1]NC[CH1]=C1([C,O])"),  #1110,0111
    Chem.MolFromSmarts("[C,N,O]1[CH1]=C([C,O])N[CH1]=C1([C,O])"),  #1010,0101
    Chem.MolFromSmarts(
        "[C,N,O]1C([C,O])=[CH1]NC([C,O])=C1([C,O])"),  #0100,0010
    Chem.MolFromSmarts("[C,N,O]1C([C,O])=[CH1]N[CH1]=C1([C,O])"),  #0110
    Chem.MolFromSmarts("[C,N,O]1C=CNC=C1"),
    Chem.MolFromSmarts('N#CC=CN'),
    Chem.MolFromSmarts('C=CNS(=O)=O'),
    Chem.MolFromSmarts('NC=CS(=O)=O'),
    Chem.MolFromSmarts('C=CNC=S')
])

newfilt = NewPatternFilter('Enol')
newfilt.SetFilterPattern(Chem.MolFromSmarts("*=C[O&H1]"))
AllFilters['Enol'] = newfilt
newfilt.SetExceptions([
    Chem.MolFromSmarts("O=C[O&H1]"),
    Chem.MolFromSmarts("[a]=[c,C][O&H1]"),
    Chem.MolFromSmarts("[a]C=C[O&H1]")
])

newfilt = NewPatternFilter('FC1')
newfilt.SetFilterPattern(Chem.MolFromSmarts("F[C&X4][N,O]"))
AllFilters['FC1'] = newfilt
newfilt.SetExceptions(
    [Chem.MolFromSmarts("Fc([a])=[N,O]")])  #Should [n,o] be aromatic or not?

newfilt = NewPatternFilter('FC2')
newfilt.SetFilterPattern(Chem.MolFromSmarts("FC=[N,O]"))
AllFilters['FC2'] = newfilt
newfilt.SetExceptions(
    [Chem.MolFromSmarts("Fc([a])=[N,O]")])  #Should [n,o] be aromatic or not?

newfilt = NewPatternFilter('Geminal1')
newfilt.SetFilterPattern(Chem.MolFromSmarts("[N&H2,O&H1][C&X4][N&H2,O&H1]"))
AllFilters['Geminal1'] = newfilt

newfilt = NewPatternFilter('Geminal2')
newfilt.SetFilterPattern(Chem.MolFromSmarts("[C&X4][N&H1][C&X4][N&H1][C&X4]"))
AllFilters['Geminal2'] = newfilt

newfilt = NewPatternFilter('Hemiacetal')
newfilt.SetFilterPattern(Chem.MolFromSmarts("[O&H1][C&X4][N,O]"))
AllFilters['Hemiacetal'] = newfilt

newfilt = NewPatternFilter('Hemiaminal1')
newfilt.SetFilterPattern(Chem.MolFromSmarts("[N&H2][C&X4]O"))
AllFilters['Hemiaminal1'] = newfilt

newfilt = NewPatternFilter('Hemiaminal2')
newfilt.SetFilterPattern(Chem.MolFromSmarts("C[N&H1][C&X4]O"))
AllFilters['Hemiaminal2'] = newfilt

newfilt = NewPatternFilter('Hemiaminal3')
newfilt.SetFilterPattern(Chem.MolFromSmarts("CN(C)[C&X4]O"))
AllFilters['Hemiaminal3'] = newfilt

newfilt = NewPatternFilter('HetHet_OO')
newfilt.SetFilterPattern(Chem.MolFromSmarts("[#8][#8]"))
AllFilters['HetHet_OO'] = newfilt

#From HetHet_NO
newfilt = NewPatternFilter('HetHet_NS')
newfilt.SetFilterPattern(Chem.MolFromSmarts("NSO"))
AllFilters['HetHet_NS'] = newfilt
newfilt.SetExceptions([
    Chem.MolFromSmarts("C=NS[C,S]"),
    Chem.MolFromSmarts("C=N[S&H1]"),
    Chem.MolFromSmarts("[N,O]=CNSC"),
    Chem.MolFromSmarts("[N,O,S]=CN[S&H1]"), nitros
])

newfilt = NewPatternFilter('HetHet_OS1')
newfilt.SetFilterPattern(Chem.MolFromSmarts("C=COSN"))
newfilt.SetExceptions([
    Chem.MolFromSmarts("S1(=O)O[C,S](=O)**1"),
    Chem.MolFromSmarts("S1(=O)O[C,S](=O)***1"),
])
AllFilters['HetHet_OS1'] = newfilt

newfilt = NewPatternFilter('HetHet_OS2')
newfilt.SetFilterPattern(Chem.MolFromSmarts("NS[O&H1]"))
newfilt.SetExceptions([
    Chem.MolFromSmarts("S1(=O)O[C,S](=O)**1"),
    Chem.MolFromSmarts("S1(=O)O[C,S](=O)***1"),
])
AllFilters['HetHet_OS2'] = newfilt

newfilt = NewPatternFilter('HetHet_OS3')
newfilt.SetFilterPattern(Chem.MolFromSmarts("[O&H1]SO[C,O,N]"))
newfilt.SetExceptions([
    Chem.MolFromSmarts("S1(=O)O[C,S](=O)**1"),
    Chem.MolFromSmarts("S1(=O)O[C,S](=O)***1"),
])
AllFilters['HetHet_OS3'] = newfilt

#This filter isn't explicitly stated but is clearly used.
newfilt = NewPatternFilter('Hetero triple bond')
newfilt.SetFilterPattern(Chem.MolFromSmarts("[N,O]C#*"))
newfilt.SetFixRoutine(FixByRemovingHeteroatoms)
AllFilters["Hetero triple bond"] = newfilt

newfilt = NewPatternFilter('HetHet_NN')
newfilt.SetFilterPattern(Chem.MolFromSmarts("N~N"))
AllFilters['HetHet_NN'] = newfilt
newfilt.SetExceptions([
    Chem.MolFromSmarts("[#6]N([#6])N=[#6]"),
    Chem.MolFromSmarts("[#6]=N[N&H1][#6]"),
    Chem.MolFromSmarts("[c,C]=N[N&H2]"),
    Chem.MolFromSmarts("[#6]=NN=[#6]"),
    Chem.MolFromSmarts("NNC=[N,O,S]"),
    Chem.MolFromSmarts('NNS=O')
])

newfilt = NewPatternFilter('HetHet_XO')
newfilt.SetFilterPattern(Chem.MolFromSmarts("[F,Cl,Br,I]O"))
AllFilters['HetHet_XO'] = newfilt

newfilt = NewPatternFilter('HetHet_XN')
newfilt.SetFilterPattern(Chem.MolFromSmarts("[F,Cl,Br,I]N"))
AllFilters['HetHet_XN'] = newfilt

newfilt = NewPatternFilter('HetHet_NO')
newfilt.SetFilterPattern(Chem.MolFromSmarts("N-,=O"))
AllFilters['HetHet_NO'] = newfilt
newfilt.SetExceptions([
    Chem.MolFromSmarts("[#6]=NO[C,S]"),
    Chem.MolFromSmarts("[#6]=N[O&H1]"),
    Chem.MolFromSmarts("[N,O,S]=[C,S]NO[C,S]"),
    Chem.MolFromSmarts("[N,O,S]=CN[O&H]"), nitros
])

newfilt = NewPatternFilter('HetHet_SS')
newfilt.SetFilterPattern(Chem.MolFromSmarts("[#16][#16]"))
AllFilters['HetHet_SS'] = newfilt

newfilt = NewPatternFilter('HetHet_Aromatic')
newfilt.SetFilterPattern(Chem.MolFromSmarts("[n,o][n,o][n,o][n,o][n,o]"))
AllFilters['HetHet_Aromatic'] = newfilt

newfilt = NewPatternFilter('HetHetHet')
newfilt.SetFilterPattern(Chem.MolFromSmarts("N[O,N]N"))
AllFilters['HetHetHet'] = newfilt

newfilt = NewPatternFilter('Hetero-sulfur without sulfone')
newfilt.SetFilterPattern(Chem.MolFromSmarts('*[S&X2][!#6]'))
AllFilters['Hetero-sulfur without sulfone'] = newfilt

newfilt = NewPatternFilter('Hetero-SR1')
newfilt.SetFilterPattern(Chem.MolFromSmarts("[O,N]1C[O,N]1"))
AllFilters['Hetero-SR1'] = newfilt

newfilt = NewPatternFilter('Hetero-SR2')
newfilt.SetFilterPattern(Chem.MolFromSmarts("C1[N,O][N,O]C1"))
AllFilters['Hetero-SR2'] = newfilt

newfilt = NewPatternFilter('Hetero-SR3')
newfilt.SetFilterPattern(Chem.MolFromSmarts("C1[N,O]C[N,O]1"))
AllFilters['Hetero-SR3'] = newfilt

newfilt = NewPatternFilter('Hetero-SR4')
newfilt.SetFilterPattern(Chem.MolFromSmarts("[C,N]=*1***1"))
AllFilters['Hetero-SR4'] = newfilt

newfilt = NewPatternFilter('Hetero-SR5')
newfilt.SetFilterPattern(Chem.MolFromSmarts("*=*1**1"))
AllFilters['Hetero-SR5'] = newfilt

newfilt = NewPatternFilter('Hetero-SR6')
newfilt.SetFilterPattern(Chem.MolFromSmarts("O=C1**C1=O"))
AllFilters['Hetero-SR6'] = newfilt

newfilt = NewPatternFilter('Hetero-SR7')
newfilt.SetFilterPattern(Chem.MolFromSmarts("O=[C,S]1*[C,S](=O)*1"))
AllFilters['Hetero-SR7'] = newfilt

#This one is a bit of a logical statement, so we do it special-like
# it's specifically:
# (IntraMolOr1 .or. IntraMolOr2) .and. (IntraMolAnd1 .or. IntraMolAnd2)
IntraMolOr1 = Chem.MolFromSmarts('[C&X4][N&H2]')  #i.e. sp3-C - NH2
IntraMolOr2 = Chem.MolFromSmarts('[N&X3][N&H2]')  #i.e. sp3-N - NH2
IntraMolAnd1 = Chem.MolFromSmarts('O=[S,C](C)C')  #i.e. keton
IntraMolAnd2 = Chem.MolFromSmarts('O=[S&H1,C&H1]C')  #
newfilt = NewFilter('Intramol')


def findIntramol(mol):
    if mol.HasSubstructMatch(IntraMolOr1)  or mol.HasSubstructMatch(IntraMolOr2)  and \
       mol.HasSubstructMatch(IntraMolAnd1) or mol.HasSubstructMatch(IntraMolAnd2):
        return True
    else:
        return False


intraOr = (IntraMolOr1, IntraMolOr2)
intraAnd = (IntraMolAnd1, IntraMolAnd2)
intraFindTemp = NewPatternFilter('')


def fixIntramol(mol):
    changed = False

    if random.random() < 0.5: mypatterns = intraOr
    else: mypatterns = intraAnd

    for pattern in mypatterns:
        if mol.HasSubstructMatch:
            intraFindTemp.SetFilterPattern(pattern)
            madeChange = FixByRemovingHeteroatoms(mol, intraFindTemp)
            changed = changed or madeChange
    return changed


newfilt.SetFixRoutine(fixIntramol)
newfilt.SetFilterRoutine(findIntramol)
AllFilters['Intramol'] = newfilt

newfilt = NewPatternFilter('Mixed1')
newfilt.SetFilterPattern(Chem.MolFromSmarts("C=CON"))
AllFilters['Mixed1'] = newfilt

newfilt = NewPatternFilter('Mixed2')
newfilt.SetFilterPattern(Chem.MolFromSmarts("C=COC(N)=[O,N]"))
AllFilters['Mixed2'] = newfilt

newfilt = NewPatternFilter('Mixed3')
newfilt.SetFilterPattern(Chem.MolFromSmarts('C=COC(O)=N'))
AllFilters['Mixed3'] = newfilt

newfilt = NewPatternFilter('Mixed4')
newfilt.SetFilterPattern(Chem.MolFromSmarts('O=C(N)[O&H1]'))
AllFilters['Mixed4'] = newfilt

newfilt = NewPatternFilter('Mixed5')
newfilt.SetFilterPattern(Chem.MolFromSmarts('O=C(N)[O&H1][C,O,N]'))
AllFilters['Mixed5'] = newfilt

newfilt = NewPatternFilter('Mixed6')
newfilt.SetFilterPattern(Chem.MolFromSmarts('[C,S](=O)O[C,S]=O'))
AllFilters['Mixed6'] = newfilt
newfilt.SetExceptions([
    Chem.MolFromSmarts("O=[C,S]1**[C,S](O1)=O"),
    Chem.MolFromSmarts("O=[C,S]1***[C,S](O1)=O")
])

newfilt = NewPatternFilter('Mixed7')
newfilt.SetFilterPattern(Chem.MolFromSmarts('OCOCO'))
AllFilters['Mixed7'] = newfilt

newfilt = NewPatternFilter(
    'Ortho1')  #Note that this differs from what they said
newfilt.SetFilterPattern(Chem.MolFromSmarts('[N,O]C(C)(C)N'))
AllFilters['Ortho1'] = newfilt

newfilt = NewPatternFilter('Ortho2')
newfilt.SetFilterPattern(Chem.MolFromSmarts('[N,O][C&H1](C)N'))
AllFilters['Ortho2'] = newfilt

newfilt = NewPatternFilter('Orthoester')
newfilt.SetFilterPattern(Chem.MolFromSmarts('C([N,O])([N,O])([N,O])'))
AllFilters['Orthoester'] = newfilt

newfilt = NewPatternFilter('Topo1-33spiro')
newfilt.SetFilterPattern(Chem.MolFromSmarts('*~1~*~*~1~2~*~*2'))
AllFilters['Topo1-33spiro'] = newfilt

newfilt = NewPatternFilter('Topo1-33fuse')
newfilt.SetFilterPattern(Chem.MolFromSmarts('*~1~*~2~*1~*2'))
AllFilters['Topo1-33fuse'] = newfilt

newfilt = NewPatternFilter('Topo1-34spiro')
newfilt.SetFilterPattern(Chem.MolFromSmarts('*~1~*~*~*~1~2~C~*2'))
AllFilters['Topo1-34spiro'] = newfilt

newfilt = NewPatternFilter('Topo1-34fuse')
newfilt.SetFilterPattern(Chem.MolFromSmarts('*~1~*~*~2~*1~*2'))
AllFilters['Topo1-34fuse'] = newfilt

newfilt = NewPatternFilter('Topo1-44spiro')
newfilt.SetFilterPattern(Chem.MolFromSmarts('*~1~*~*~*~1~2~*~*~*2'))
AllFilters['Topo1-44spiro'] = newfilt

newfilt = NewPatternFilter('Topo1-44fuse')
newfilt.SetFilterPattern(Chem.MolFromSmarts('*~1~*~*~2*~1~*~*2'))
AllFilters['Topo1-44fuse'] = newfilt

newfilt = NewPatternFilter('Topo1-44bridge')
newfilt.SetFilterPattern(Chem.MolFromSmarts('*~1~*~2~*~*~1~*~2'))
AllFilters['Topo1-44bridge'] = newfilt

newfilt = NewPatternFilter('Non-aromatic nitro')
newfilt.SetFilterPattern(Chem.MolFromSmarts('[!a][N+](=O)[O-]'))
AllFilters['Non-aromatic nitro'] = newfilt

newfilt = NewPatternFilter('Non-aromatic halogen')
newfilt.SetFilterPattern(Chem.MolFromSmarts('[!c][F,Cl,Br,I]'))
AllFilters['Non-aromatic halogen'] = newfilt

newfilt = NewPatternFilter('Hetero-aromatic halogen')
newfilt.SetFilterPattern(Chem.MolFromSmarts('[!c]c[F,Cl,Br,I]'))
AllFilters['Hetero-aromatic halogen'] = newfilt

newfilt = NewPatternFilter('Double bond in 4 ring')
newfilt.SetFilterPattern(Chem.MolFromSmarts('*~1~*=,#,:*~*1'))
AllFilters['Double bond in 4 ring'] = newfilt
newfilt.SetExceptions([Chem.MolFromSmarts('[A]1[A][a]~[a]1')])

newfilt = NewPatternFilter('Double bond in 3 ring')
newfilt.SetFilterPattern(Chem.MolFromSmarts('*~1~*=,#,:*1'))
AllFilters['Double bond in 3 ring'] = newfilt

newfilt = NewPatternFilter('Polycyclic1')
newfilt.SetFilterPattern(Chem.MolFromSmarts('*@N(@*)@[O,N]'))
AllFilters['Polycyclic1'] = newfilt

newfilt = NewPatternFilter('Polycyclic2')
newfilt.SetFilterPattern(Chem.MolFromSmarts('*@N(@*)@C@=*'))
AllFilters['Polycyclic2'] = newfilt

newfilt = NewPatternFilter('Beta keto carboxyl')
newfilt.SetFilterPattern(Chem.MolFromSmarts('[O&H1]C(=O)[C&X4]C(=O)[!N]'))
AllFilters['Beta keto carboxylate'] = newfilt

############################################3
#Fixing routines for various filters
for filtname in [
        'aminal', 'C=N', 'Enamine', 'HetHet_NN', 'HetHet_OO', 'HetHet_NO',
        'Ortho1', 'Ortho2', 'HetHet_SS', 'Enol', 'HetHet_NS', 'HetHet_OS1',
        'HetHet_OS3', 'HetHet_OS2', 'HetHet_Aromatic', 'HetHetHet',
        'Non-aromatic nitro', 'Non-aromatic halogen', 'Beta keto carboxylate',
        'Orthoester'
]:
    AllFilters[filtname].SetFixRoutine(FixByRemovingHeteroatoms)

for filtname in [
        'Hetero-SR4', 'Hetero-SR5', 'Hetero-SR6', 'Hetero-SR7',
        'Double bond in 3 ring', 'Double bond in 4 ring', 'allene', 'acidtaut'
]:
    AllFilters[filtname].SetFixRoutine(FixBySaturating)


def FixTopo144Bridge(mol, filter):
    changed = False
    matches = filter.FilterWithExceptions(mol)
    while len(matches) > 0:
        atoms = list(GetIAtoms(matches[0], mol))
        bonds = list(GetAtomIBonds(matches[0], mol))
        if not (bonds[0].HasProp('group') or bonds[2].HasProp('group')):
            mol.RemoveBond(bonds[0].GetBeginAtomIdx(),
                           bonds[0].GetEndAtomIdx())
            mol.RemoveBond(bonds[2].GetBeginAtomIdx(),
                           bonds[2].GetEndAtomIdx())
            try:
                mol.AddBond(atoms[1].GetIdx(), atoms[4].GetIdx(), bondorder[1])
            except RuntimeError as e:
                print e
                break
            changed = True
        else:
            break
        matches = filter.FilterWithExceptions(mol)
    return changed


AllFilters['Topo1-44bridge'].SetFixRoutine(FixTopo144Bridge)
AllFilters['Topo1-33fuse'].SetFixRoutine(RemoveBondNumber(3))
AllFilters['Topo1-34fuse'].SetFixRoutine(RemoveBondNumber(4))
AllFilters['Topo1-44fuse'].SetFixRoutine(RemoveBondNumber(4))
