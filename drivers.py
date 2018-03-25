#!usr/bin/env python
#-*- coding: utf-8 -*-
'''
 This forms an middle layer between the main algorithm steps and the lowlevel functions
'''
from molfails import *
from rdkithelpers import *
import mutate
import filters
import mprms
import random

###### mutation probabilities:
bondflip=0.8
atomflip=0.8
ringadd=0.1
AddFreq=0.5 #Actual probability: (.8*.7)*.5=.28
DelFreq=0.8 #Actual probability: .224
RingRemove=0.2 #actual probability=.8*.3=.24
MutateStereo=False
StereoFlip=0.2


def DriveMutations(lib):
    nDups, nExcp, nCand=(0,0,0)
    # 1. CROSSOVERS:
    newmols=[]
    for i in xrange( min(mprms.nCross, len(lib)-1) ):
        try:
            if i<mprms.nCross*mprms.EdgeRatio and mprms.EdgeLen>0:
                mol1=random.choice(lib[:mprms.EdgeLen])
                mol2=random.choice(lib+newmols)
            else:
                mol1=random.choice(lib+newmols)
                mol2=random.choice(lib+newmols)
            candidate=mutate.Crossover(mol1,mol2)
            Finalize(candidate)
            newmols.append(candidate)
            nCand+=1
        except MutateFail: nExcp+=1
    nbefore=len(newmols)
    newmols=RemoveDuplicates(newmols)
    nDups+=nbefore-len(newmols)

    # 2. MUTATIONS
    for i in xrange(mprms.nMut):
        if i<mprms.nMut*mprms.EdgeRatio and mprms.EdgeLen>0:
            candidate=Chem.Mol(random.choice(lib[:mprms.EdgeLen]))
        else:
            candidate=Chem.Mol(random.choice(lib+newmols))
        try:
            candidate=MakeMutations(candidate)
            newmols.append(candidate)
            nCand+=1
        except MutateFail: nExcp+=1
    nbefore=len(newmols)
    newmols=RemoveDuplicates(newmols)
    nDups+=nbefore-len(newmols)

    return newmols

def DriveFilters(lib):

    # filter by setting the failed attribute True
    for mol in lib:
        changed, failed = filters.FixAndFilter(mol)
        #if changed: cn.Finalize(mol)
        mol.SetBoolProp('failed', bool(failed))

    # effective filter step. 
    newmols=filter(lambda mol:not mol.GetBoolProp('failed'), lib)
    return newmols

def ExtendPool(pool, lib, newmols):
    return lib+newmols



############ EXTRA FUNCTIONS: #############

##############################################################
# Simple utility to remove duplicate structures
# Tests for identical molecules by generating smiles strings
# Take the molecule with more information generated
def GetScore(mol):
    score=0
    for prop in ['filtered', 'hasstructure', 'Objective', 'tautomerized']:
        try:
            score += bool(mol.GetProp(prop))
        except KeyError:
            pass
    if not mol.HasProp('selected'): mol.SetBoolProp('selected', False)
    return score

def RemoveDuplicates(lib):
    lib.sort( key=lambda x: x.GetProp('isosmi') )
    i=1
    while i<len(lib):
        if lib[i].GetProp('isosmi')==lib[i-1].GetProp('isosmi'):
            iscore=GetScore(lib[i])
            i1score=GetScore(lib[i-1])

            if iscore>=i1score:
                lib[i].SetProp('selected',
                    lib[i].GetProp('selected')+lib[i-1].GetProp('selected') )
                lib.pop(i-1)
            else:
                lib[i-1].SetProp('selected',
                    lib[i].GetProp('selected')+lib[i-1].GetProp('selected') )
                lib.pop(i)

        else: i+=1
    return lib

############ MUTATIONS INTERFACE: ############

# Mutation driver
#@captureMolExceptions
def MakeMutations(candidate):
    candidate=SingleMutate(candidate)
    #candidate,trimmed = TrimAtoms(candidate)
    #if MutateStereo: FlipStereo(candidate)
    try:
        Finalize(candidate)
    except ValueError:
        print Chem.MolToSmiles(candidate)

    return candidate

def SingleMutate(candidateraw):
    #############################################################
    # I still need to test if the candidate is REALLY mutated
    # when actually a RWMol is mutated!. 
    # maybe we have to return the candidate or something?
    #############################################################

    if candidateraw is None: raise ValueError('candidate is none')
    else: candidate = Chem.RWMol(candidateraw)
    #global nAdd,nAddFail,nRemove,nRemoveFail,nBreak,nBreakFail,\
    #       nNewRing,nNewRingFail,nFlip,nFlipFail,nAtomType,nAtomTypeFail,\
    #       nNoMutation
    nFlip,     nFlipFail     =(0,0)
    nAtomType, nAtomTypeFail =(0,0)
    nNewRing,  nNewRingFail  =(0,0)
    nBreak,    nBreakFail    =(0,0)
    nAdd,      nAddFail      =(0,0)
    nRemove,   nRemoveFail   =(0,0)

    parent=candidate.GetProp('isosmi')
    ResetProps(candidate)
    change=False
    candidate.SetProp('parent',parent)

    # 0. Change Backbone
    #If there's a mutatable backbone, call custom routines to mutate it
    #if MutateBackbone:
    #    try:
    #        change=MutateBackbone(candidate)
    #        if change:
    #            cn.Finalize(candidate)
    #            return candidate
    #    except mt.MutateFail: pass

    # 1. bond-flip:
    bonds=list( GetBonds(candidate, notprop='group') )
    if random.random()<bondflip and len(bonds)>0:
        nFlip+=1
        try:
            try:
                mutate.FlipBond(candidate,random.choice(bonds))
                change=True
                Finalize(candidate)
            except: raise MutateFail
        except MutateFail:
            nFlipFail+=1

    # 2. Flip atom identity
    atoms=filter(CanChangeAtom, candidate.GetAtoms())
    if random.random()<atomflip and len(atoms)>0:
        nAtomType+=1
        try:
            mutate.SwitchAtom(candidate,random.choice(atoms))
            change=True
            try: Finalize(candidate)
            except: raise MutateFail
        except MutateFail:
            nAtomTypeFail+=1

    # 3. add a ring - either a new aromatic ring, or bond
    # Aromatic ring addition disabled - probably not necessary
    if random.random()<ringadd:
        nNewRing+=1
        try:
            mutate.AddBond(candidate)
            try: Finalize(candidate)
            except: raise MutateFail
            #return candidate
        except MutateFail:
            nNewRingFail+=1

    # 4. Try to remove a bond to break a cycle
    if random.random() < RingRemove:
        bondringids = candidate.GetRingInfo().BondRings()
        # need to flatten bondringids:
        if bondringids:
            # fancy manner to flatten a list of tuples with possible duplicates
            bondringids = set.intersection(*map(set, bondringids))
            bonds=GetIBonds(bondringids, candidate, notprop='group')
        else: bonds=[]
        nBreak+=1
        if len(bonds) != 0:
            mutate.RemoveBond(candidate,random.choice(bonds))
            try: Finalize(candidate)
            except: raise MutateFail
            #return candidate
        else:
            nBreakFail+=1

    # 5. add an atom
    atoms = GetAtoms(candidate, notprop='protected')
    bonds = GetBonds(candidate, notprop='group'    )
    if ( random.random() < AddFreq and len(atoms)+len(bonds)>0 ):
        nAdd+=1
        try:
            mutate.AddAtom(candidate,random.choice(atoms+bonds))
            try: Finalize(candidate)
            except: raise MutateFail
            #return candidate
        except MutateFail:
            nAddFail+=1

    # 6. remove an atom
    atoms= filter(CanRemoveAtom, candidate.GetAtoms())
    if len(atoms)>1 and random.random() < DelFreq:
        nRemove+=1
        try:
            try:
                mutate.RemoveAtom(candidate,random.choice(atoms))
                Finalize(candidate)
            except: raise MutateFail
            #return candidate
        except MutateFail:
            nRemoveFail+=1

    # Finally: # necessary? # I believe this is not properly implemented
    # so I skip this error raising
    #if not change:
    #    nNoMutation+=1
    #    raise MutateFail()
    if type(candidate)==Chem.RWMol: candidate=candidate.GetMol()
    return candidate




