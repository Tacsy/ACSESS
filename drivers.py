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

    Finalize(candidate)

    return candidate

def SingleMutate(candidate):
    #global nAdd,nAddFail,nRemove,nRemoveFail,nBreak,nBreakFail,\
    #       nNewRing,nNewRingFail,nFlip,nFlipFail,nAtomType,nAtomTypeFail,\
    #       nNoMutation
    nFlip,     nFlipFail     =(0,0)
    nAtomType, nAtomTypeFail =(0,0)

    parent=candidate.GetProp('isosmi')
    ResetProps(candidate)
    change=False
    candidate.SetProp('parent',parent)

    #If there's a mutatable backbone, call custom routines to mutate it
    #if MutateBackbone:
    #    try:
    #        change=MutateBackbone(candidate)
    #        if change:
    #            cn.Finalize(candidate)
    #            return candidate
    #    except mt.MutateFail: pass

    #Try a bond-flip:
    bonds=list( GetBonds(candidate, notprop='group') )
    if random.random()<bondflip and len(bonds)>0:
        nFlip+=1
        try:
            mutate.FlipBond(candidate,random.choice(bonds))
            change=True
            Finalize(candidate)
        except MutateFail:
            nFlipFail+=1

    #Flip atom identity
    atoms=filter(CanChangeAtom, candidate.GetAtoms())
    if random.random()<atomflip and len(atoms)>0:
        nAtomType+=1
        try:
            mutate.SwitchAtom(candidate,random.choice(atoms))
            change=True
            Finalize(candidate)
        except MutateFail as mmf:
            nAtomTypeFail+=1

    return candidate




