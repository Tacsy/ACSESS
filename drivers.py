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
from output import logtime

###### mutation probabilities:
p_BondFlip=0.8
p_AtomFlip=0.8
p_RingAdd=0.1
p_AddFreq=0.5 #Actual probability: (.8*.7)*.5=.28
p_DelFreq=0.8 #Actual probability: .224
p_RingRemove=0.2 #actual probability=.8*.3=.24
MutateStereo=False
StereoFlip=0.2

##### workflow switches:
GenStruc      = 0
startFilter   = 2
startTautomer = 10
startGenStruc = 20
KeepNoGeomPool= True

debug=False

def SetIterationWorkflow(gen):
    global GenStruc
    Tautomerizing = (gen>=startTautomer)
    #if Tautomerizing: cn.CanonicalTautomer=True

    Filter = (gen>= startFilter)
    GenStruc = (GenStruc and gen>= startGenStruc)
    return Tautomerizing, Filter, GenStruc


@logtime()
def DriveMutations(lib):
    nDups, nExcp, nCand=(0,0,0)
    # 1. CROSSOVERS:
    print "crossovers...",
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
            if type(candidate)==Chem.RWMol: candidate=candidate.GetMol()
            Finalize(candidate)
            newmols.append(candidate)
            nCand+=1
        except MutateFail: nExcp+=1
    newmols=filter(Sane, newmols)

    nbefore=len(newmols)
    newmols=RemoveDuplicates(newmols)
    nDups+=nbefore-len(newmols)
    if debug: print "nDups after CX:", nDups

    # 2. MUTATIONS
    print "mutations...",
    for i in xrange(mprms.nMut):
        if i<mprms.nMut*mprms.EdgeRatio and mprms.EdgeLen>0:
            candidate=Chem.Mol(random.choice(lib[:mprms.EdgeLen]))
        else:
            candidate=Chem.Mol(random.choice(lib+newmols))
        try:
            candidate=MakeMutations(candidate)
            if type(candidate)==Chem.RWMol:
                candidate=candidate.GetMol()
            try:
                Finalize(candidate)
            except Exception as e:
                raise MutateFail(candidate)
            newmols.append(candidate)
            nCand+=1
        except MutateFail: nExcp+=1
    newmols=filter(Sane, newmols)
    nbefore=len(newmols)
    newmols=RemoveDuplicates(newmols)
    nDups=nbefore-len(newmols)
    if debug: print "nDups after MUs:", nDups
    if debug:
        with open('mutatelib','w') as f:
            for mol in newmols: f.write(Chem.MolToSmiles(mol)+'\n')

    return newmols

@logtime()
def DriveFilters(lib, Filter=True, GenStrucs=False):
    if not (Filter or GenStrucs): return lib
    print "filtering...", 

    # filter by setting the failed attribute True
    for mol in lib:
        changed, failed = filters.FixAndFilter(mol)
        if debug and failed:
            print "changed:{}, failed:{}, mol:{}".format(changed, failed, Chem.MolToSmiles(mol))
        mol.SetBoolProp('failed', bool(failed))
        mol.SetBoolProp('filtered', True)
        if failed in ['unknown', 'False'] or failed==True:
            failed='unknown'
        if debug and failed:
            mol.SetProp('failedfilter', failed)
            print failed
            filterFile.write(Chem.MolToSmiles(mol) + '  ' + failed + '\n')
    # effective filter step. 
    newmols=filter(lambda mol:not mol.GetBoolProp('failed'), lib)
    if debug: print "    nFiltered:", len(lib)-len(newmols)
    filterFile.flush()
    return newmols

def DrivePoolFilters(pool, Filtering, GenStrucs, Tautomerizing, gen):
    global GenStruc
    pool = filter(bool, pool)
    if Tautomerizing:
        CanonicalTautomer=True
        if gen==startTautomer:
            print 'tautomerizing pool...',
            pool = filter(Tautomerize, pool)
    if GenStrucs and gen==startGenStrucs and not KeepNoGeomPool:
        print 'restarting and filtering pool...',
        pool=[mol for mol in pool if mol.GetBoolProp('hasstructure')]
        DriveFilters(pool, Filtering, GenStruc)
    if (Filtering and gen==startFilter) or (
        GenStruc and gen==startGenStruc and KeepNoGeomPool):
        print 'pool-',
        pool = DriveFilters(pool, Filtering, GenStrucs)
    return pool


@logtime()
def ExtendPool(pool, lib, newmols):
    pool=RemoveDuplicates(lib+newmols)
    return pool

@logtime()
def DriveSelection(pool, subsetSize):
    print "selecting...",
    #1. select maximin algorithm.
    if hasattr(mprms,'metric') :
        from distance import Maximin
    else:
        from similarity import Maximin

    #2. select
    answer = Maximin(pool, mprms.subsetSize)
    if type(answer)==tuple:
        siml=answer[0]
        lib=answer[1]
    else:
        lib=answer

    return lib



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
    if not mol.HasProp('selected'): mol.SetIntProp('selected', 0)
    return score

def RemoveDuplicates(lib):
    lib.sort( key=lambda x: x.GetProp('isosmi') )
    i=1
    while i<len(lib):
        #print len(lib),
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
    if random.random()<p_BondFlip and len(bonds)>0:
        if debug: print "1",
        nFlip+=1
        try:
            Chem.Kekulize(candidate, True)
            bonds=list( GetBonds(candidate, notprop='group') )
            mutate.FlipBond(candidate,random.choice(bonds))
            Finalize(candidate)
            if aromatic: Chem.SetAromaticity(mol)
            return candidate
        except MutateFail:
            if aromatic: Chem.SetAromaticity(mol)
            nFlipFail+=1

    # 2. Flip atom identity
    atoms=filter(CanChangeAtom, candidate.GetAtoms())
    if random.random()<p_AtomFlip and len(atoms)>0:
        if debug: print "2",
        nAtomType+=1
        try:
            mutate.SwitchAtom(candidate,random.choice(atoms))
            try: 
                Finalize(candidate)
            except:
                raise MutateFail
            return candidate
        except MutateFail:
            nAtomTypeFail+=1

    # 3. add a ring - either a new aromatic ring, or bond
    # Aromatic ring addition disabled - probably not necessary
    if random.random()<p_RingAdd:
        if debug: print "3",
        nNewRing+=1
        try:
            mutate.AddBond(candidate)
            try: Finalize(candidate)
            except: raise MutateFail
            return candidate
        except MutateFail:
            nNewRingFail+=1

    # 4. Try to remove a bond to break a cycle
    if random.random() < p_RingRemove:
        if debug: print "4",
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
            except Exception as e:
                print "error with RemoveBond with mol:", Chem.MolToSmiles(True)
                print e
                raise MutateFail
            return candidate
        else:
            nBreakFail+=1

    # 5. add an atom
    atoms = GetAtoms(candidate, notprop='protected')
    bonds = GetBonds(candidate, notprop='group'    )
    if ( random.random() < p_AddFreq and len(atoms)+len(bonds)>0 ):
        if debug: print "5",
        nAdd+=1
        try:
            mutate.AddAtom(candidate,random.choice(atoms+bonds))
            try: Finalize(candidate)
            except: raise MutateFail
            return candidate
        except MutateFail:
            nAddFail+=1

    # 6. remove an atom
    if len(atoms)>1 and random.random() < p_DelFreq:
        if debug: print "6",
        nRemove+=1
        Chem.Kekulize(candidate, True)
        atoms= filter(CanRemoveAtom, candidate.GetAtoms())
        try:
            try:
                mutate.RemoveAtom(candidate, random.choice(atoms))
                Finalize(candidate)
            except Exception as e:
                print "e:", e, "mol:", Chem.MolToSmiles(candidate) 
                raise MutateFail
            if aromatic: Chem.SetAromaticity(mol)
            return candidate
        except MutateFail:
            nRemoveFail+=1
        if aromatic: Chem.SetAromaticity(mol)

    # Finally: # necessary? # I believe this is not properly implemented
    # so I skip this error raising
    #if not change:
    #    nNoMutation+=1
    #    raise MutateFail()
    if Chem.MolToSmiles(candidate)=='O=c1[nH][nH]ccc1=S': print "HERE!"
    #if Chem.MolToSmiles(candidate)=='ccc(cc(c)CC)COC': print "HERE!b"
    if type(candidate)==Chem.RWMol: candidate=candidate.GetMol()
    return candidate



