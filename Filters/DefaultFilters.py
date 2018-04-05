#!/usr/bin/env python
#-*- coding: utf-8 -*-

from filters import NewFilter, NewPatternFilter
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkithelpers import *
from molfails import *
import random
try:
    from GraphPlanar import IsPlanar
except:
    print 'Could not import GraphPlanarity.so'
    print 'Compile by typing "make" in base of source directory.'
    IsPlanar=lambda x:True

DefaultFilters={}
MxAtm=0
maxWeight=0
maxRings=8
maxRingSize=8
RingSizeExceptions=1
UseRuleOf10=False

############ FILTERS

# BREDT VIOLATIONS
#AllFilters["Bredt's rule"]=NewFilter("Bredt's rule")
#BredtViolation=Chem.MolFromSmarts('[R&x2]@;=,:[R&x3](@[R&x2])@[R&x2]')
#BredtViolation=NewPatternFilter('bredt violation')
#BredtViolation.SetFilterPattern(Chem.MolFromSmarts('[R]@;=,:[R&x3](@[R])@[R]'))
#BredtViolation.SetExceptions([ Chem.MolFromSmarts('[R]@[R&x3](@[R])@[x3,x4]')])#,
#                           #    Chem.MolFromSmarts('[R]@;=,:[R&x3]@[R&x3]')])
#AllFilters["Bredt's rule"].SetFilterRoutine(BredtsRule)
#AllFilters["Bredt's rule"].SetFixRoutine(FixBredt)


#######################################
# Filters available to all flavors    #
DefaultFilters['Too big']=NewFilter('Too big')
def TooBig(mol):
    global maxWeight
    if MxAtm>0:
        if mol.GetNumHeavyAtoms()>MxAtm: return True
    if maxWeight>0:
        if Descriptors.MolWt(mol)>maxWeight:
            return True
    return False
DefaultFilters['Too big'].SetFilterRoutine(TooBig)

def CutMoreRings(mol):
    RI = mol.GetRingInfo()
    if not IsPlanar(mol): return True
    nrings,sa,sb=SSSR(mol,force=True)
    if UseRuleOf10 and RuleOf10(mol)>10: return True
    if RI.NumRings()>maxRings: return True
    return False

def CutRings(mol):
    import mutate
    #non-planar graphs and highly cyclic: keep breaking rings until graph
    #is planar or has an allowed number of rings
    changed=False
    itry=0
    while CutMoreRings(mol):
        itry+=1
        bondringids = mol.GetRingInfo().BondRings()
        # need to flatten bondringids:
        # fancy manner to flatten a list of tuples with possible duplicates
        bondringids = set.intersection(*map(set, bondringids))
        bonds=GetIBonds(bondringids, mol, notprop='group')
        #bonds=list(mol.GetBonds(OEAndBond(OEBondIsInRing(),
        #                                  mt.BndNotInGroup)))
        if len(bonds)==0: raise MutateFail()
        mutate.RemoveBond(mol, random.choice(bonds))
        if mol.HasProp('ringcount'): mol.ClearProp('ringcount')
        changed=True
        if itry >= mt.MAXTRY: raise MutateFail()
    return changed

DefaultFilters['Non-planar graph (Euler critereon)']=NewFilter('Non-planar graph (Euler critereon)')
def NonPlanarEuler(mol):
    # NOTE. Somehow TEMPO crashes on this 
    Radical=False
    nNodes=mol.GetNumHeavyAtoms()
    if nNodes >= 3 and mol.GetNumBonds()>3*nNodes-6 and not Radical:
        return 'Non-planar graph (Euler critereon)'
DefaultFilters['Non-planar graph (Euler critereon)'
          ].SetFilterRoutine(NonPlanarEuler)
DefaultFilters['Non-planar graph (Euler critereon)'].SetFixRoutine(CutRings)


DefaultFilters['Non-planar graph (Boyes)'
          ]=NewFilter('Non-planar graph (Boyes)')
def NotPlanarBoyes(mol):
    #print "In NotPlanarBoyes", 
    return not IsPlanar(mol)
DefaultFilters['Non-planar graph (Boyes)'
          ].SetFilterRoutine(NotPlanarBoyes)
DefaultFilters['Non-planar graph (Boyes)'].SetFixRoutine(CutRings)

DefaultFilters['Too many rings']=NewFilter('Too many rings')
def TooManyRings(mol):
    # Ring features
    if mol.GetNumAtoms()>maxRings:
        nrings,sa,sb=SSSR(mol)
        if sum(nrings)>maxRings:
            return 'Too Many Rings'
    return False
DefaultFilters['Too many rings'].SetFilterRoutine(TooManyRings)
DefaultFilters['Too many rings'].SetFixRoutine(CutRings)


fname='SSSR ring bigger than max allowed size'
DefaultFilters[fname]=NewFilter(fname)
def BiggestRing(mol):
    if mol.GetNumAtoms()>maxRingSize:
        nrings,sa,sb=SSSR(mol)
        if sum(nrings[maxRingSize-2:]) > RingSizeExceptions:
            failname='SSSR ring bigger than max allowed size'
            return failname
DefaultFilters[fname].SetFilterRoutine(BiggestRing)
def CutBiggestRings(mol):
    ringatoms, ringbonds=SSSR_GetRings(mol)
    bondids={bond.GetIdx():bond for bond in mol.GetBonds()}
    changed=False
    while sum(mol.GetData('ringcounts')[maxRingSize-2:])>RingSizeExceptions:
        toobig=set()
        smallenough=set()
        for ring in ringbonds:
            if len(ring)>maxRingSize: toobig.update(ring)
            else: smallenough.update(ring)
        canremove=toobig-smallenough
        removelist = [ bondids[id] for id in canremove
                            if not bondids[id].HasData('mygroup') ]
        if len(canremove)==0: raise MutateFail()
        else: mol.DeleteBond(  random.choice(removelist) )
        changed=True
        ringatoms,ringbonds=SSSR_GetRings(mol,True)
    return changed
DefaultFilters[fname].SetFixRoutine(CutBiggestRings)

LookUpFilter=NewFilter('Compound not in lookup table')
def lu(mol):
    smi=Chem.MolToSmiles(mol, True)
    return (smi not in lookUpTable)
LookUpFilter.SetFilterRoutine(lu)

LipinskiFilter=NewFilter('Lipinski violation')
def lp_routine(mol):
    return LipinskiRuleOf5(mol)>1
def LipinskiRuleOf5(mol):
    PropCalc(mol)
    ofs.flush()
    return mol.GetData('Lipinski violations')
LipinskiFilter.SetFilterRoutine(lp_routine)


RuleOf10Filter=NewFilter('Rule of 10')
def r10f_routine(mol):
    return RuleOf10(mol)>10
def RuleOf10(mol):
    OEPerceiveChiral(mol)
    ringcount,na,nb=SSSR(mol)
    nStereos=OECount(mol,OEIsChiralAtom())+OECount(mol,OEIsChiralBond())
    return sum(ringcount)+nStereos
RuleOf10Filter.SetFilterRoutine(r10f_routine)

SAScoreFilter=NewFilter("SA-Score synthetic accessibility")
def sascore_filt(mol):
    import os
    score=sa.SAScore(mol)
    if sa.SARescale:
        condition = score>SAScore
    else:
        condition = score<SAScore
    if condition:
        return 'SAScore: '+str(score)
    else:
        return False
SAScoreFilter.SetFilterRoutine(sascore_filt)

