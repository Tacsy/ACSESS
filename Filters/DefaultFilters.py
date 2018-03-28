#!/usr/bin/env python
#-*- coding: utf-8 -*-

from filters import NewFilter, NewPatternFilter
from rdkit import Chem
from rdkithelpers import *
from molfails import *
import random

AllFilters={}

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
MyFilters={}
MyFilters['Too big']=NewFilter('Too big')
tempmol=oe.OEMol()
def TooBig(mol):
    if MxAtm>0:
        if OECount(mol,OEIsHeavy())>MxAtm: return True
    if MaxWeight>0:
        if OECalculateMolecularWeight(mol)>MaxWeight: return True
    return False
MyFilters['Too big'].SetFilterRoutine(TooBig)

def CutMoreRings(mol):
    OEFindRingAtomsAndBonds(mol)
    if not IsPlanar(mol): return True
    nrings,sa,sb=SSSR(mol,force=True)
    if UseRuleOf10 and RuleOf10(mol)>10: return True
    if sum(nrings)>MaxRings: return True
    return False

def CutRings(mol):
    #non-planar graphs and highly cyclic: keep breaking rings until graph
    #is planar or has an allowed number of rings
    changed=False
    itry=0
    while CutMoreRings(mol):
        itry+=1
        bonds=list(mol.GetBonds(OEAndBond(OEBondIsInRing(),
                                          mt.BndNotInGroup)))
        if len(bonds)==0: raise MutateFail()
        mt.RemoveBond(mol, random.choice(bonds))
        if mol.HasData('ringcount'): mol.DeleteData('ringcount')
        changed=True
        if itry >= mt.MAXTRY: raise MutateFail()
    return changed

LookUpFilter=NewFilter('Compound not in lookup table')
def lu(mol):
    smi=oe.OECreateCanSmiString(mol)
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

