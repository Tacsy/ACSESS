#!/usr/bin/env python
#-*- coding: utf-8 -*-
from rdkit import Chem
from rdkit.Chem import AllChem

from molfails import MutateFail
from rdkithelpers import *

debug=False
'''
This filter.py contains all the possible filters in the previous version of
ACSESS, in order to make the code in the cleaner way.

previous filters included:
    Filters.py
    GDBFilters.py
    DruglikeFilters.py
    OEFilters.py
    SMUFilters.py
    Filters_BC.py
    QiuFilter.py

Note from Jos:
    I don't think it is possible to include and combine every filter into this module
    since this module would become unmanageable. We could make a filters subfolder to 
    order them. 
'''
'''
Some clearification:
    all filter-related function are acting on single molecule and return a
    bool value to indicate whether a molecule should be filtered out (True)
    or not (False)
'''

ActiveFilters = dict()
FilterFlavor = 'GDB'
extrafilters = []
MAXTRY=10

############################################################
#       Functions from Filter.py
############################################################

def Init():
    global ActiveFilters, FilterFlavor

    if FilterFlavor=='GDB':
        from Filters import GDBFilters as FilterFlavor
    else:
        raise TypeError('FilterFlavor not recognized')
    # so every filterflavor-module has a dictionary which is formatted as:
    # filterfunction might be a filterclass with a __call__ attribute though.
    # {'filtername1':filterfunction1, 'filtername2':filterfunction2, ... }
    ActiveFilters.update(FilterFlavor.AllFilters)
    if extrafilters:
        from Filters.ExtraFilters import ExtraFilters
        for extrafilter in extrafilters:
            ActiveFilters[extrafilter]=ExtraFilters[extrafilter]
        print "{} added".format(extrafilter)
    print "ActiveFilters:", ActiveFilters

    return

def FixAndFilter(mol):
    changed=False
    failure=False
    for _ in xrange(MAXTRY):
        # in old version: here prepare for 2D filters
        # run through all filters:
        for filtername in ActiveFilters:
            failure = ActiveFilters[filtername](mol)
            if failure:
                if debug: print "in F&F:", failure
                try:
                    success= ActiveFilters[filtername].Fix(mol)
                except MutateFail:
                    success=False
                    changed=True
                if not success:
                    return changed, failure
                else:
                    changed=True
                    Finalize(mol)
                    break
        if not failure: break #i.e. all filters passed without problems
    return changed, failure

##################################
# These use to be in Classes.py: #
##################################

class NewFilter(object):
    def __init__(self, name):
        self.name=name
        self.HasFix=False

    def __call__(self, mol):
        filtered=self.function(mol)
        if filtered:
            if type(filtered)==bool or filtered==1: return self.name
            else: return filtered
        else: return 0

    def __repr__(self):
        return "NewFilter: {}\n".format(self.name)

    def Fix(self,mol):
        if self.HasFix:
            return self.fixroutine(mol)
        else: return False

    def SetFilterRoutine(self,filterroutine):
        self.function=filterroutine

    def SetFixRoutine(self,fixroutine):
        self.fixroutine=fixroutine
        self.HasFix=True

class NewPatternFilter(NewFilter):
    def __init__(self, name):
        self.name=name
        self.HasExceptions=False
        self.HasFix=False

    def __call__(self, mol):
        if self.HasExceptions:
            if len(self.FilterWithExceptions(mol))>0:
                return self.name
            else: return False
        else:
            #if self.pattern.SingleMatch(mymol): return self.name
            if mol.HasSubstructMatch(self.pattern): return self.name
            else: return False

    def __repr__(self):
        #return "{}\n".format(self.name, self.pattern.GetPattern().GetTitle())
        return "{} {}".format(self.name, Chem.MolToSmarts(self.pattern))

    def Fix(self,mymol):
        if self.HasFix:
            return self.fixroutine(mymol,self)
        else: return False

    def SetFilterRoutine(self,*args,**kwargs): raise NotImplementedError()

    def SetFilterPattern(self,pattern):
        self.pattern=pattern

    def SetExceptions(self,exc):
        self.HasExceptions=True
        try: len(exc)
        except TypeError: exc=[exc]
        self.MyExceptions=exc

    def FilterWithExceptions(self, mol):
        #if oe.OEHasExplicitHydrogens(mymol) and not Radical:
        #    print 'WARNING: Explicit hydrogens detected in FilterWithExceptions'
        #    print '     Filters may not function properly'
        #    print mymol.GetData('isosmi')

        #matches=list( self.pattern.Match(mymol,True) )
        matches = list( mol.GetSubstructMatches(self.pattern))
        if not self.HasExceptions: return matches

        # RDKit automattically return a tuple of matches, each already a tuple.
        matches=[ set(match) for match in matches ]

        # Remove matches that are substructures of exceptions
        # NOT TESTED:
        for exception in self.MyExceptions:
            exmatches= mol.GetSubstructMatches(exception)
            for exmatch in exmatches: # for each exception substructure found:
                matches= [match for match in matches if not match.issubset( exmatches ) ]
                if len(matches)==0: break
            if len(matches)==0: break
        return matches

# ensure molecule has specific pattern
def CheckSubstructure(mol, patterns):
    # get mol and return bool for filter or not
    # patterns: list of pattern for substructure search 
    dumpMol = False
    # for substructure search
    for pattern in xrange(len(patterns)):
        if mol.HasSubstructMatch(pattern):
            dumpMol = True
            break

    return dumpMol

############################################################
#       Functions from QiuFilter.py
############################################################

# ensure molecule do not have ring with ring size N
def CheckRingSize(mol, ringSize):
    # get mol and return bool for filter or not
        
    dumpMol = False
    # calculate smallest set of rings (SSR)
    ssr = Chem.GetSymmSSSR(mol)
    for i in range(len(ssr)):
        if len(list(ssr[i])) == ringSize:
            dumpMol = True

    return dumpMol

# ensure molecule do not have a bond with bondorder N
def CheckBondOrder(mol, bondOrder):
    # get mol and return bool for filter or not  
    # bondOrder should be double, 1.0 for single, 1.5 for aromatic, 2.0 for
    # double, 3.0 for triple

    dumpMol = False
    # loop over bonds in molecule and check bond order
    for bond in mol.GetBonds():
        if bond.GetBondTypeAsDouble() == bondOrder:
            dumpMol = True

    return dumpMol
