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

    Every filterflavor-module has a dictionary which is formatted as:
    filterfunction might be a filterclass with a __call__ attribute though.
    {'filtername1':filterfunction1, 'filtername2':filterfunction2, ... }
'''

ActiveFilters = dict()
FilterFlavor = 'GDB'
extrafilters = []
MAXTRY=10
SAScore=0.0

############################################################
#       Functions from Filter.py
############################################################

def Init():
    global ActiveFilters, FilterFlavor

    # 1. load default filters into ActiveFilters:
    from Filters import DefaultFilters
    ActiveFilters.update(DefaultFilters.DefaultFilters)
    if SAScore:
        ActiveFilters['SAScore']=SAScoreFilter
        SAScore.Init()

    # 2. load the main set of filters of either GDB/Druglike/etc.
    if FilterFlavor=='GDB':
        from Filters import GDBFilters as FilterFlavor
    else:
        raise TypeError('FilterFlavor not recognized')
    ActiveFilters.update(FilterFlavor.AllFilters)

    # 3. load extra specific filters
    if extrafilters:
        from Filters.ExtraFilters import ExtraFilters
        for extrafilter in extrafilters:
            if extrafilter in ['Qiu', 'qiu']:
                ActiveFilters['qiu1']=ExtraFilters['qiu1']
                ActiveFilters['qiu2']=ExtraFilters['qiu2']
            else:
                ActiveFilters[extrafilter]=ExtraFilters[extrafilter]
        print "{} added".format(extrafilter)
    for AcFil in sorted(ActiveFilters.keys()): print AcFil

    return

def FixAndFilter(mol):
    changed=False
    #failure=False
    for i in xrange(MAXTRY):
        # in old version: here prepare for 2D filters
        # run through all filters:
        for filtername in ActiveFilters:
            failure = ActiveFilters[filtername](mol)
            if filtername=='allene':
                print 'allene', type(ActiveFilters[filtername]), 
                print Chem.MolToSmiles(mol)
                print failure
            mol.SetBoolProp('failed', bool(failure))
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
        if i==MAXTRY-1 or (not failure): #i.e. all filters passed without problems
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
            match=mol.HasSubstructMatch(self.pattern)
            if match:
                #if self.name=='allene':print "match:", match
                return self.name
            else:
                return False

    def __repr__(self):
        #return "{}\n".format(self.name, self.pattern.GetPattern().GetTitle())
        return "{} {}".format(self.name, Chem.MolToSmarts(self.pattern))

    def Fix(self,mymol):
        if self.HasFix:
            return self.fixroutine(mymol,self)
        else: return False

    #def SetFilterRoutine(self,*args,**kwargs): raise NotImplementedError()

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

# THIS FUNCTION IS NOT USED YET!
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

