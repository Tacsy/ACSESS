#!/usr/bin/env python
#-*- coding: utf-8 -*-
debug = False
from rdkit import Chem
from rdkit.Chem import AllChem

from molfails import MutateFail
from rdkithelpers import *

import mprms
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
GeomFilter = lambda x: False
MAXTRY = 10
SAScore = 0.0
verbose = True

############################################################
#       Functions from Filter.py
############################################################


def Init():
    global ActiveFilters, FilterFlavor, GeomFilter

    # 1. load default filters into ActiveFilters:
    from Filters import DefaultFilters
    ActiveFilters.update(DefaultFilters.DefaultFilters)
    if SAScore:
        from Filters import SAS
        ActiveFilters['SAScore'] = DefaultFilters.SAScoreFilter
        #SAS.Init()

    # 2. load the main set of filters of either GDB/Druglike/etc.
    if FilterFlavor == 'GDB':
        from Filters import GDBFilters as FilterFlavor
        from Filters.GDBFilters import GeomFilter
    else:
        raise TypeError('FilterFlavor not recognized')
    ActiveFilters.update(FilterFlavor.AllFilters)

    # 3. load extra specific filters
    if extrafilters:
        from Filters.ExtraFilters import ExtraFilters
        for extrafilter in extrafilters:
            if extrafilter in ['Qiu', 'qiu']:
                ActiveFilters['qiu1'] = ExtraFilters['qiu1']
                ActiveFilters['qiu2'] = ExtraFilters['qiu2']
                print "qiu filters added".format(extrafilter)
            else:
                if extrafilter in ExtraFilters:
                    ActiveFilters[extrafilter] = ExtraFilters[extrafilter]
                    print "{} added from ExtraFilters".format(extrafilter)
                else:
                    ActiveFilters[extrafilter] = NewFilter(extrafilter)
                    ActiveFilters[extrafilter].SetFilterRoutine(getattr(mprms, extrafilter))
                    print "{} added from user input".format(extrafilter)
    if debug:
        for AcFil in sorted(ActiveFilters.keys()):
            print AcFil

    return


def FixAndFilter(mol):
    if mol is None: return None, True
    elif not type(mol) == Chem.RWMol:
        print "mol to RWmol"
        mol = Chem.RWMol(mol)

    # 1. set on aromaticity
    try:
        Chem.SetAromaticity(mol)
    except Exception as e:
        if verbose:
            print "didn't manage to set aromaticity for:",
            print Chem.MolToSmiles(mol)

    # 2. filter
    changed, filt = FixFilters(mol)

    # 3. Switch off aromaticity:
    try:
        Chem.Kekulize(mol, True)
    except Exception as e:
        if verbose:
            print "didn't manage to kekulize:",
            print Chem.MolToSmiles(mol),
        filt = 'unkekulizable'

    # 4. Sanitize if Fixroutines changed the molecule
    if debug and filt:
        print "changed:{}, failed:{}, mol:{}".format(changed, filt,
                                                     Chem.MolToSmiles(mol))
    if changed: Finalize(mol)

    # 5. Set molprops
    mol.SetBoolProp('filtered', True)
    if type(filt) is bool: filt = {True: 'unknown', False: ''}[filt]
    mol.SetProp('failedfilter', filt)
    return changed, filt


def FixFilters(mol):
    changed = False
    failure = False
    #if SetSulfurState(mol): changed=True
    #First, try geometry filter
    if mol.HasProp('hasstructure'):
        failure = GeomFilter(mol)
        if failure: return changed, failure

    lastfailure = ''
    for i in xrange(MAXTRY):
        #Run through all filters
        for ft in ActiveFilters:
            failure = ActiveFilters[ft](mol)
            if failure:  #try to fix the problem
                if debug:
                    print '{} failure {} with {}'.format(
                        ft, failure, Chem.MolToSmiles(mol)),
                # 1. Fixes are based on kekulized forms of the molecules:
                try:
                    Chem.Kekulize(mol, True)
                except ValueError:
                    #raise MutateFail
                    success = False
                    changed = True
                if debug: print Chem.MolToSmiles(mol)
                # 2. Fix:
                try:
                    success = ActiveFilters[ft].Fix(mol)
                except (MutateFail, ValueError) as e:
                    success = False
                    changed = True
                # 3. Force set back to Aromatic. It that fails->fail
                try:
                    Chem.SetAromaticity(mol)
                except ValueError:
                    success = False
                    changed = True
                # 4. if error->fail. otherwise new change->refilter(break)
                if not success: return changed, failure
                else:
                    changed = True
                    # 4b. To get a new change it should stand the test:
                    try:
                        Finalize(mol)
                    except ValueError:
                        return changed, failure

                    if debug: print "\tFixed!", Chem.MolToSmiles(mol)
                    break
        if i == MAXTRY - 1 or (not failure):
            return changed, failure


def NewFixAndFilter(mol):
    changed = False
    if mol is None: return None, True
    elif not type(mol) == Chem.RWMol: mol = Chem.RWMol(mol)
    for i in xrange(MAXTRY):
        for filtername in ActiveFilters:
            try:
                failure = ActiveFilters[filtername](mol)
            except AttributeError:
                print "filtername:", filtername
                print "mol:", Chem.MolToSmiles(mol)
                raise SystemExit('error')
            mol.SetBoolProp('failed', bool(failure))
            if failure:
                if debug: print "in F&F:", failure
                try:
                    success = ActiveFilters[filtername].Fix(mol)
                except MutateFail:
                    success = False
                    changed = True
                if not success:
                    return changed, failure
                else:
                    changed = True
                    Finalize(mol)
                    break
        if i == MAXTRY - 1 or (
                not failure):  #i.e. all filters passed without problems
            return changed, failure


##################################
# These use to be in Classes.py: #
##################################


class NewFilter(object):
    def __init__(self, name):
        self.name = name
        self.HasFix = False

    def __call__(self, mol):
        filtered = self.function(mol)
        if filtered:
            if type(filtered) == bool or filtered == 1:
                return self.name
            else: return filtered
        else: return 0

    def __repr__(self):
        return "NewFilter: {}\n".format(self.name)

    def Fix(self, mol):
        if self.HasFix:
            return self.fixroutine(mol)
        else:
            return False

    def SetFilterRoutine(self, filterroutine):
        self.function = filterroutine

    def SetFixRoutine(self, fixroutine):
        self.fixroutine = fixroutine
        self.HasFix = True


class NewPatternFilter(NewFilter):
    def __init__(self, name):
        self.name = name
        self.HasExceptions = False
        self.HasFix = False

    def __call__(self, mol):
        if self.HasExceptions:
            if len(self.FilterWithExceptions(mol)) > 0:
                return self.name
            else:
                return False
        else:
            match = mol.HasSubstructMatch(self.pattern)
            if match:
                #if self.name=='allene':print "match:", match
                return self.name
            else:
                return False

    def __repr__(self):
        #return "{}\n".format(self.name, self.pattern.GetPattern().GetTitle())
        return "{} {}".format(self.name, Chem.MolToSmarts(self.pattern))

    def Fix(self, mymol):
        if self.HasFix:
            if debug:
                print "in PatternFix:", self.name, self.fixroutine.__name__
            return self.fixroutine(mymol, self)
        else:
            return False

    #def SetFilterRoutine(self,*args,**kwargs): raise NotImplementedError()

    def SetFilterPattern(self, pattern):
        self.pattern = pattern

    def SetExceptions(self, exc):
        self.HasExceptions = True
        try:
            len(exc)
        except TypeError:
            exc = [exc]
        self.MyExceptions = exc

    def FilterWithExceptions(self, mol):
        matches = list(mol.GetSubstructMatches(self.pattern))
        if not self.HasExceptions: return matches

        # RDKit automattically return a tuple of matches, each already a tuple.
        matches = [set(match) for match in matches]

        # Remove matches that are substructures of exceptions
        # NOT TESTED:
        for exception in self.MyExceptions:
            exmatches = mol.GetSubstructMatches(exception)
            for exmatch in exmatches:  # for each exception substructure found:
                matches = [
                    match for match in matches if not match.issubset(exmatches)
                ]
                if len(matches) == 0: break
            if len(matches) == 0: break
        return matches


# FUNCTIONS BELOW ARE NOT USED YET!
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


def CheckLibrary(mol, library):
    # get mol and return bool for filter or not
    # library: list of molecule string SMILES to check whether or not the
    # molecue is in the library or not
    dumpMol = False
    # library search:
    if Chem.MolToSmiles(mol) in library:
        dumpMol = True

    return dumpMol
