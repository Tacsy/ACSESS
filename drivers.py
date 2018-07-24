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
import random, sys
from output import logtime, StartTimer, EndTimer, stats
import output

###### mutation probabilities:
p_BondFlip = 0.8
p_AtomFlip = 0.8
p_RingAdd = 0.1
p_AddFreq = 0.5  #Actual probability: (.8*.7)*.5=.28
p_DelFreq = 0.8  #Actual probability: .224
p_RingRemove = 0.2  #actual probability=.8*.3=.24
p_AddAroRing = 0.5  #actual probability is rather low sing free double bonds aren't prevalent
p_AddFusionRing = 0.5  #actual probability is rather low sing free double bonds aren't prevalent
MutateStereo = False
StereoFlip = 0.2

##### workflow switches:
startFilter = 2
startTautomer = None
startGenStruc = 0
KeepNoGeomPool = True

##### other variables
maxPool = 1000

debug = False


def SetIterationWorkflow(gen):
    Tautomerizing = (startTautomer and gen >= startTautomer)
    Filter = (gen >= startFilter)
    GenStruc = (gen >= startGenStruc)
    return Tautomerizing, Filter, GenStruc


@logtime()
def DriveMutations(lib):
    global stats
    nDups, nExcp, nCand = (0, 0, 0)
    # 1. CROSSOVERS:
    print "crossovers...",
    sys.stdout.flush()
    newmols = []
    for i in xrange(min(mprms.nCross, len(lib) - 1)):
        try:
            if i < mprms.nCross * mprms.EdgeRatio and mprms.EdgeLen > 0:
                mol1 = random.choice(lib[:mprms.EdgeLen])
                mol2 = random.choice(lib + newmols)
            else:
                mol1 = random.choice(lib + newmols)
                mol2 = random.choice(lib + newmols)
            candidate = mutate.Crossover(mol1, mol2)
            if type(candidate) == Chem.RWMol: candidate = candidate.GetMol()
            Finalize(candidate, aromatic=False)
            newmols.append(candidate)
            stats['nCand'] += 1
        except MutateFail:
            stats['nExcp'] += 1
    newmols = filter(Sane, newmols)

    nbefore = len(newmols)
    newmols = RemoveDuplicates(newmols)
    stats['nDups'] += nbefore - len(newmols)

    # 2. MUTATIONS
    print "mutating...",
    sys.stdout.flush()
    for i in xrange(mprms.nMut):
        if i < mprms.nMut * mprms.EdgeRatio and mprms.EdgeLen > 0:
            candidate = Chem.Mol(random.choice(lib[:mprms.EdgeLen]))
        else:
            candidate = Chem.Mol(random.choice(lib + newmols))
        try:
            candidate = MakeMutations(candidate)
            if type(candidate) == Chem.RWMol:
                candidate = candidate.GetMol()
            try:
                Finalize(candidate, aromatic=False)
            except Exception as e:
                raise MutateFail(candidate)
            newmols.append(candidate)
            stats['nCand'] += 1
        except MutateFail:
            stats['nExcp'] += 1
    newmols = filter(Sane, newmols)
    nbefore = len(newmols)
    newmols = RemoveDuplicates(newmols)
    stats['nDups'] = nbefore - len(newmols)
    if debug:
        with open('mutatelib', 'w') as f:
            for mol in newmols:
                f.write(Chem.MolToSmiles(mol) + '\n')

    return newmols


###############################################
# Drive filtering
def DriveFilters(lib, Filtering, GenStrucs):
    global stats

    def torwmol(mol):
        if not type(mol) == Chem.RWMol: return Chem.RWMol(mol)
        else: return mol

    lib = map(torwmol, lib)

    #nSizeCut,nDups,nFilt=(0,0,0)
    nbefore = len(lib)
    #lib=[mol for mol in lib if mol.GetBoolProp('filtered') or not fl.TooBig(mol) ]
    #nSizeCut += nbefore-len(lib)
    if not (Filtering or GenStrucs): return lib
    ####### Parallel filtering ######
    if Filtering and mprms.mpi:
        import parallel as pl
        if GenStrucs:
            pl.ScatterFixFilterStruc(lib)
        else:
            pl.ScatterFixFilter(lib)
    #################### Serial filtering #####################
    elif Filtering:
        print "filtering...",
        sys.stdout.flush()
        StartTimer('Filters')
        for mol in lib:
            if not mol.HasProp('filtered'):
                filters.FixAndFilter(mol)
                if debug:
                    if not mol.HasProp('failedfilter'): print "Jos Error",
                    else: print "ff:", mol.GetProp('failedfilter'),
            else:
                assert mol.GetBoolProp('filtered') == True
                mol.SetProp('failedfilter', '')
        nbefore = len(lib)
        lib = RemoveDuplicates(lib)
        stats['nDups'] += nbefore - len(lib)
        EndTimer('Filters')

        #### Generate 3d structures if requested ###
        if GenStrucs:
            StartTimer('GenStrucs')
            print "genstrucs..."
            sys.stdout.flush()
            for i, mol in enumerate(lib):
                if sys.stdout.isatty():
                    sys.stdout.write("{}%\r".format(int(100 * i / len(lib))))
                    sys.stdout.flush()
                if (mol.HasProp('hasstructure')
                        or mol.GetProp('failedfilter')):
                    continue

                from helpers import xyzfromrdmol
                try:
                    xyz = xyzfromrdmol(mol)
                    mol.SetProp('hasstructure', xyz)
                    if filters.GeomFilter(mol) == 'SAV':
                        mol.SetProp('failedfilter', 'SAV')
                except NoGeom:  #if structure can't be generated, filter it
                    mol.SetProp('failedfilter', 'failed to generate geometry')
                except MutateFatal:
                    from pprint import pprint as pp
                    print "Mutate Fatal Error with:", pp(mol.__dict__)
                    print "smiles:", Chem.MolToSmiles(mol), "end smi"
                    mol.SetProp('failedfilter', 'failed to generate geometry')
                    raise
            EndTimer('GenStrucs')

    ## Remove filtered compounds from library ##
    if Filtering:
        for mol in lib:
            try:
                failed = mol.GetProp('failedfilter')
            except KeyError:
                print "no failed filter:", Chem.MolToSmiles(mol)
            if failed:
                stats['nFilt'] += 1
                try:
                    smi = Chem.MolToSmiles(mol)
                except NotImplementedError:
                    smi = mol._smiles
                output.filterFile.write(smi + '  ' + failed + '\n')
        lib = [mol for mol in lib if not mol.GetProp('failedfilter')]
        output.filterFile.flush()
    return lib


def DrivePoolFilters(pool, Filtering, GenStruc, Tautomerizing, gen):
    pool = filter(bool, pool)
    if Tautomerizing:
        CanonicalTautomer = True
        if gen == startTautomer:
            print 'tautomerizing pool...',
            sys.stdout.flush()
            pool = filter(Tautomerize, pool)
    if GenStruc and gen == startGenStruc and not KeepNoGeomPool:
        print 'restarting and filtering pool...',
        sys.stdout.flush()
        pool = [mol for mol in pool if mol.GetBoolProp('hasstructure')]
        DriveFilters(pool, Filtering, GenStruc)
    if (Filtering and gen == startFilter) or (GenStruc and gen == startGenStruc
                                              and KeepNoGeomPool):
        print 'pool-',
        pool = DriveFilters(pool, Filtering, GenStruc)
    return pool


@logtime()
def ExtendPool(pool, lib, newmol):
    def key(x):
        if not x.HasProp('selected'): x.SetIntProp('selected', 0)
        return x.GetIntProp('selected')

    from random import shuffle
    if len(pool) + len(newmol) <= maxPool:
        pool += newmol
        return RemoveDuplicates(pool)
    else:
        newpool = RemoveDuplicates(
            newmol + lib)  #this is the minumum pool content
        if len(newpool) < maxPool:
            shuffle(pool)
            oldpool = sorted(set(pool) - set(lib), key=key, reverse=True)
            newpool = newpool + oldpool[:maxPool - len(newpool)]
        return RemoveDuplicates(newpool)


@logtime()
def DriveSelection(pool, subsetSize):
    print "selecting...",
    sys.stdout.flush()
    #1. select maximin algorithm.
    if mprms._similarity:
        from similarity import FPMaximin
        lib = FPMaximin(pool, mprms.subsetSize)
    else:
        from distance import Maximin
        lib = Maximin(pool, mprms.subsetSize)
    return lib


############ EXTRA FUNCTIONS: #############


##############################################################
# Simple utility to remove duplicate structures
# Tests for identical molecules by generating smiles strings
# Take the molecule with more information generated
def GetScore(mol):
    score = 0
    for prop in ['filtered', 'hasstructure', 'Objective', 'tautomerized']:
        try:
            score += bool(mol.GetProp(prop))
        except KeyError:
            pass
    if not mol.HasProp('selected'): mol.SetIntProp('selected', 0)
    return score


def RemoveDuplicates(lib):
    lib.sort(key=lambda x: x.GetProp('isosmi'))
    i = 1
    while i < len(lib):
        #print len(lib),
        if lib[i].GetProp('isosmi') == lib[i - 1].GetProp('isosmi'):
            iscore = GetScore(lib[i])
            i1score = GetScore(lib[i - 1])

            if iscore >= i1score:
                lib[i].SetIntProp('selected', lib[i].GetIntProp('selected') +
                                  lib[i - 1].GetIntProp('selected'))
                lib.pop(i - 1)
            else:
                lib[i - 1].SetIntProp('selected',
                                      lib[i].GetIntProp('selected') +
                                      lib[i - 1].GetIntProp('selected'))
                lib.pop(i)

        else:
            i += 1
    return lib


############ MUTATIONS INTERFACE: ############


# Mutation driver
#@captureMolExceptions
def MakeMutations(candidate):
    ''' The mutations, based on not-aromatic SMILES'''
    # 1. Kekulize:
    try:
        Chem.Kekulize(candidate, True)
    except ValueError:
        print "MakeMutation. Kekulize Error:", Chem.MolToSmiles(candidate)
        raise MutateFail(candidate)
    # 2. Mutate
    candidate = SingleMutate(candidate)
    # 3. SetAromaticity again
    try:
        if aromatic: Chem.SetAromaticity(candidate)
        Finalize(candidate)
    except ValueError:
        print "MakeMutation. SetAromaticity Error:", Chem.MolToSmiles(
            candidate)
        raise MutateFail(candidate)
    return candidate


def SingleMutate(candidateraw):
    #############################################################
    # I still need to test if the candidate is REALLY mutated
    # when actually a RWMol is mutated!.
    # maybe we have to return the candidate or something?
    #############################################################

    if candidateraw is None: raise ValueError('candidate is none')
    else: candidate = Chem.RWMol(candidateraw)
    global stats

    parent = candidate.GetProp('isosmi')
    ResetProps(candidate)
    change = False
    candidate.SetProp('parent', parent)

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
    bonds = list(GetBonds(candidate, notprop='group'))
    if random.random() < p_BondFlip and len(bonds) > 0:
        if debug: print "1",
        stats['nFlip'] += 1
        change = True
        try:
            Chem.Kekulize(candidate, True)
            bonds = list(GetBonds(candidate, notprop='group'))
            mutate.FlipBond(candidate, random.choice(bonds))
            Finalize(candidate, aromatic=False)
        except MutateFail:
            stats['nFlipFail'] += 1

    # 2. Flip atom identity
    atoms = filter(CanChangeAtom, candidate.GetAtoms())
    if random.random() < p_AtomFlip and len(atoms) > 0:
        if debug: print "2",
        stats['nAtomType'] += 1
        change = True
        try:
            mutate.SwitchAtom(candidate, random.choice(atoms))
            try:
                Finalize(candidate, aromatic=False)
            except:
                raise MutateFail
        except MutateFail:
            stats['nAtomTypeFail'] += 1

    # 3. add a ring - either a new aromatic ring, or bond
    # Aromatic ring addition disabled - probably not necessary
    if random.random() < p_RingAdd:
        if debug: print "3",
        stats['nNewRing'] += 1
        try:
            mutate.AddBond(candidate)
            try:
                Finalize(candidate, aromatic=False)
            except:
                raise MutateFail
            return candidate
        except MutateFail:
            stats['nNewRingFail'] += 1

    # 4. Try to remove a bond to break a cycle
    if random.random() < p_RingRemove:
        if debug: print "4",
        bondringids = candidate.GetRingInfo().BondRings()
        # need to flatten bondringids:
        if bondringids:
            # fancy manner to flatten a list of tuples with possible duplicates
            bondringids = set.intersection(*map(set, bondringids))
            bonds = GetIBonds(bondringids, candidate, notprop='group')
        else:
            bonds = []
        stats['nBreak'] += 1
        if len(bonds) != 0:
            mutate.RemoveBond(candidate, random.choice(bonds))
            try:
                Finalize(candidate, aromatic=False)
            except Exception as e:
                print "error with RemoveBond with mol:", Chem.MolToSmiles(True)
                print e
                raise MutateFail
            return candidate
        else:
            stats['nBreakFail'] += 1

    # 5. add an atom
    atoms = GetAtoms(candidate, notprop='protected')
    bonds = GetBonds(candidate, notprop='group')
    if (random.random() < p_AddFreq and len(atoms) + len(bonds) > 0):
        if debug: print "5",
        stats['nAdd'] += 1
        try:
            mutate.AddAtom(candidate, random.choice(atoms + bonds))
            try:
                Finalize(candidate, aromatic=False)
            except:
                raise MutateFail
            return candidate
        except MutateFail:
            stats['nAddFail'] += 1

    # 6. remove an atom
    if len(atoms) > 1 and random.random() < p_DelFreq:
        if debug: print "6",
        inismi = Chem.MolToSmiles(candidate)
        stats['nRemove'] += 1
        Chem.Kekulize(candidate, True)
        atoms = filter(CanRemoveAtom, candidate.GetAtoms())
        try:
            try:
                mutate.RemoveAtom(candidate, random.choice(atoms))
                Finalize(candidate, aromatic=False)
            except Exception as e:
                print "initial mol:", inismi
                print "e:", e, "mol:", Chem.MolToSmiles(candidate)
                raise MutateFail
            return candidate
        except MutateFail:
            stats['nRemoveFail'] += 1

    # 7. Add aromatic ring to a bond
    if random.random() < p_AddAroRing:
        freesinglebonds = GetFreeBonds(candidate, order=1, sides=True)
        #print "freesinglebonds:", freesinglebonds
        freedoublebonds = GetFreeBonds(candidate, order=2, notprop='group')
        triplebonds = filter(lambda bond: bond.GetBondType() == bondorder[3],
                             candidate.GetBonds())
        correctbonds = freedoublebonds + triplebonds + freesinglebonds
        if len(correctbonds) > 1:
            #print "n freedoublebonds:", len(correctbonds)
            if debug: print "7",
            stats['nAddArRing'] += 1
            Chem.Kekulize(candidate, True)
            try:
                candidate = mutate.AddArRing(candidate,
                                             random.choice(correctbonds))
                try:
                    Finalize(candidate, aromatic=False)
                except Exception as e:
                    print "exception in AddArRing with", Chem.MolToSmiles(
                        candidate)
                    print e
                    raise MutateFail
                return candidate
            except MutateFail:
                stats['nAddArRingFail'] += 1

    # 8. Add aromatic ring to two rings
    if random.random() < p_AddFusionRing:
        try:
            p = Chem.MolFromSmarts('[h]@&=*(@*)@[h]')
            matches = candidate.GetSubstructMatches(p)
        except RuntimeError:
            stats['nAddFusionRingFail'] += 1
        else:
            if matches:
                stats['nAddFusionRing'] += 1
                match = random.choice(matches)
                try:
                    candidate = mutate.AddFusionRing(candidate, match)
                    return candidate
                except MutateFail:
                    stats['nAddFusionRingFail'] += 1

    if not change:
        stats['nNoMutation'] += 1
        raise MutateFail()
    #if type(candidate)==Chem.RWMol: candidate=candidate.GetMol()
    return candidate
