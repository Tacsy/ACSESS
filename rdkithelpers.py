#!/usr/bin/env python
#-*- coding: utf-8 -*-
debug=False

from rdkit import Chem
import traceback
from molfails import MutateFail

########### Some global variables 
bondorder={ 1:Chem.BondType.SINGLE,
            2:Chem.BondType.DOUBLE,
            3:Chem.BondType.TRIPLE }
aromatic=False

############ FROM CANONICAL #############################
def Finalize(mol, CanonicalTautomer=False, aromatic=aromatic):
    ''' This function makes sure that a new mutated/crossover/fixfilter molecule:
        - corrects the implicit valence
        - removes data from parent molecules.
        - optionally converts the molecule to its canonical tautomer
        - checks if it is a valid molecule, i.e. Sanitize step. (without aromaticity)
    '''
    mol.UpdatePropertyCache()
    if type(mol)==Chem.RWMol:
        RW=True
    else:RW=False
    ResetProps(mol)
    if CanonicalTautomer:
        Tautomerize(mol)
    try: Sanitize(mol, aromatic)
    except Exception as e:
        print "Error in Finalize with", Chem.MolToSmiles(mol, False), e,
        if debug:
            for item in traceback.extract_stack():
                print item
    #Chem.SetAromaticity(mol)
    return

def ResetProps(mol):
    ''' makes a fresh molecule without properties inherited from there parents
        by deleting all listed properties and resetting the SMILES string. '''
    isosmi = Chem.MolToSmiles(mol, True)
    for prop in [ 'filtered', 
            'hasstructure', 
            'tautomerized', 
            'minimized', 
            'selected', 
            'failed',
            'failedfilter',
            'Objective']:
        mol.ClearProp(prop)
    mol.SetProp('isosmi',isosmi)
    return

def Sane(mol, *args, **kwargs):
    ''' A SanitizeCheck function; function returns False if molecule is not sane '''
    try:
        Sanitize(mol, *args, **kwargs)
        return True
    except:
        return False

#def ToRWMol(mol):
#    if type(mol)==Chem.RWMol: return mol
    

############ NEW RDKIT HELPERS: #########################
flatten = lambda X:tuple(set(i for x in X for i in x))
def MolAtIsInGroup(mol, groupid):
    requirement = lambda atom:atom.HasProp('group') and atom.GetProp('group')==groupid
    return filter(requirement, mol.GetAtoms())
def MolBondIsInGroup(mol, groupid):
    requirement = lambda bond:bond.HasProp('group') and bond.GetProp('group')==groupid
    return filter(requirement, mol.GetBonds())
CanRemoveAtom=lambda atom:not atom.HasProp('protected') and not atom.HasProp('fixed')
CanChangeAtom=lambda atom:(not atom.HasProp('group')) or atom.HasProp('grouprep')

def GetSmallestRingSize(atom):
    return min([ i for i in range(1,20) if atom.IsInRingSize(i) ])

# Three Indices based getters:
def GetIAtoms(indices, mol, notprop=None):
    if notprop:
        requirement = lambda atom:atom.GetIdx() in indices and not atom.HasProp(notprop)
    else:
        requirement = lambda atom:atom.GetIdx() in indices
    return filter(requirement, mol.GetAtoms())
def GetIBonds(indices, mol, notprop=None):
    if notprop:
        requirement = lambda bond:bond.GetIdx() in indices and not bond.HasProp(notprop)
    else:
        requirement = lambda bond:bond.GetIdx() in indices
    return filter(requirement, mol.GetBonds())
def GetAtomIBonds(atomids, mol):
    ''' gets bonds based on the atom indices '''
    bonds=[]
    for bond in mol.GetBonds():
        if all( index in atomids for index in ( bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())):
            bonds.append(bond)
    return bonds

def GetAtoms(mol, notprop=None):
    if notprop:
        requirement = lambda atom:not atom.HasProp(notprop)
    else:
        requirement = lambda atom:True
    return filter(requirement, mol.GetAtoms())
def GetBonds(mol, notprop=None):
    if notprop:
        requirement = lambda bond:not bond.HasProp(notprop)
    else:
        requirement = lambda bond:True
    return filter(requirement, mol.GetBonds())

def GetXAtoms(mol, num=6, notprop=None):
    '''gets the number of atoms in mol with atomic num. default carbon, 6 and not has prop notprop
    for every atom it is checked if it has property indicated with the string notprop
    '''
    if notprop:
        requirement = lambda atom:atom.GetAtomicNum()==num and not atom.HasProp(notprop)
    else:
        requirement = lambda atom:atom.GetAtomicNum()==num
    return filter( requirement, mol.GetAtoms())

def GetFreeBonds(mol, order=None, notprop=None):
    ''' This function is used by aromatic ring addition. It returns the bonds with specified order
    and on both bond.atoms at least one hydrogen. '''
    # 1. select bonds with required bondorder
    if order:
        IsOrder = lambda bond:bond.GetBondType()==bondorder[order]
    else:
        IsOrder = lambda x:True
    ordbonds  = filter(IsOrder, mol.GetBonds())
    # 2. select bonds which have at least on H at each atom
    HasHs= lambda bond:all( atom.GetImplicitValence() for atom in (bond.GetBeginAtom(), bond.GetEndAtom()))
    withHbonds= filter(HasHs, ordbonds)
    # 3. test if not has prop notprop
    if notprop:
        bonds = filter(lambda bond:not bond.HasProp(notprop), withHbonds)
    else:
        bonds = withHbonds
    return bonds


########## Set List properties

def SetListProp(mol, name, iterable):
    ''' Since rdkit molecules can only store single values, list properties 
    are stored as strings. Values are separated by a space
    NB: float precision is hardcoded to 20 decimals'''
    string = ' '.join(['{:.20f}'.format(x) for x in iterable])
    mol.SetProp(name, string)
def GetListProp(mol, name):
    string = mol.GetProp(name)
    return map(float, string.split())

########## function potentially useful to avoid aromaticity

def Sanitize(mol, aromatic=False):
    '''The rdkit sanitize step with the option to switch off aromaticity'''
    if aromatic:
        Chem.SanitizeMol(mol)
    else:
        Chem.SanitizeMol(mol,
            sanitizeOps=Chem.SANITIZE_ALL^Chem.SANITIZE_SETAROMATICITY
            )
        #sanitizeOps=Chem.SANITIZE_ALL^Chem.SANITIZE_KEKULIZE^\
        #Chem.SANITIZE_SETAROMATICITY^Chem.SANITIZE_CLEANUP^\
        #Chem.SANITIZE_CLEANUPCHIRALITY)
    return

########## TAUTOMERIZING:

def Tautomerize(mol):
    try:
        if mol.GetBoolProp('tautomerized'): return
    except KeyError: pass
    smi1 = Chem.MolToSmiles(mol)
    from molvs import Standardizer
    s = Standardizer()
    try:
        s.standardize(mol)
    except ValueError as e:
        MutateFail(mol)
        return False
    #from molvs.tautomer import TautomerCanonicalizer
    #t = TautomerCanonicalizer()
    #t.canonicalize(mol)
    mol.SetBoolProp('tautomerized', True)
    smi2 = Chem.MolToSmiles(mol)

    if not smi1==smi2: print "tautomerized:", smi1, 'to:', smi2
    return True

####### Complex Ring Functions
ringsearch = {}
def SSSR(mol, force=False):

    if (not force) and mol.HasProp('ringcounts'):
        ringcounts = map(int, mol.GetProp('ringcounts').split())
        sharedatoms= mol.GetIntProp('sharedatoms')
        sharedbonds= mol.GetIntProp('sharedbonds')
        return (ringcounts,)+(sharedatoms, sharedbonds)
    #(map(int, mol.GetProp('ringcounts').split()),)+mol.GetIntProp('sharedcounts'

    ## SSSR ring analysis
    RI = mol.GetRingInfo()

    AssignedAtoms=set()
    AssignedBonds=set()
    SharedAtoms=set()
    SharedBonds=set()

    #nringatom= oe.OECount(mol, AtInRing )
    #nringbond= oe.OECount(mol, BondInRing )
    nringatom =len(flatten(RI.AtomRings()))
    nringbond =len(flatten(RI.BondRings()))

    nRings = [0] * max(nringatom ,8)

    #Loop over all possible ring sizes
    for i in xrange(3,nringatom+1):
        if (len(AssignedAtoms)==nringatom and
            len(AssignedBonds)==nringbond): break

        if not ringsearch.has_key(i):
            ringsearch[i] = Chem.MolFromSmarts('*~1'+'~*'*(i-1)+'1')

        #Find all instances of ring size i
        #matches=ringsearch[i].Match(mol,True)
        matches=mol.GetSubstructMatches(ringsearch[i])
        for match in matches:
            atomids=set(match)
            bondids=set([ bond.GetIdx() for bond in GetAtomIBonds(match, mol) ])
            #bondids=set( [bond.target.GetIdx()
            #              for bond in match.GetBonds() ] )

            #Count this ring only if some of its atoms or bonds
            #have not already been assigned to a smaller ring
            if not (atomids.issubset(AssignedAtoms) and
                    bondids.issubset(AssignedBonds)):
                nRings[i-3]+=1
                SharedAtoms.update(
                    atomids.intersection(AssignedAtoms))
                SharedBonds.update(
                    bondids.intersection(AssignedBonds))
                AssignedAtoms.update(atomids)
                AssignedBonds.update(bondids)

    mol.SetProp('ringcounts', " ".join(map(str, nRings)))
    mol.SetIntProp('sharedatoms',len(SharedAtoms))
    mol.SetIntProp('sharedbonds',len(SharedBonds))
    return nRings,len(SharedAtoms),len(SharedBonds)

#Smallest set of smallest rings analysis
#Just returns counts of various ring sizes. While SSSR in general
#is not invariant, these counts are.
ringsearch={}
def SSSR_GetRings(mol,force=False):

    if not force and mol.HasProp('numSSSRrings'):
        ringatoms=[]
        ringbonds=[]
        numring=mol.GetProp('numSSSRrings')
        for i in xrange(numring):
            ringatoms.append( set(map(int, mol.GetProp('AtomSSSR_'+str(i)).split())))
            ringbonds.append( set(map(int, mol.GetProp('BondSSSR_'+str(i)).split())))
        return ringatoms,ringbonds

    ## SSSR ring analysis
    RI = mol.GetRingInfo()
    #oe.OEFindRingAtomsAndBonds(mol)
    AssignedAtoms=set()
    AssignedBonds=set()
    SharedAtoms=set()
    SharedBonds=set()

    nringatom =len(flatten(RI.AtomRings()))
    nringbond =len(flatten(RI.BondRings()))
    nRings = [0] * max(nringatom,8)

    #Loop over all possible ring sizes
    ringatoms=[]
    ringbonds=[]
    for i in xrange(3,nringatom+1):
        if (len(AssignedAtoms)==nringatom and
            len(AssignedBonds)==nringbond): break
        if not ringsearch.has_key(i):
            ringsearch[i] = Chem.MolFromSmarts('*~1'+'~*'*(i-1)+'1')

        #Find all instances of ring size i
        matches=mol.GetSubstructMatches(ringsearch[i])
        for match in matches:
            if len(AssignedAtoms)==nringatom and \
               len(AssignedBonds)==nringbond: break
            atomids=set(match)
            bondids=set([ bond.GetIdx() for bond in GetAtomIBonds(mol, match) ])

            #Count this ring only if some of its atoms or bonds
            #have not already been assigned to a smaller ring
            if not (atomids.issubset(AssignedAtoms) and
                    bondids.issubset(AssignedBonds)):
                ringatoms.append( atomids )
                ringbonds.append( bondids )
                nRings[i-3]+=1
                SharedAtoms.update(
                    atomids.intersection(AssignedAtoms))
                SharedBonds.update(
                    bondids.intersection(AssignedBonds))
                AssignedAtoms.update(atomids)
                AssignedBonds.update(bondids)

    mol.SetProp('ringcounts', " ".join(map(str, nRings)))
    mol.SetIntProp('sharedatoms',len(SharedAtoms))
    mol.SetIntProp('sharedbonds',len(SharedBonds))
    mol.SetIntProp('numSSSRrings', len(ringatoms) )
    for i in xrange(len(ringatoms)):
        mol.SetProp('AtomSSSR_'+str(i), " ".join(map(str,list(ringatoms[i]))))
        mol.SetProp('BondSSSR_'+str(i), " ".join(map(str,list(ringbonds[i]))))
    return ringatoms, ringbonds


