#!/usr/bin/env python
#-*- coding: utf-8 -*-

from rdkit import Chem
import traceback
debug=True

bondorder={ 1:Chem.BondType.SINGLE,
            2:Chem.BondType.DOUBLE,
            3:Chem.BondType.TRIPLE }

############ FROM CANONICAL #############################
def Finalize(mol):
    if type(mol)==Chem.RWMol:
        RW=True
    else:RW=False
    try: Chem.SanitizeMol(mol)
    except Exception as e:
        if debug:
            for item in traceback.extract_stack(): print item
        print "Error in Finalize with", Chem.MolToSmiles(mol), e
    #Chem.SetAromaticity(mol)
    ResetProps(mol)
    return

def ResetProps(mol):
    # should 'isosmi' be included?
    isosmi = Chem.MolToSmiles(mol, True)
    for prop in [ 'filtered', 'hasstructure', 'tautomerized', 'minimized', 'selected', 'failed' ]:
        mol.ClearProp(prop)
    mol.SetProp('isosmi',isosmi)
    return

def Sane(mol):
    try:
        Chem.SanitizeMol(mol)
        return True
    except:
        return False

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
            bondids=set([ bond.GetIdx() for bond in GetAtomIBonds(match) ])

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


