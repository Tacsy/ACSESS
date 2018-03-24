#!/usr/bin/env python
#-*- coding: utf-8 -*-

import random
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkithelpers import *
from molfails import MutateFail
import mprms

MAXTRY=100
###############################################
# Initialization the module
###############################################

def MutateInit():
    elements = mprms.elements
    halogens = [9,17,35,53]
    global elements, halogens
    
    #do intersection and difference 
    Elements = set(elements) - set(halogens)
    Halogens = set(elements) & set(halogens)
    #change from set to list
    Elements = list(Elements)
    Halogens = list(Halogens)
    

   
###############################################################################
#                            CrossOver Methods                                #
###############################################################################

def GetFragment(mol):
    #print "in get fragment with:", Chem.MolToSmiles(mol)
    tried=[]
    for itry in xrange(MAXTRY):
        #frag=mol.CreateCopy()
        frag=Chem.Mol(mol)
        frag.ClearComputedProps()
        ResetProps( frag )

        
        #Cut an acyclic bond to break molecule in 2
        # get the ring bond ids. this gives a tuple of tuple for every ring
        ringbondids = frag.GetRingInfo().BondRings()
        # flatten this tuple to a 1D set
        ringbondids = list(set([ bond for subset in ringbondids for bond in subset ] ))
        # get all the bonds that are NOT ringbonds or already tried
        bonds = [ bond for bond in mol.GetBonds() if bond.GetIdx() not in tried+ringbondids ]
        if len(bonds)==0: break
        # choose a bond
        delbond=random.choice(bonds)
        tried.append( delbond.GetIdx() )
        # get atom indices of the bond
        i, j = ( delbond.GetBeginAtomIdx(), delbond.GetEndAtomIdx() )

        # convert mol to an editable mol to be able to remove the bond
        #em = Chem.EditableMol(frag)
        em = Chem.RWMol(frag)
        em.RemoveBond(i,j)
        try:
            Chem.SanitizeMol(em)
        except ValueError: continue
        # get new fragments as two new normal Mol objects
        frags = Chem.GetMolFrags(em.GetMol(), asMols=True)
        molfrags = [Chem.MolToSmiles(x,True) for x in frags]
        if len(frags)<2:
            #print "only one frag"
            continue

        #Check to make sure fragments have at least 2 atoms
        if any( frag.GetNumAtoms()<2 for frag in frags):
            #print "not enough atoms on frag"
            continue
        print  molfrags

        # convert to RWMols?
        frags = map(Chem.RWMol, frags)
        return frags

    # If here, we failed to find a good bond to cut
    #print "no nice bond to cut"
    raise MutateFail()


def Crossover(m1, m2):
    #Fragment molecules
    m1fs=GetFragment(m1)
    m2fs=GetFragment(m2)
    #print "in Crossover:", m1fs, m2fs

    #pick fragments for crossover checking:
    # 1. Molecular weight
    Choices=[]
    maxWeight = mprms.maxWeight
    if maxWeight>0:
        w1=[ Descriptors.MolWt(f) for f in m1fs]
        w2=[ Descriptors.MolWt(f) for f in m2fs]
        for i,j in ( (i,j) for i in xrange(2) for j in xrange(2)):
            if w1[i]+w2[j] < maxWeight+50.0: Choices.append( (i,j) )
    # 2. Number of Atoms
    MxAtm=mprms.MxAtm
    if MxAtm>0:
        a1=[f.GetNumAtoms() for f in m1fs]
        a2=[f.GetNumAtoms() for f in m2fs]
        for i,j in ( (i,j) for i in xrange(2) for j in xrange(2)):
            if a1[i]+a2[j] < MxAtm+4: Choices.append( (i,j) )

    if len(Choices)==0:
        raise MutateFail()
    else:
        choice=random.choice(Choices)
        mol1=m1fs[ choice[0] ]
        mol2=m2fs[ choice[1] ]

    # now mol2 has to be connected to mol1
    # but what happens here?
    mol1.SetProp('parent1', m1.GetProp('isosmi'))
    mol2.SetProp('parent2', m2.GetProp('isosmi'))

    '''
    if mol1.HasProp('groups'):
        m1groups=[ groupid for groupid in mol1.GetProp('groups')
                   if len(MolAtIsInGroup(mol1, groupid))>0 ]
        mol1.ClearProp('groups')
    else:
        m1groups=[]
    igrp=1
    if mol2.HasProp('groups'):
        for groupnum in mol2.GetProp('groups'):
            m2groups=mol2.GetProp('groups')
            if len(MolAtIsInGroup(mol2, groupnum))==0: continue
            if groupnum in m1groups: #need to reassign group id
                while igrp in m1groups+m2groups: igrp+=1
                m1groups.append(igrp)
                for atom in MolAtIsInGroup(mol2, groupnum):
                    atom.SetProp('group',igrp)
                for bond in MolBondIsInGroup(mol2, groupnum):
                    bond.SetProp('group',igrp)
                m2groups.append(igrp)
                m2groups.remove(groupnum)
            else:
                m1groups.append(groupnum)
        if len(m2groups)>0: mol2.SetProp('groups',m2groups)
    if len(m1groups)>0: mol1.SetProp('groups',m1groups)
    '''

    ###################################
    #Append molecule 2 to molecule 1. 
    newmol = Chem.CombineMols(mol1, mol2)
    newids = Chem.GetMolFrags(newmol)

    # Now, bond the two fragments to create the child molecule.
    # We choose a random pair of possible atoms to do the attachment;
    possibleA1 = [ atom for atom in GetIAtoms(newids[0], newmol) if EmptyValence(atom)>0 ]
    possibleA2 = [ atom for atom in GetIAtoms(newids[1], newmol) if EmptyValence(atom)>0 ]

    if len(possibleA1)==0 or len(possibleA2)==0:
        #print "no possible atoms!"
        raise MutateFail()
    newmolRW = Chem.RWMol(newmol)
    atom1 = random.choice(possibleA1)
    atom2 = random.choice(possibleA2)
    newmolRW.AddBond( atom1.GetIdx(), atom2.GetIdx(), Chem.BondType.SINGLE )

    print "new mol!", Chem.MolToSmiles(newmolRW)
    return newmolRW.GetMol()


#########################################
# Returns free valence for an atom
# Obviously, there's a problem if it's negative
## After implicit hydrogens are assigned with the valence model,
## this can be replaced by the implicit hydrogen count
# MaxValence={6:4, 7:3, 8:2, 9:1, 15:3, 17:1, 16:2, 35:1, 50:3, 51:2}
#Jos added hydrogen in this list. I don't know if that is a good idea but for 
# openshell a lot of H's will end up in the SMILES strings
MaxValence={1:1, 6:4, 7:3, 8:2, 9:1, 15:3, 17:1, 16:2, 35:1, 50:3, 51:2}
def EmptyValence(atom):

    # Sulfur can have up to 2 double bonds to oxygen
    # (if it's not aromatic)
    if (atom.GetAtomicNum()==16 and atom.HasProp('grouprep')
        ) and (atom.GetProp('grouprep')=='sulfone' or
               atom.GetProp('grouprep')=='sf5'):
            maxv=6

    elif (atom.GetAtomicNum and atom.HasProp('grouprep') and
         atom.GetProp('grouprep')=='nitro'): maxv=4

    else:
        try:
            maxv = MaxValence[atom.GetAtomicNum()]
        except KeyError:
            print "Error in EmptyValence"
            print atom.GetAtomicNum()
            raise


    return maxv - atom.GetExplicitValence()

###############################################################################
#                            Mutation Methods                                 #
###############################################################################

#################################################################
# Flip a bond order
bondorder={ 1:Chem.BondType.SINGLE,
            2:Chem.BondType.DOUBLE,
            3:Chem.BondType.TRIPLE }

#@captureMolExceptions
def FlipBond(mol,bond):
    if bond.HasProp('group'): raise MutateFatal(mol,
        'In-group bond passed to FlipBond')
    oldorder=int(bond.GetBondTypeAsDouble())
    atom1=bond.GetBeginAtom()
    atom2=bond.GetEndAtom()
    full=False
    if EmptyValence(atom1)==0: full=True
    if EmptyValence(atom2)==0: full=True
    if oldorder==1 and full: raise MutateFail()
    elif oldorder==1: add=1
    elif oldorder==3: add=-1
    elif oldorder==2 and full: add=-1
    else: add=random.choice( (-1,1) )
    bond.SetBondType(bondorder[oldorder+add])


#@captureMolExceptions
def SwitchAtom(mol,atom):
    if atom.HasProp('group') and not atom.HasProp('grouprep'):
        raise MutateFatal(mol, 'protected or grouped atom passed to switchatom')

    changed=False

    # # # # # # # # # # # # #
    # from Filters import OptSulfone
    # Special rules for sulfurs that can optionally be sulfones
    #if atom.GetAtomicNum()==16 and len(OptSulfone)>0:
    #    SulfCand=set()
    #    for group in OptSulfone:
    #        for match in mol.GetSubstructMatches(group, True):
    #            for mat in GetIAtoms(match, mol):
    #                if mat.target.IsSulfur(): SulfCand.add(mat.target)
    #    if atom in SulfCand and random.random()>0.4:
    #        if atom.HasProp('grouprep'):
    #            RemoveGroup(mol,atom.GetProp('mygroup'))
    #            return
    #        elif not atom.HasProp('mygroup'):
    #            MakeSulfone(mol,atom)
    #            return

    #If it's a group representative, switching it will delete the
    #rest of the group; since they're hard to create, they're also
    #hard to destroy
    #if atom.HasProp('grouprep':
    #    if random.random()>0.3:
    #        RemoveGroup(mol,atom.GetProp('mygroup'))
    #        changed=True
    #    else:
    #        raise MutateFail()

    # # # # # # # # # # # # # #
    #Special cases for Nitro groups and halogens,
    #which can only be on aromatic rings
    #neighbors=list(atom.GetAtoms())
    #if (  atom.GetExplicitValence()==1
    #      and neighbors[0].IsAromatic() ):

        # Add a halogen, if possible
    #    if (len(Halogens)>0
    #          and random.random()>0.6
    #          and neighbors[0].GetAtomicNum()==6):
    #        CanAddHalogen=True
    #        nextneighbors=neighbors[0].GetAtoms()
    #        for nb in nextneighbors:
    #            if (  nb.GetAtomicNum() != 6
    #                  and nb.GetIdx() != atom.GetIdx()):
    #                CanAddHalogen=False
    #                break
    #        if CanAddHalogen:
    #            atom.SetAtomicNum( random.choice(Halogens))
    #            return

        # If not a halogen, maybe make a nitro group
    #    if random.random() < 0.15:
    #        MakeNitro(mol,atom)
    #        return

    # # # # # # # # #
    # Regular atom switching
    atnum=atom.GetAtomicNum()
    elems= [ el for el in elements if el !=atnum ]

    for itry in range(50):
        cand=random.choice(elems)
        if MaxValence[cand]>=atom.GetExplicitValence():
            atom.SetAtomicNum(cand)
            return

    if not changed:
        raise MutateFail() # if here, mutation failed ...
    return
