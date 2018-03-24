#!/usr/bin/env python
#-*- coding: utf-8 -*-


from rdkit import Chem

############ FROM CANONICAL #############################
def Finalize(mol):
    Chem.SanitizeMol(mol)
    ResetProps(mol)
    Chem.SetAromaticity(mol)
    return

def ResetProps(mol):
    # should 'isosmi' be included?
    isosmi = Chem.MolToSmiles(mol, True)
    for prop in [ 'filtered', 'hasstructure', 'tautomerized', 'minimized', 'selected', 'failed' ]:
        mol.ClearProp(prop)
    mol.SetProp('isosmi',isosmi)
    return

############ NEW RDKIT HELPERS: #########################
def MolAtIsInGroup(mol, groupid):
    requirement = lambda atom:atom.HasProp('group') and atom.GetProp('group')==groupid
    return filter(requirement, mol.GetAtoms())
def MolBondIsInGroup(mol, groupid):
    requirement = lambda bond:bond.HasProp('group') and bond.GetProp('group')==groupid
    return filter(requirement, mol.GetBonds())
CanRemoveAtom=lambda atom:not atom.HasProp('protected') and not atom.HasProp('fixed')

def GetSmallestRingSize(atom):
    return min([ i for i in range(1,20) if atom.IsInRingSize(i) ])

# Two Indices based getters:
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
    return filter(requirement, m.GetBonds())

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
def GetAtomIBonds(atomids, mol):
    bonds=[]
    for bond in mol.GetBonds():
        if all( index in atomids for index in ( bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())):
            bonds.append(bond)
    return bonds

def GetXAtoms(mol, num=6, notprop=None):
    '''gets the number of atoms in mol with atomic num. default carbon, 6 and not has prop notprop
    for every atom it is checked if it has property indicated with the string notprop
    '''
    if notprop:
        requirement = lambda atom:atom.GetAtomicNum()==num and not atom.HasProp(notprop)
    else:
        requirement = lambda atom:atom.GetAtomicNum()==num
    return filter( requirement, mol.GetAtoms())

