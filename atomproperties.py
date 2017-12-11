#!/usr/bin/env python
#-*- coding: utf-8 -*-
import mprms






############################################################
# Assign polarizabilities to every atom.
# Units are Angstrom^3/2. We assume all valences are filled.
#
# This will work without explicit hydrogens, but it will
# obviously not assign properties to implicit hydrogens.
#
# Values are from Table I of
# Miller and Savchik, JACS 101(24) 7206-7213, 1979.
# dx.doi.org/10.1021/ja00518a014
############################################################
def AssignAtomicPolarizability(mol):
    
    for atom in mol.GetAtoms():
        #get degree for specific molecule
        nBond = atom.GetTotalDegree()

        #Hydrogen (H)
        if atom.GetAtomicNum() == 1:
            atom.SetDoubleProp('polarizability', 0.314)

        #Carbon (C)
        elif atom.GetAtomicNum() == 6:
            if nBond == 4:
                atom.SetDoubleProp('polarizability', 1.294)
            elif nBond == 2:
                GetAtoms.SetDoubleProp('polarizability', 1.393)
            elif nBond == 3:
                if atom.GetNumExplicitHs() + atom.GetNumImplicitHs() > 0:
                    atom.SetDoubleProp('polarizability', 1.428)
                else:
                    '''
                    in this part I employ a different logic than that in
                    previous version of ACSESS and may be consult with Aaron for
                    further details
                    '''
                    cross = True
                    for nbor in atom.GetNeighbors():
                        if not nbor.GetIsAromatic():
                            cross = False
                            break
                    if cross:
                        atom.SetDoubleProp('polarizability', 1.800)
                    else:
                        atom.SetDoubleProp('polarizability', 1.428)

        #Nitrogen (N)
        elif atom.GetAtomicNum() == 7:
            if atom.GetIsAromatic():
                if nBond == 2:
                    atom.SetDoubleProp('polarizability', 1.262)
                else:
                    atom.SetDoubleProp('polarizability', 1.220)
            else:
                if nBond == 1:
                    atom.SetDoubleProp('polarizability', 1.304)
                else:
                    atom.SetDoubleProp('polarizability', 1.435)
        
        #Oxygen (O)
        elif atom.GetAtomicNum() == 8:
            if atom.GetIsAromatic():
                atom.SetDoubleProp('polarizability', 1.099)
            else:
                if nBond == 1:
                    atom.SetDoubleProp('polarizability', 1.216)
                else:
                    atom.SetDoubleProp('polarizability', 1.290)

        #Sulfur (S)
        elif atom.GetAtomicNum() == 16:
            if atom.IsAromatic():
                atom.SetDoubleProp('polarizability', 2.982)
            elif nBond == 2:
                atom.SetDoubleProp('polarizability', 3.496)
            else:
                atom.SetDoubleProp('polarizability', 3.967)

        #Halogens 
        elif atom.GetAtomicNum() == 9:
            atom.SetDoubleProp('polarizability',1.046)
        elif atom.GetAtomicNum() == 15:
            atom.SetDoubleProp('polarizability',3.000)
        elif atom.GetAtomicNum() == 17:
            atom.SetDoubleProp('polarizability',3.130)
        elif atom.GetAtomicNum() == 35:
            atom.SetDoubleProp('polarizability',5.577)
        elif atom.GetAtomicNum() == 53:
            atom.SetDoubleProp('polarizability',8.820)
                
	#Iridium (I)
	#This param value was obtained by fitting the above known tau values
	#In general polarizability increases with atomic number so we used
        #linear fit to get the value
	#This is a crudest approx so could be wrong!
	elif atom.GetAtomicNum() == 77:
            atom.SetData('polarizability',12.77)

        else:
            raise KeyError('No polarizabilities for atomic number'
                           +str(atom.GetAtomicNum()))   

############################################################
# Topological steric effect index (TSEI)
# Cao and Liu, J Chem Inf Comput Sci 44(2), 678-687, 2004.
# dx.doi.org/10.1021/ci034266b

# Covalent radii (in angstroms)
# These are the 2008 values from CSD, except for carbon which is from Cao
#
# This returns numbers a little different than those returned by Marvin.
# If this is to be used for serious business, we should figure out why.
# First check would probably be covalent radius values
# Second would be path termination length
############################################################
def AssignTSEI(mol):
    return




############################################################
# Calculate autocorrelation vector for atomic properties
# Method first described in Moreau and Broto. Nouv. J. Chim.1980, 4, 757-764.
# A more recent reference is dx.doi.org/10.1002/poc.610061008
############################################################
def AutoCorreationVector(mol, maxBonds, props):
    return





############################################################
# Calculation of logp and molar refractivity as given by Crippen in 
# J Chem Inf Comp Sci 1999 39(5) 868-873                   
# dx.doi.org/10.1021/ci990307l                             
############################################################

'''
this part from previous version of ACSESS code is now implemented in 
rdkit package, with same paper as the reference. there is no need for 
us to implement these functions here
'''

