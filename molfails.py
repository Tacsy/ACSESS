#!/usr/bin/env python
#-*- coding: utf-8 -*-
import os
from rdkit import Chem

#These exceptions are used only for signaling 
#which means no molecular information is carried with them

class NoGeom(Exception):
    pass

class NoConvergence(Exception):
    pass

#Subexceptions of the MolFail exception can carry lots of information
#about a troublesome molecule with them.
#
#If the exception isn't caught, it will dump a lot of information to stdout,
#and write the offending molecule to failMol.smi

class MolFail(Exception):
    def __init__(self, mol = None, message = None):
        self.message = message
        self.mol = mol

    def __str__(self):
        fatalformat = '{0:>4} {1:>8} {2:>17} {3:>11} {4:>13} {5}\n'
        outstring= '-------------- FATAL ERROR --------------\n\n'
        if self.message is not None:
            outstring += self.message+'\n\n'

        if self.mol is not None:
            outstring += 'Type: ' + str(type(self.mol))
            outstring += '\nCanonical smiles: ' + Chem.MolToSmiles(self.mol)
            outstring += '\nIsomeric smiles: ' + Chem.MolToSmiles(self.mol,
                          isomericSmiles = True)
            outstring += '\n\n************Atoms**************\n'
            outstring += fatalformat.format('Idx','Element',
                                            'ExplicitValence',
                                            '#Hydrogens','FormalCharge','Data')
            for atom in self.mol.GetAtoms():
                outstring += fatalformat.format(atom.GetIdx(),
                                                atom.GetAtomicNum(),
                                                atom.GetExplicitValence(),
                                                atom.GetTotalNumHs(
                                                    includeNeighbors = False),
                                                atom.GetFormalCharge(),
                                                atom.GetPropsAsDict(
                                                    includePrivate = False))

                                                
            
            '''
            Need to rewrite when more familiar with bond in rdkit

            outstring+= '\n\n************Bonds**************\n'
            for bond in self.mol.GetBonds():
                outstring += str(bond.GetBgn()) + '-' + str(bond.GetEnd())+\
                     ' order:'+str(bond.GetOrder()) + str(bond.GetData())+'\n'
            outstring+='\nData: '+str(self.mol.GetData())
            outstring+='\nCoords: '+str(self.mol.GetCoords())
            if self.mol.HasData('hasstructure') and \
                   self.mol.GetData('hasstructure') and \
                   SaveStruc:
                geommol=self.mol.GetData('hasstructure')
                outstring+= '\nStructureCoords:'+str(geommol.GetCoords())

            ofs=oe.oemolostream()
            ofs.SetFormat(oe.OEFormat_OEB)
            ofs.open('FailMol.oeb')
            oe.OEWriteMolecule(ofs,self.mol)
            ofs.close()
            outstring+= '\nMolecule written to FailMol.oeb'
            '''

        return outstring


################################
# Other exception possibilities
################################

class MutateFail(MolFail):
    pass

class BadMolecule(MolFail):
    pass

class MutateFatal(MolFail):
    pass

class GroupAlreadyAssigned(MolFail):
    pass

class ChemistryError(MolFail):
    pass
