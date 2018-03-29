#!/usr/bin/env python
import oechem as oe
import numpy as np
from Helpers import *
SARescale=False

#Needs to be called first to read in database fragfile
def SAInit(fragfile='/home/jlt84/RD-ACSESS/ACSESS/Filters/SAFrag2.P'):
    global frags
    frags=Depickle(fragfile)
    print 'Reading fragment library ...'
    for key in frags:
        frags[key]=np.log10(frags[key])


#Get complexity score
isterm=lambda atom: (not atom.GetAtomicNum()==1) and atom.GetDegree()==1
spirosearch=Chem.MolFromSmarts('*@*(@*)(@*)@*')
def ComplexityScore(mol,termPenalty=False,debug=False):
    nrings,sa,sb=SSSR(mol)

    nspiro=len( list( mol.GetSubstructMatches(spirosearch) ) )

    nmacro=sum(nrings[6:])

    natoms=mol.NumAtoms()
    sizePenalty=natoms**1.005-natoms

    if termPenalty==1:
        nTerm=len(filter(isterm, mol.GetAtoms()))
        sizePenalty += nTerm**1.035-nTerm
    elif termPenalty==2:
        nTerm=(1.0*len(filter(isterm, mol))/ mol.GetNumHeavyAtoms()
        sizePenalty += 16.0 * nTerm**3.0

    nstereo=oe.OECount( mol, oe.OEIsChiralAtom() )
    if debug: print sa,nspiro,nmacro,sizePenalty,nstereo

    return np.log10( nmacro+1 ) + sizePenalty + np.log10(nstereo+1) +\
        np.log10( sa-nspiro+1 ) + np.log10( nspiro+1)
        

#Get fragment score
def FragScore(mol,warn=False):
    global frags
    myscore=0.0
    for nfrag,fragment in enumerate(ExtendedConnectivityFragments(mol)):
        if frags.has_key(fragment):
            myscore+=frags[fragment]
        else:
            if warn: print 'novel fragment: '+fragment
            myscore-=0.1

    nfrag+=1
    return myscore/(1.0*nfrag)


#Return total score
def SAScore(mol,termPenalty=False):
    sascore = FragScore(mol)-ComplexityScore(mol,termPenalty)
    if SARescale:
        # need to transform "raw" value into scale between 1 and 10
        samin = -4.0
        samax = 2.5
        sascore = 11. - (sascore - samin + 1) / (samax - samin) * 9.
        # smooth the 10-end
        if sascore > 8.:
            sascore = 8. + np.log(sascore + 1. - 9.)
        if sascore > 10.:
            sascore = 10.0
        elif sascore < 1.:
            sascore = 1.0
    return sascore
        



############################
# Routines for getting fragments
AtomSetFunctor = lambda ids: oe.PyAtomPredicate( lambda a: a.GetIdx() in ids )

atommap=oe.OEAtomArray(500)
subsetmol=oe.OEGraphMol()
isheavy=oe.IsHeavy()

def ExtendedConnectivityFragments(m,levels=(1,4),
                                  getfrag=False,neutral=True):
    mol=m.CreateCopy()

    if neutral:
        for atom in mol.GetAtoms(): atom.SetFormalCharge(0)
        oe.OEAssignMDLHydrogens(mol)

    oe.OEAddExplicitHydrogens(mol)
    neighbors={}

    #Build neighbor list for each atom
    for atom in mol.GetAtoms():
        neighbors[ atom.GetIdx() ]=set( a.GetIdx() for a in atom.GetAtoms() )


    #Enumerate all fragments for each atom
    for atom in mol.GetAtoms(oe.OEIsHeavy()):
        myid=atom.GetIdx()
        myfrag=set( (myid , ) )

        #Add next level of bonded atoms
        for ilevel in xrange(*levels):
            newset=set()
            for nID in myfrag:
                newset.update( neighbors[nID] )
            newset=newset-myfrag
            if len(newset)==0: break
            myfrag.update(newset)
            subsetmol.Clear()

            #Isolate new fragment
            oe.OESubsetMol( subsetmol, mol, AtomSetFunctor( myfrag) ,
                            False, False, atommap)


            #Make "star" atoms representing edge connectivity
            atommap[ myid ].SetData('center',True)
            for aid in newset:
                dest=atommap[aid]
                if isheavy(dest):
                    if dest.IsAromatic():
                        dest.SetAtomicNum(18)
                    else:
                        dest.SetAtomicNum(89)
                    
            if getfrag: yield subsetmol
            else: yield oe.OECreateSmiString(subsetmol ,
                                             oe.OEOFlavor_CAN_Hydrogens |
                                             oe.OEOFlavor_CAN_ExtBonds |
                                             oe.OEOFlavor_CAN_Canonical )




#Run as a program with -getfrag to analyze a database
import sys
if __name__=='__main__':


    if len(sys.argv)==1:
        print 'USAGE: SAScore.py [molecule file]'

    elif sys.argv[1]=='-getfrag':
        from Classes import Counter
        from Helpers import basename,Enpickle

        print 'Molecules | unique fragments | total fragments'
        frags=Counter( [] )
        for i,mol in enumerate( oe.oemolistream( sys.argv[2] ).GetOEGraphMols()):
            frags.AddList( ExtendedConnectivityFragments( mol ) )

            if i%250==0: print '\r',i,len(frags),sum(frags.values() ),
            sys.stdout.flush()

        Enpickle(frags,basename(sys.argv[2])+'.P.gz')
           

    else:
        fragfile='/home/dbchem/av65/newgarlics/SAFrag.P.gz'
        
        infile=oe.oemolistream(sys.argv[1])

        print 'Reading fragment database from '+fragfile+'...',
        sys.stdout.flush()
        SAInit(fragfile)

        print 'done.'

        print 'Reading molecules from file :',sys.argv[1]
        print 'Writing results to file:     ',basename( sys.argv[1]
                                                   )+'_SAScore.smi'
        outfile=oe.oemolostream( basename( sys.argv[1] )+'_SAScore.smi')

        sys.stdout.flush()
        
        for imol,mol in enumerate(infile.GetOEGraphMols()):
            mol.SetTitle( str(SAScore(mol)) )
            oe.OEWriteMolecule(outfile,mol)
            if imol%20==0:
                print '\rCalculated',imol+1,'scores ...',
                sys.stdout.flush()

        print 'Finished.'

        outfile.close()
