#!/usr/bin/env python
#-*- coding: utf-8 -*-
from helper import GetScratchDir, Compute3DCoords
from rdkit import Chem

'''
this module include various previous modules that correspond to property
calculation in order to make the code in a cleaner way, including:
    CNDOCalc.py
    CNDOGapCalc.py
    MMGBSACalc.py
    TBGapCalc.py
    CNDOCalc_Old.py
    TBVerticalGapCalc.py
    ZINDOCalc.py
    Docking.py
    GetPartialCharge.py
    Huckel.py
    Gaussian.py
'''

#General Gaussian calculation driver
def GaussianCalc(mol, method, basis=None, maxcyc=None, solvent=False, 
                 optimize=False, memory=None):
    #some preset names and lists
    SemiEmpirical = ("AM1","PM6","PM7","CNDO","INDO","MINDO3","ZINDO")
    AbInitio = ("HF","B3LYP")

    infile = 'gau.com'
    outfile = 'gau.log'
    exe = 'g16' #can be set to g09

    inloc = GetScratchDir() + '/' + infile
    outloc = GetScratchDir() + '/' + outfile

    gauinp = open(inloc, 'w')

    if memory is not None:
        print >> gauinp, '%mem=' + str(memory)
    
    if method.upper() in AbInitio:
        if basis is None:
            raise ValueError, "Basis required for ab initio calculation"
        theory =  method + '/' + basis
    elif method.upper() in SemiEmpirical:
        if basis is not None:
            raise ValueError, "No basis required for semiempirical calculation"
        theory = method

    if solvent: theory = theory + ' ' + 'SCRF=PCM'

    if optimize:
        if maxcyc is not None: 
            theory += ' OPT(MaxCycles=' + str(maxcyc) +')'
        else:
            theory += ' OPT'

    print >> gauinp, '#p ' + theory
    print >> gauinp, '\n GAUSSIAN RUN CALLED BY ACSESS \n\n'
    print >> gauinp, str(Chem.GetFormalCharge(Chem.AddHs(mol))) + ' 1'

    molCoords = Compute3DCoords(mol)

    for line in molCoords:
        coords = '{:<3}{:<10}{:<10}{:<10}'.format(line[3],line[0],line[1],line[2])

        print >> gauinp, coords

    print >> gauinp, '\n\n'

    gauinp.close()

    #run calculation and parse output
    os.system('cd ' + GetScratchDir() + ';' + exe + ' < ' + inloc + ' > ' + outloc)
    #parse using cclib 
    #calcdata = ccopen(outloc).parse()

    return calcdata

