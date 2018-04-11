#!/usr/bin/env python
#-*- coding: utf-8 -*-

'''
This is the main function of the ACSESS
'''
debug=True
##############################
# Import statements
##############################
import sys, random
from rdkit import Chem
sys.path.append('.')
import mprms
import init
import drivers as dr
import output
from output import stats
import objective
from helpers import DumpMols, FinishSelection
from distance import AveNNDistance
from similarity import NNSimilarity
iterhead="\n-------------------- Iteration {0} ----------------\n"

def initiate():
    ##############################
    # 1. Check input
    ##############################
    init.Initialize() #Does most of the input verification
    ##############################
    # 2. Get starting library
    ##############################
    global startiter, lib, pool
    startiter, lib, pool= init.StartLibAndPool(mprms.restart)
    return

def evolve():
    global startiter, lib, pool, iterhead

    ###################################################
    ##########                              ###########
    ##########           MAIN LOOP          ###########
    ##########                              ###########
    ###################################################
    for gen in xrange(startiter, mprms.nGen):
        # 0. Decide workflow for current iteration:
        Tautomerizing, Filtering, GenStrucs= dr.SetIterationWorkflow(gen)
 
        # 1. PRELOGGING
        print iterhead.format(gen)
        stats.update({'gen':gen, 'nPool':len(pool), 'nLib':len(lib)})
        if debug:
            print "startiter:", startiter
            print "len lib:", len(lib)
            print "len pool:", len(pool)
 
        # 2.MUTATIONS AND CROSSOVERS
        newlib = dr.DriveMutations(lib)
 
        # 3. FILTERS
        newlib = dr.DriveFilters(newlib, Filtering, GenStrucs)

        # 3b. When necessary FilterPool is switching is just set on
        pool   = dr.DrivePoolFilters(pool, Filtering, GenStrucs, Tautomerizing, gen)

        # 4. OBJECTIVE
        if mprms.optimize:
            print "acsess.py len(newlib, pool):{} {}".format(len(newlib), len(pool))
            pool = objective.EvaluateObjective(newlib + pool, gen)
        else:
            pool = dr.ExtendPool(pool, lib, newlib)
 
        # 5. SELECTION
        if mprms.optimize:
            oldN = len(pool)
            lib, pool = objective.SelectFittest(pool, mprms.subsetSize, gen)
            stats['nUnFit']=oldN - len(pool)
        elif len(pool)>mprms.subsetSize:
            #lib = random.sample(pool, mprms.subsetSize)
            lib = dr.DriveSelection(pool, mprms.subsetSize)
        else:
            lib = [ mol for mol in pool ]
        FinishSelection(lib)
        if len(lib)==0: raise RuntimeError('no molecules left')

        # 6. DIVERSITY IS:
        #siml = NNSimilarity(lib)
        if hasattr(mprms, 'metric'):
            siml = AveNNDistance(lib)
        else:
            siml = NNSimilarity(lib)
        print '\nLIBRARY DIVERSITY: ',siml
 
        # 7. POSTLOGGING
        if debug:
            with open('mylib','w') as f:
                for mol in lib: f.write(Chem.MolToSmiles(mol)+'\n')
        if gen % mprms.writeInterval==0 or gen==mprms.nGen-1:
            DumpMols(lib, gen)
            DumpMols(pool)
        stats['diversity']=siml
        output.PrintTimings()
        output.PrintStat()

    output.PrintTotalTimings()
    print "DONE"
    return

if __name__=="__main__":
    class RunACSESS(object):
        def __init__(self):
            initiate()
        def evolve(self):
            evolve()
    run = RunACSESS()
    run.evolve()

