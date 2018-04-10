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
    init.ReadMPRMS() #Does most of the input verification
    ##############################
    # 2. Open output
    ##############################
    global openmode
    if hasattr(mprms, 'restart') and mprms.restart: 
        openmode = 'a'
    else:
        mprms.restart=False
        openmode = 'w'
    InitiateFileHandlers(openmode)
    ##############################
    # 3. Get starting library
    ##############################
    global startiter, lib, pool
    startiter, lib, pool= init.StartLibAndPool(mprms.restart)
    return

def InitiateFileHandlers(openmode):
    global convergeFile, statsFile, coordsStdDev, statformat, convergeformat
    dr.filterFile = open('filters.dat', openmode)
    convergeFile = open('convergence.dat', openmode)
    output.statsFile = open('stats.dat', openmode)
    objective.simstats=open('Fitness.dat','a')
    coordsStdDev = open('stddev.dat', openmode)
    statformat='{0:>7} {1:>8} {2:>8} {3:>10} {4:>10} {5:>10} {6:>11} {7:>7} {8:>9}'
    convergeformat='{0:>8} {1:>13} {2:>12} {3:>13} {4:>10} {5:>10}'
    if not mprms.restart:
        print >> output.statsFile, statformat.format("#Gen","TooBig","Mutants",
                                      "FailedMut","Undiverse","Filtered",
                                      "Duplicates","Unfit","PoolSize")
        print >> convergeFile, convergeformat.format(
            '#- Round','-- Diversity','-- Max Atoms',
            '-- SubsetSize','-- Filters','-- 3D Geom')
    return

#initiate()


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
            lib, pool = objective.SelectFittest(pool, mprms.subsetSize, gen)
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

