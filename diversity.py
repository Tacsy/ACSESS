#!/usr/bin/env python
#-*- coding: utf-8 -*-

'''
This is the main function of the ACSESS
'''
debug=True
##############################
# Import statements
##############################
import sys
from rdkit import Chem
sys.path.append('.')
import mprms
import init
import drivers as dr

##############################
# Check input
##############################

init.ReadMPRMS() #Does most of the input verification

if mprms.restart: 
    openmode = 'a'
else:
    openmode = 'w'

##############################
# Open output
##############################

filterFile = open('filters.dat', openmode)
convergeFile = open('convergence.dat', openmode)
statsFile = open('stats.dat', openmode)
coordsStdDev = open('stddev.dat', openmode)

statformat='{0:>7} {1:>8} {2:>8} {3:>10} {4:>10} {5:>10} {6:>11} {7:>7} {8:>9}'
convergeformat='{0:>8} {1:>13} {2:>12} {3:>13} {4:>10} {5:>10}'
if not mprms.restart:
    print >> statsFile, statformat.format("#Gen","TooBig","Mutants",
                                      "FailedMut","Undiverse","Filtered",
                                      "Duplicates","Unfit","PoolSize")
    print >> convergeFile, convergeformat.format(
        '#- Round','-- Diversity','-- Max Atoms',
        '-- SubsetSize','-- Filters','-- 3D Geom')


iterhead="\n-------------------- Iteration {0} ----------------\n"

##############################
# Get starting library
##############################

startiter, lib, pool= init.StartLibAndPool(mprms.restart)

###################################################
##########                              ###########
##########           MAIN LOOP          ###########
##########                              ###########
###################################################

for gen in xrange(startiter, mprms.nGen):
    print iterhead.format(gen)

    # PRELOGGING
    if debug:
        print "startiter:", startiter
        print "len lib:", len(lib)
        print "len pool:", len(pool)

    # MUTATIONS
    newlib = dr.DriveMutations(lib)

    # FILTERS
    newlib = dr.DriveFilters(newlib)

    # Select Diverse Set
    pool = dr.ExtendPool(pool, lib, newlib)
    if len(pool)>mprms.subsetSize:
        raise NotImplementedError(0)
        #lib = Maximin(pool)
    else:
        lib = [ mol for mol in pool ]

    if debug:
        with open('mylib','w') as f:
            for mol in lib: f.write(Chem.MolToSmiles(mol)+'\n')

print "DONE"
