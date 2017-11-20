#!/usr/bin/env python
#-*- coding: utf-8 -*-




##############################
# Check input
##############################

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





###################################################
##########                              ###########
##########           MAIN LOOP          ###########
##########                              ###########
###################################################

for gen in xrange(startiter, mprms.nGens):

