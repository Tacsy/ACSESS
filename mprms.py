#!/usr/bin/env python
#-*- coding: utf-8 -*-

"""
This is the parameter file containing essential parameters for ACSESS to parse
"""

##############
### input ####
##############

# input seedfile, can be filename, or None
seedFile = 'mols.smi'

##############
### output ###
##############

# write library every 2 iterations
writeInterval = 2

##############
#### main ####
##############

restart = False
# run for 1000 generations
nGen = 25
# create 40 mutants each generation
nMut = 40
# create 50 crossovers each generation
nCross = 50
# subset size after each generation of 70
subsetSize = 70
# start generating structures at interation 100
startGenStruc = 100
# start filetering compounds at interation 100
startFilter = 100
# preferentially mutate the most diverse 10 compounds
edgeLen = 10
# first 10% of mutations will be only of most diverse 10 compounds
edgeRatio = 0.1
# cap molecular weight at 500 Daltons
maxWeight = 500

##############
### mutate ###
##############

# mutation frequencies
p_BondFlip = 0.7
p_AtomFlip = 0.7
p_RingAdd = 0.2
p_AddFreq = 0.5  # actual probability: (0.8*0.7)*0.5 = 0.28
p_DelFreq = 0.8  # actual probability: 0.224
p_RingRemove = 0.3 # actual probability: 0.3*0.8 = 0.24

# allowed elements during mutation (by atomic number)
elements = [6,7,8,9,16,17] # C, N, O, F, S, Cl, specifically

##############
#### spec ####
##############

# run MPI (program must be run with MPI)
mpi = False
# do mutations with MPI (only available if mpi = True)
mutateMPI = False
# set every molecule to its canonical tautomerization state
canonicalTautomer = True
# generate 3D structures
genStruc = True
# open filters
qiuFilter = False
