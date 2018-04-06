#!/usr/bin/env python

"""
Make maximally diverse library of compounds weight 500 Daltons or less
"""
#metric='mqn' #Use Reymond's Molecular Quantum Numbers metric
metric='autocorr2D'
writeInterval=2 #Write library every 10 iterations (i.e., it10.oeb, it20.oeb...)
aromatic=False

#parameters
nGen=30   #Run for 30 generations
nMut=400  #Create 400 mutants each generation
nCross=50 #Create 50 crossovers each generation
subsetSize=100 #Subset size of 100
MaxPool=1000
mpi=False       # Run MPI (PROGRAM MUST BE RUN WITH MPI)
mutateMPI=False
canonicalTautomer=True #Set each molecule to its canonical tautomerization state
genStruc=True #Generate 3D Structures 
startGenStruc=25 #Start generating structures at iteration 8
startFilters=1 #Start filtering compounds at iteration 0
edgeLen=10    #Preferentially mutate the most diverse 10 compounds
edgeRatio=0.1 #First 10% of mutations will be only of most diverse 10 compounds
maxWeight=500 #Cap molecular weight at 500 Daltons
rseed=787

# mutation frequencies
p_BondFlip = 0.7
p_AtomFlip = 0.7
p_RingAdd = 0.2
p_AddFreq = 0.5  # actual probability: (0.8*0.7)*0.5 = 0.28
p_DelFreq = 0.8  # actual probability: 0.224
p_RingRemove = 0.3 # actual probability: 0.3*0.8 = 0.24
p_AddAroRing = 0.0

#Allowed elements (by atomic number)
elements=[6,7,8,9] #Allowed elements
#extrafilters=['qiu']
