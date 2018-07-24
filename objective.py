#!/usr/bin/env python
""" Routines for chemical space optimizations. """
import numpy as np
import sys
import os
import mprms
from rdkit import Chem
from rdkithelpers import *
from drivers import RemoveDuplicates
import distance
import output

# fixed global variables
_array_or = np.vectorize(lambda x, y: x or y)
_sumhelper = None

# global variables
minimize = False  #Minimizing or maximizing the objective function?
compGeom = False  #Does objective function require 3d structure?
TargetScore = None
FixedCutoff = False
InitialCutoff = 0.0
TakeFittest = 1
SaturateAt = 0.8
NeighborhoodMaximin=True
NeighborhoodFactor = 0.75
fitnessfunction = None
CINDES_interface = False
debug = False

#Initialize module
def Init():
    global minsign, CINDES_interface, CINDES_run
    if minimize: minsign = 1.0
    else: minsign = -1.0
    if compGeom and not SaveStruc:
        raise ValueError('Objective function requires geometry, ' +
                         'but SaveStruc is false.')

    if CINDES_interface:
        global qc
        import QCindes as qc
    elif callable(fitnessfunction):
        pass
    else:
        raise NotImplementedError(
            'only CINDES objectives are supported currently')

    if debug:
        if restart: mode = 'a'
        else: mode = 'w'
    return

#Drive computation of Objective Function
def ComputeObjectives(mols_tocalc, gen=0):
    if CINDES_interface:
        print 'calculating via CINDES program'
        qc.calculate(mols_tocalc, gen=gen)
    elif callable(fitnessfunction):
        for mol in mols_tocalc:
            value = fitnessfunction(mol)
            mol.SetDoubleProp('Objective', float(value))
    else:  #Serial
        raise NotImplementedError(
            'only CINDES objectives are supported currently')


#Make sure all objective function values are calculated
def UpdateObjective(mylib, gen=0):
    print "UpdateObjective len(mylib):{}".format(len(mylib)),
    mols_tocalc = [mol for mol in mylib if not mol.HasProp('Objective')]
    print "n mols_tocalc:", len(mols_tocalc)

    print 'Computing Objective...'
    ComputeObjectives(mols_tocalc, gen)

    newvals = {
        Chem.MolToSmiles(m, True): m.GetProp('Objective')
        for m in mols_tocalc
    }
    return


######################################################################
# Optimization routine: rank all compounds in terms of objective function
# Remove compounds below the threshold
# Returns both picked compounds (as library) and unpicked compounds (as pool)
def EvaluateObjective(totallib, gen):
    nGen = mprms.nGen
    NumIn = len(totallib)

    # 0. Remove Duplicates
    totallib = RemoveDuplicates(totallib)

    # 1. do actual objective evaluation
    UpdateObjective(totallib, gen)

    # 2. Sort compounds by fitness
    totallib.sort(
        key=lambda x: x.GetDoubleProp('Objective'), reverse=not minimize)

    # 3. Cut compounds below cutoff
    if gen > SaturateAt * nGen:
        cutoff = TargetScore
    else:
        cutoff = InitialCutoff + (TargetScore - InitialCutoff) * (
            gen * 1.0) / (nGen * SaturateAt)
    newpool = [ mol for mol in totallib
                if mol.GetDoubleProp('Objective') * minsign <= cutoff * minsign
              ]

    # 4. Log results 
    NumOut = len(newpool)
    print "current cutoff:", cutoff, 'minsign:', minsign
    if NumOut == 0:
        raise ValueError('No compounds left after cutoff!')
    print 'Removed', NumIn - NumOut, 'compounds using objective function threshold'

    output.obstats.clear()
    output.obstats['NumIn']=NumIn
    output.obstats['NumOut']=NumOut

    return newpool


def SelectFittest(newpool, subsetSize, gen):

    # only do selection when there are enough molecules
    if len(newpool) > subsetSize + TakeFittest:
        # 1. we split newpool in two parts: mostfit and pool
        # Note: we assume that newpool is sorted based on the objective value.
        # we select the mostfit automatically and do a diversity selection of pool
        mostfit = newpool[0:TakeFittest]
        pool = newpool[TakeFittest:]
        
        if NeighborhoodMaximin:
            newlib  = NeighborhoodMaximin(pool, subsetSize)
        elif False:
            # here I could call a function optimizing diversity and the objective value by optimizing
            # a fitnessfunction: fitness(x) = c1 * Div(x) + c2 * Prop(x) with
            # - c1 + c2 = 1
            # - Div and Prop are normalized: Take care of outliers(inf / -inf values!) / negative+positive values
            pass
        else:
            raise NotImplementedError('currently only NeighborhoodMaximin is implemented')
    else:
        newlib = newpool
        mostfit = []

    #Print statisticals
    fvals = np.array(
        [mol.GetDoubleProp('Objective') for mol in newlib + mostfit])
    print "fvals:", fvals
    output.obstats.update({
        'gen': gen,
        'fvals': fvals,
    })
    output.PrintObjectiveStat()

    newlib += mostfit
    newlib.sort(
        key=lambda x: x.GetDoubleProp('Objective'), reverse=not minimize)
    return newlib, newpool


def NeighborhoodMaximin(pool, subsetSize):
    global _sumhelper
    nSwap = 0

    # 2. we assign molecular coordinates to the molecules in pool
    distance.ScatterCoords(pool)
    coords = np.array([GetListProp(mol, 'coords') for mol in pool])

    # 3. Do some initializations:
    scores = np.ma.array([mol.GetDoubleProp('Objective') for mol in pool])
    if _sumhelper is None:
        _sumhelper = np.ones(len(coords[0]))

    # 4 Calculate the diversity measure for the pure diversity subset
    # 4.1 First Select a pure diversity based sample subset
    picks = distance.Maximin(coords, subsetSize)
    templib = [pool[i] for i in picks]
    if distance.NormCoords:
        coords = coords / distance.GetStdDevs(templib)
    # 4.2 calculate the average distance
    output.StartTimer('NN DIST CALC')
    AveDistSqr = distance.AveNNDistance(templib)
    output.EndTimer('NN DIST CALC')
    print 'average diversity value of pure diversity subset:', AveDistSqr
    # 4.3 calculate the average objective value
    aveobj = sum( m.GetDoubleProp('Objective') for m in templib) / (len(templib)) 
    print 'Average objective value of pure diversity  subset:', aveobj

    ############################ NEIGHBORHOOD MAXIMIN ##################
    #Discard the original subset; instead, pick the BEST SCORING COMPOUND
    #within the neighborhood of each compound

    # make a masked array so we don't pick already picked ones 
    pickmask = np.zeros(len(pool), dtype=np.bool)
    for i in picks:
        pickmask[i] = True
    # mask every value larger than TargetScore.
    targetmask = np.ma.getmask(
        np.ma.masked_greater(scores * minsign, TargetScore * minsign))

    output.StartTimer("OBJECTIVE MXMN")
    print 'Optimizing library ...',
    newlib = []

    ######### MAIN LOOP #########
    for ipick in picks:
        myscore = pool[ipick].GetDoubleProp('Objective')
        #Skip compounds already at target
        if myscore * minsign <= TargetScore * minsign:
            newlib.append(pool[ipick])
            continue

        #Mask compounds outside of current neighborhood
        #or that have already been picked
        # 1D difference vector ( positives and negatives )
        diffvec = coords - coords[ipick]
        # total sum of squared distances
        distsqr = np.dot(diffvec * diffvec, _sumhelper)
        # get boolean array with False's for every distancevalue greater than 0.75*AveDist 
        neighbormask = np.ma.getmask(
            np.ma.masked_greater(distsqr, NeighborhoodFactor * AveDistSqr))
        # combine maks with the pickmask
        neighbor_pick_mask = _array_or(neighbormask, pickmask)
        if debug: print len(neighbormask) - sum(neighbormask),

        #If any compounds in the neighborhood hit the target, pick
        #the closest one
        # get dist's of mols not already picked, in neighborhood
        # and objective value above cutoff value
        distsqr = np.ma.masked_array(
            distsqr, mask=_array_or(targetmask, neighbor_pick_mask))
        if distsqr.count() > 0:
            # change pick for best pick:
            inewpick = np.argmin(distsqr)
            newlib.append(pool[inewpick])
            # adjust pickmask
            pickmask[inewpick] = True
            pickmask[ipick] = False
            nSwap += 1
            continue

        #If there is no compound hitting the target, pick the
        #best one in the neigbhorhood
        # get only scores in neighborhood not already picked
        scores.mask = neighbor_pick_mask
        # and get the best value even if not fullfilling cutoff
        inewpick = np.argmin(minsign * scores)
        # if inewpick is different from the current ipick:
        if scores[inewpick] * minsign < myscore * minsign:
            newlib.append(pool[inewpick])
            pickmask[inewpick] = True
            pickmask[ipick] = False
            nSwap += 1
        else:
            newlib.append(pool[ipick])

    #####
    #Done with optimizing maximin
    newlib.sort(
        key=lambda x: x.GetDoubleProp('Objective'), reverse=not minimize)
    print 'swapped', nSwap, '/', len(newlib), 'compounds'
    print 'Average libary score after optimization: ',

    output.EndTimer('OBJECTIVE MXMN')
    output.obstats['nSwap']=nSwap
    return newlib
