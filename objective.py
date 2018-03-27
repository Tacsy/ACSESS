#!/usr/bin/env python

""" Routines for chemical space optimizations. """
import numpy as np
import sys
import os
import mprms
from rdkit import Chem

minimize=False #Minimizing or maximizing the objective function?
compGeom=False #Does objective function require 3d structure?
TargetScore=None
TargetScoreOss=None
NeighborhoodMaximin=True
FixedCutoff=False
InitialCutoff=0.0
TakeFittest=0
SaturateAt=0.8
NeighborhoodFactor=0.75
array_or=np.vectorize(lambda x,y:x or y)
sumhelper=None
qscatter=False
CINDES_interface=False
CINDES_run=None
debug=False
OBDatabase=None

#Initialize module
def Init():
    global minsign, CINDES_interface, CINDES_run
    if minimize: minsign = 1.0
    else:        minsign =-1.0
    if compGeom and not SaveStruc:
        raise ValueError('Objective function requires geometry, ' +
                         'but SaveStruc is false.')

    if CINDES_interface:
        import QCindes as qc
        CINDES_run = qc.Init()
        print CINDES_run
    else:
        raise NotImplementedError('only CINDES objectives are supported currently')

    if debug:
       if restart: mode='a'
       else: mode='w'
    return


###############################################################
# Output setup
fitnessformat='{0:>7} {1:>8} {2:>8} {3:>10} {4:>10} {5:>10} {6:>10}'
    
#Drive computation of Objective Function
def ComputeObjectives(mols_tocalc, gen=0):
    if CINDES_interface:
        import QCindes as qc
        global CINDES_run
        print 'calculating via CINDES program'
        qc.calculate(mols_tocalc, CINDES_run, gen=gen)
    else: #Serial
        raise NotImplementedError('only CINDES objectives are supported currently')
                
#Make sure all objective function values are calculated
def UpdateObjective(mylib, gen=0):
    mols_tocalc=[ mol for mol in mylib if not mol.HasProp('Objective') ]

    print 'Computing Objective...'
    ComputeObjectives(mols_tocalc, gen)
    
    newvals={Chem.MolToSmiles(m, True):m.GetProp('Objective') for m in mols_tocalc}
    return

######################################################################
# Optimization routine: rank all compounds in terms of objective function
# Remove compounds below the threshold
# Returns both picked compounds (as library) and unpicked compounds (as pool)
def EvaluateObjective(totallib, gen):
    nGen     = mprms.nGen
    NumIn    = len(totallib)

    # 1. do actual objective evaluation
    UpdateObjective(totallib, gen)
    # Sort compounds by fitness
    totallib.sort(key=lambda x: x.GetDoubleProp('Objective'),reverse=not minimize)
         
    # 2. Cut compounds below cutoff
    if gen>SaturateAt*nGen:
        cutoff=TargetScore
    else:
        cutoff=InitialCutoff + (TargetScore-InitialCutoff)*( gen*1.0)/(nGen*SaturateAt)
    newpool = [ mol for mol in totallib if mol.GetDoubleProp('Objective')*minsign <= cutoff*minsign ]
    NumOut=len(newpool)
    if NumOut==0:
        raise ValueError('No compounds left after cutoff!')
    print 'Removed',NumIn-NumOut,'compounds using objective function threshold'

    return newpool

def SelectFittest(newpool, SubsetSize):
    if len(newpool)>SubsetSize+TakeFittest:
        mostfit = newpool[0:TakeFittest]
        pool = newpool[TakeFittest:]   ###at this point pool is updated so that only compounds satisfying fitness stay
        nSwap=0

        dt.ScatterCoords( pool )

        coords=np.array( [ mol.GetProp('coords') for mol in pool ] )
        scores=np.ma.array( [mol.GetDoubleProp('Objective') for mol in pool] )
	if sumhelper is None: sumhelper=np.ones(len(coords[0]))

        #Allows you to do pure Fitness based selection
        if SubsetSize<=0: return mostfit,[],newpool

        #Select a sample subset
        picks=dt.FastMaximin(coords,SubsetSize)
        
        templib=[ pool[i] for i in picks ]
        if dt.NormCoords:
            coords=coords/dt.GetStdDevs( templib )

        oput.StartTimer('NN DIST CALC')
        AveDistSqr=dt.AveNNDistance(templib)
        oput.EndTimer('NN DIST CALC')

        print 'Average libary score before optimization: ',sum(
            m.GetDoubleProp('Objective') for m in templib )/(
            1.0*len(templib)),'(diversity:',AveDistSqr,')'
        
        ############################ NEIGHBORHOOD MAXIMIN ##################
        #Discard the original subset; instead, pick the BEST SCORING COMPOUND
        #within the neighborhood of each compound        
        pickmask=np.zeros( len(pool), dtype=np.bool)
        for i in picks: pickmask[i]=True
        targetmask = np.ma.getmask( np.ma.masked_greater(
                    scores*minsign , TargetScore*minsign )) 

        oput.StartTimer("OBJECTIVE MXMN")
        print 'Optimizing library ...',
        newlib=[]

        ######### MAIN LOOP #########
        for ipick in picks:
            myscore=pool[ipick].GetDoubleProp('Objective')
            #Skip compounds already at target
            if myscore*minsign <= TargetScore*minsign:
                newlib.append( pool[ipick] )
                continue

            #Mask compounds outside of current neighborhood
            #or that have already been picked
            diffvec = coords - coords[ipick]
            distsqr = np.dot(diffvec*diffvec,sumhelper)
            neighbormask = np.ma.getmask( np.ma.masked_greater(
                                         distsqr,NeighborhoodFactor*AveDistSqr))
            neighbor_pick_mask=array_or(neighbormask,pickmask)
            if debug: print len(neighbormask)-sum(neighbormask),

            #If any compounds in the neighborhood hit the target, pick
            #the closest one
            distsqr=np.ma.masked_array(distsqr,mask=array_or(targetmask,
                                                          neighbor_pick_mask))
            if distsqr.count()>0:
                inewpick=np.argmin(distsqr)
                newlib.append( pool[inewpick] )
                pickmask[ inewpick ] = True
                pickmask[ ipick ] = False
                nSwap+=1
                continue

            #If there is no compound hitting the target, pick the
            #best one in the neigbhorhood
            scores.mask = neighbor_pick_mask
            inewpick = np.argmin(minsign*scores)
            if scores[inewpick]*minsign < myscore*minsign:
                newlib.append( pool[inewpick] )
                pickmask[ inewpick ] = True
                pickmask[ ipick ] = False
                nSwap+=1
            else:
                newlib.append(pool[ipick])

        #####
        #Done with optimizing maximin
        newlib.sort(key=lambda x: x.GetDoubleProp('Objective'),reverse=not minimize)
        print 'swapped',nSwap,'/',len(newlib),'compounds'
        print 'Average libary score after optimization: ',
        
        oput.EndTimer('OBJECTIVE MXMN')
    else:
       newlib=newpool
       mostfit=[]

    #Print statisticals
    fvals1= np.array( [ mol.GetDoubleProp('Objective') for mol in newlib+mostfit] )
    print "fvals1:", fvals1

    print>>simstats,fitnessformat.format(gen,NumIn,NumOut,
                                         np.average(fvals1),
                                         np.min(fvals1),
                                         np.max(fvals1),0.0)
    
    simstats.flush()

    newlib += mostfit
    newlib.sort( key=lambda x:x.GetDoubleProp('Objective'), reverse=not ob.minimize)
    return newlib, newpool

