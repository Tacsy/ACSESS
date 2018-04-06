#!/usr/bin/env python
#-*- coding: utf-8 -*-
import numpy as np
from pcadimreduction import PCA 
from scipy.spatial.distance import cdist
import random

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkithelpers import *

import parallel as pl
import mongoserver
import output

metric=None
NormCoords=False

'''
this module include various previous modules that corresponds to chemical space
representation and selection, including:
    CellDiversity.py
    ChemGPS.py
    Distance.py
    Coords.py
    PCAmodule.py
    Kohonen.py
    BuildRandomKoho.py
    MakeSubsetDistance.py
'''

############################################################
#           Functions from Coords.py
############################################################

# Initialize coords
def Init():
    global Coords

    # set coordination system
    if metric == 'mqn':
        from molproperty import CalcMQNs as Coords
    elif metric == 'autocorr':
        from molproperty import MoreauBrotoACVec as Coords
    elif metric == 'autocorr2D':
        from molproperty import AutoCorr2D as Coords
    elif metric is None:
        print 'No coordinate system set!, so Similarity measures will be used!'
    else:
        raise KeyError('Unknown metric specified in parameter file: ' + metric)

# Get coordinates
def SetCoords(mol):
    # set the coordinate of the molecule and return as a np.array object
    # if not, calculate first and then assign it as a property of the molecule
    if mol.HasProp('coords'):
        coord = GetListProp(mol, 'coords')
    else:
        coord = Coords(mol)
        coord = np.nan_to_num(coord)
        SetListProp(mol, 'coords', coord)
    return np.array(coord)

# Drive MPI coordinate calculation
def ScatterCoords(mols):
    
    if not pl.mpi:
        output.StartTimer('COORDS')
        for mol in mols:
            if not mol.HasProp('coords'):
                SetCoords(mol)
        output.EndTimer('COORDS')

        return None

    needCalc = [mol for mol in mols if not mol.HasProp('coords')]
    if len(needCalc) == 0:
        return None

    # Check if values are already computed
    if mprms.UseMongo:
        import mongoserver as mongo

        needSMI = [Chem.MolToSmiles(mol) for mol in needCalc]
        myMongo = mongo.LookupDB(metric, needSMI)
        for v, m in zip(myMongo, needCalc):
            if v is not None:
                SetListProp(m, 'coords', v)

        oldn = len(needCalc)
        needCalc = [mol for mol in mols if not mol.HasProp('coords')]
        print 'Used %d memorized coordinate values' %(oldn - len(needCalc))

        if len(needCalc) == 0:
            return None

    output.StartTimer('COORDS')
    print 'Scattering chemical space coordinate calculation ...', len(needCalc)
    pl.MyTask.SetFunction(MPICoordCalc)
    coords = pl.MyTask.RunMPI(needCalc)

    for coord, mol in zip(coords, needCalc):
        SetListProp(mol, 'coords', coord)

    # Update values in mongo database
    if mprms.UseMongo:
        needSMI = [Chem.MolToSmiles(mol) for mol in needCalc]
        mongo.UpdateDB(metric, {s: GetListProp(m, 'coords') for s, m in zip(needSMI, needCalc)})
        print '%d new memorized coordinate values' %len(needCalc)

    output.EndTimer('COORDS')

# Function for scattering coordinate calculations
#@pl.MPIScatter
def MPICoordCalc(smistring):
    mol = Chem.MolFromSmiles(smistring)

    return SetCoords(mol)



############################################################
#           Functions from Distance.py
############################################################

def GetStdDevs(mols):
    if type(mols) == np.ndarray:
        coords = mols
    else:
        ScatterCoords(mols)
        coords = np.array([SetCoords(mol) for mol in mols])
    std_dev = np.std(coords,axis=0)
    for i in xrange(len(std_dev)):
        if abs(std_dev[i]) < 1e-10: 
            std_dev[i] = 1.0

    return std_dev


def HandleMolCoords(mols,std_dev=None,norm=True):
    #assemble distance vectors
    if type(mols) == np.ndarray:
        print 100
        coords = mols
        passMols = False
    else:
        print 200
        passMols = True
        ScatterCoords(mols)
        coords = np.array([SetCoords(mol) for mol in mols])
    
    #return coordinates according to normalization or not
    if not norm:
        return passMols, coords
    else:
        #normalize coordinates using given std_dev
        if std_dev is not None:
            avgs = np.average(coords,axis=0)
            if not (type(std_dev)==np.ndarray and std_dev.ndim==1):
                std_dev=GetStdDevs(std_dev)

            return passMols, (coords-avgs)/std_dev
        #normalize coordinates by their std_dev
        else:
            try:
                avgs = np.average(coords,axis=0)
            except TypeError:
                print coords
                raise
            std_dev = GetStdDevs(mols)
            return passMols, (coords-avgs)/std_dev


def Maximin(mols,nMol,firstpick=None,startCoords=None,verbose=False):
    '''
    Maximin maximum diversity selection
    mols can be a numpy array containing the coordinates,
    or a list of RDKit molecules

    Use scipy.cdist to handle big array calculations for speed
    Prunes molecules that have no chance of being selected from 
    the set every 1000 steps
    '''
    passMols, coords = HandleMolCoords(mols)
    #if # of mols is smaller than # of mols selected, just keep all of them
    if len(mols) <= nMol:
        if not passMols:
            return range(len(mols))
        else:
            return mols
    #if # of mols is larger
    minDist = np.array([np.infty]*len(mols))
    if startCoords is not None:
        if verbose:
            print 'Calculating starting distances ...'
        for i in xrange(0,len(startCoords),10):
            minDist = np.minimum(minDist,
                            np.min(cdist(startCoords[i:i+10,:],coords),axis=0))
        assert len(minDist) == len(coords)
        minDist = np.array(minDist)
        picks = [np.argmax(minDist)]
    else:
        if firstpick is None:
            firstpick = random.randint(0,len(mols)-1)
        picks = [firstpick]

    lastcoord = coords[picks[0]:picks[0]+1]
    indices = range(coords.shape[0])
    
    #####################
    # select the subset #
    #####################
    for i in xrange(nMol-1):
        dists = cdist(lastcoord,coords)[0]
        minDist = np.minimum(dists,minDist)
        nextpick = minDist.argmax()
        if minDist[nextpick] == 0:
            break
        picks.append(indices[nextpick])
        lastcoord = coords[nextpick:nextpick+1]
        #remove redundant molecules from the array
        if i%1000 == 0:
            oldindices = indices
            oldcoords = coords
            oldminDist = minDist
            indices = []
            coords = []
            minDist = []
            for j in xrange(len(oldminDist)):
                if oldminDist[j]>0.0:
                    indices.append(oldindices[j])
                    coords.append(oldcoords[j])
                    minDist.append(oldminDist[j])
            coords = np.array(coords)
            minDists = np.array(minDist)
            if verbose:
                print i, len(minDist)
                sys.stdout.flush()

    #minimum distance went to 0
    if nMol > 1 and i < nMol-2:
        nleft = nMol - 1 - i
        remaining = set(xrange(len(mols))) - set(picks)
        picks = picks + random.sample(remaining,nleft)

    if passMols:
        return [mols[i] for i in picks]
    else:
        return picks

def SplitSpace(ids,coords):
    '''
    Split sapce along first PCA coordinate,
    return indices of compounds in the two sets
    '''
    coords = np.array([coords[i] for i in ids])
    pcaDecomp = PCA(coords)

    pc1 = np.real(np.dot(coords,pcadecomp.evecs[:,0]))
    avg = sum(pc1)/len(ids)
    std_dev = np.sqrt(sum((pc1-avg)**2)/len(ids))

    iPlus = []
    iMinus = []
    iBound = set()
    for i, pc, coord in zip(xrange(len(pc1)), pc1, coords):
        if pc > avg:
            iPlus.append(ids[i])
        else:
            iMinus.append(ids[i])

        if abs(pc-avg) <= 0.1*std_dev:
            iBound.add(i)

    return iPlus, iMinus, iBound

def PCAMaximin(mols,nMol,nsplit=None,verbose=False):
    '''
    Maximin maximum diversity selection -
    Distribute over multiple nodes by splitting the space using PCA
    '''
    #need to have at least 2 nodes
    if nsplit is not None:
        numNodes = 2**nsplit
    else:
        if not pl.mpi:
            raise NotImplementedError()
        numNodes = pl.MyTask.size
        nsplit = int(np.floor(np.log2(numNodes)))
        numNodes = 2**nsplit

    if numNodes < 2:
        print "Using regular maximin instead of PCAmaximin", numNodes, nsplit
        return Maximin(mols, nMol)

    passMols, coords = HandleMolCoords(mols)

    if len(mols) <= nMol:
        if not passMols:
            return range(len(mols))
        else:
            return mols
    
    ids = [range(len(mols))]

    #split molecules along PCA axis
    boundaryIDs = set()
    for i in xrange(nsplit):
        new_ids = []
        for idset in ids:
            mPlus, mMinus, mBound = SplitSpace(idset,coords)
            new_ids.append(mPlus)
            new_ids.append(mMinus)
            boundaryIDs.update(mBound)
        ids = new_ids

    if verbose:
        print 'PCA/Maximin: split system into regions of sizes: ',
        for idn in ids:
            print len(idn),
        print 

    #to avoid boundary problems, quickly pick a subset of the things
    #on the boundaries, ensuring that picks will not cluster there
    boundaryIDs = list(boundaryIDs)
    cb = np.array([coords[i] for i in boundaryIDs])

    if verbose:
        print 'PCA/Maximin: using fast maximin to select compounds on region bound aries.'
    b_set = FastMaximin(cb,int(nMol*len(boundaryIDs)*1.0/len(mols))) 
    pickset = [boundaryIDs[i] for i in b_set]
    startcoords = np.array([cb[i] for i in b_set])

    if verbose:
        print 'PCA/Maximin: selected', len(pickset), 'out of', len(boundaryIDs),\
              'compounds on region boundaries'
            
    #run maximin on each set
    print 'Scattering PCA-segmented maximin over', numNodes, 'nodes.'
    numToPick = int(1.1*nMol/numNodes)
    toSend = [(np.array([coords[i] for i in idn]), numToPick, startcoords) for idn in minDists]

    pl.MyTask.SetFunction(MPI_PCA_Maximin)
    picks = pl.MyTask.RunMPI(toSend)

    #compile results
    newpicks = []
    newmols = []
    for ids, result in zip(ids, picks):
        for idx in result:
            if passMols:
                newmols.append(mols[ids[idx]])
            else:
                newpicks.append(ids[idx])
    if passMols:
        newmols += [mols[i] for i in pickset]
    else:
        newpicks += pickset

    if passMols:
        return newmols
    else:
        return newpicks

#@pl.MPIScatter
def MPI_PCA_Maximin(MPISEND):
    coords, nMol, startcoords = MPISEND
    picks = Maximin(coords, nMol, startcoords=startcoords)
    return picks


def AveNNDistance(mols,getsqrt=False,std_dev=None):
    '''
    Calculate diversity function (average nearest-neighbor distance)
    '''
    passMols, coords = HandleMolCoords(mols, std_dev=std_dev)

    N = len(coords)
    #determine if distance matrix needs to be chunked or not
    if N > 2000:
        chunksize = 2000*2000/N
        nchunk = N/chunksize + 1
    else:
        chunksize = N
        nchunk = 1

    nndist = 0.0
    for i in xrange(nchunk):
        if i*chunksize >= N:
            break
        dists = cdist(coords, coords[i*chunksize:(i+1)*chunksize])
        for j in xrange(min(chunksize,dists.shape[1])):
            dists[i*chunksize+j,j]= np.inf
        r = np.min(dists,axis=0)
        if getsqrt:
            nndist += sum(r)
        else:
            nndist += np.dot(r,r)

    return nndist/len(mols)

def ScatterAveNNDistance(mols,getsqrt=False):
    '''
    MPI version of calculate diversity function (average nearest-neighbor
    distance)
    '''
    passMols, coords = HandleMolCoords(mols)
    
    N = len(coords)
    #determine if distance matrix needs to be chunked or not
    if N > 2000:
        chunksize = (2000*2000)/N
        nchunk = N/chunksize + 1
    else:
        chunksize = N
        nchunk = 1
    
    nndist = 0.0
    nNode = pl.MyTask.size
    toScatter = [(i*nchunk/nNode, (i+1)*nchunk/nNode, chunksize, getsqrt, coords)
                for i in xrange(nNode)]

    print 'Scattering average nearest neighbor distance calculation ...'
    pl.MyTask.SetFunction(MPIAveNN)

    nndist = sum(pl.MyTask.RunMPI(toScatter))

    return nndist/len(mols)

#@pl.MPIScatter
def MPIAveNN(args):
    startchunk, endchunk, chunksize, getsqrt, coords = args
    nndist = 0.0
    for i in xrange(startchunk,endchunk):
        if i*chunksize >= len(coords):
            break
        dists = cdist(coords, coords[i*chunksize:(i+1)*chunksize])
        for j in xrange(min(chunksize, dists.shape[i])):
            dists[i*chunksize+j,j] = np.inf
        r = np.min(dists,axis=0)
        if getsqrt:
            nndist += sum(r)
        else:
            nndist += np.dot(r,r)
    
    return nndist


