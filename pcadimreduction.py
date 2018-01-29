#!/usr/bin/env python
#-*- coding: utf-8 -*-
from numpy import array, dot, argsort
from numpy.linalg import eig
import numpy as np

from rdkit import Chem

# Container to hold PCA results
class PCADecompose():
    def __init__(self, means, evecs, evals, std_devs=None, forceReal=True):
        self.means = means
        if forceReal:
            self.evecs = np.real(evecs)
            self.evals = np.real(evals)
        else:
            self.evecs = evecs
            self.evals = evals
        self.ndims = len(self.evecs)

        if std_devs is not None:
            self.norm = True
            self.std_devs = std_devs
        else:
            self.norm = False

    def __call__(self, *args, **kwargs):
        return self.Project(*args, **kwargs)

    def NormCenter(self, data):
        dp = data - self.means
        if self.norm:
            dp = dp/self.std_devs
        return dp
    
    def Project(self, data, ndims=None):
        if issubclass(type(data), Chem.Mol):
            coords = array([mol.GetProp('coords') for mol in data])
            meanpoint = self.NormCenter(coords)
            return dot(meanpoint, self.evecs[:, :ndims])
        else:
            if not issubclass(type(data), np.ndarray):
                data = np.array(data)
            if ndims == None:
                ndims = self.ndims
            meanpoint = self.NormCenter(data)
            return dot(meanpoint, self.evecs[:, :ndims])


#X is an numpy array where each row is a data point
#Returns: (loadings, offsets, Projection function) from a PCADecompose object
def PCA(data, nVectors=None, norm=False):
    if type(data) is not np.ndarray:
        try:
            X = np.array([mol.GetProp('coords') for mol in data])
        except AttributeError:
            X = np.array(data)
    else:
        X = data
    nData, nDims = X.shape
    means = X.mean(axis = 0)

    # center data by mean
    X_ave = X - mean
    
    # normalize data by stddev
    if norm:
        std_devs = X.std(axis = 0)
        for i in xrange(len(std_devs)):
            if abs(std_devs[i]) < 1e-10:
                std_devs[i] = 1.0
        
        X_ave = X_ave/std_devs

    # covariance matrix
    Cov = dot(X_ave.T, X_ave) / (1.0 * nData)
    
    # diagonalize covariance matrix
    evals, evecs = eig(Cov)

    # rank eigenvectors, eigenvalues in order
    sorter = argsort(-evals)
    
    # return a PCADecompose object
    kwargs = {}
    if norm: kwargs['std_devs'] = std_devs
    return PCADecompose(means, evecs[:, sorter], evals[:, sorter], **kwargs)
