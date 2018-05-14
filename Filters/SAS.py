#!/usr/bin/env python
#-*- coding: utf-8 -*-

import os, math
from rdkit.six.moves import cPickle
from rdkit.Chem import AllChem
from rdkit.six import iteritems
from rdkit import Chem
import copy

SARescale = False
_fscores = None


def ReadFragScores(name='fpscores'):
    import gzip
    global _fscores
    #generate the full path filename
    if name == "fpscores":
        name = os.path.join(os.path.dirname(__file__), name)

    _fscores = cPickle.load(gzip.open('%s.pkl.gz' % name))
    outDict = {}
    for i in _fscores:
        for j in range(1, len(i)):
            outDict[i[j]] = float(i[0])
    _fscores = outDict


def NumBridgeheadsAndSpiro(mol, ri=None):
    nSpiro = AllChem.CalcNumSpiroAtoms(mol)
    nBridgehead = AllChem.CalcNumBridgeheadAtoms(mol)

    return nBridgehead, nSpiro


def CalcSAScore(rmol):
    if _fscores is None:
        ReadFragScores()
    mol = copy.deepcopy(rmol)
    #Chem.SanitizeMol(mol) # gives crashes!

    #fragment score
    fp = AllChem.GetMorganFingerprint(
        mol, 2)  #<- 2 is the *radius* of the circular fingerprint
    fps = fp.GetNonzeroElements()
    score1 = 0.0
    nf = 0
    for bitId, v in iteritems(fps):
        nf += v
        sfp = bitId
        score1 += _fscores.get(sfp, -4) * v
    score1 /= nf

    #features score
    nAtoms = mol.GetNumAtoms()
    nChiralCenters = len(
        Chem.FindMolChiralCenters(mol, includeUnassigned=True))
    ri = mol.GetRingInfo()
    nBridgehead, nSpiro = NumBridgeheadsAndSpiro(mol, ri)
    nMacrocycles = 0
    for x in ri.AtomRings():
        if len(x) > 8:
            nMacrocycles += 1

    sizePenalty = nAtoms**1.005 - nAtoms
    stereoPenalty = math.log10(nChiralCenters + 1)
    spiroPenalty = math.log10(nSpiro + 1)
    bridgePenalty = math.log10(nBridgehead + 1)
    macrocyclePenalty = 0.0
    # -----------------------------
    # This differs from the paper, which defines:
    #   macrocyclePenalty = math.log10(nMacrocycles+1)
    # This form generates better results when 2 or more macrocycles are present
    if nMacrocycles > 0:
        macrocyclePenalty = math.log10(2)

    score2 = 0.0 - sizePenalty - stereoPenalty - spiroPenalty - bridgePenalty - macrocyclePenalty

    # correction for the fingerprint density
    # not in the original publication
    # to make highly symmetrical molecules easier to synthesize
    score3 = 0.0
    if nAtoms > len(fps):
        score3 = math.log(float(nAtoms) / len(fps)) * 0.5

    sascore = score1 + score2 + score3

    # need to transform "raw" value into scale between 1 and 10
    minv = -4.0
    maxv = 2.5
    sascore = 11.0 - (sascore - minv + 1) / (maxv - minv) * 9.0
    # smooth the 10-end
    if sascore > 8.0:
        sascore = 8.0 + math.log(sascore - 8.0)
    if sascore > 10.0:
        sascore = 10.0
    elif sascore < 1.0:
        sascore = 1.0

    return sascore
