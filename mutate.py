#!/usr/bin/env python
#-*- coding: utf-8 -*-

import random
from rdkit import Chem
import mprms









###############################################
# Initialization the module
###############################################

def MutateInit():
    elements = mprms.elements
    halogens = [9,17,35,53]
    global Elements, Halogens
    
    #do intersection and difference 
    Elements = set(elements) - set(halogens)
    Halogens = set(elements) & set(halogens)
    #change from set to list
    Elements = list(Elements)
    Halogens = list(Halogens)
    

   
###############################################################################
#                            Mutation Methods                                 #
###############################################################################


###############################################
#returns free valence for an atom
#bbviously, there's a problem if it's negative
#after implicit hydrogens are assigned with the valence model,
#this can be replaced by the implicit hydrogen count
###############################################



