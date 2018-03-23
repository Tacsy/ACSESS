#!/usr/bin/env python
#-*- coding: utf-8 -*-
'''
 This forms an middle layer between the main algorithm steps and the lowlevel functions
'''
import mutate
import filters


def DriveMutations(lib):
    newlib=lib
    return newlib

def DriveFilters(lib):

    # filter by setting the failed attribute True
    for mol in lib:
        changed, failed = filters.FixAndFilter(mol)
        #if changed: cn.Finalize(mol)
        mol.SetBoolProp('failed', failed)

    # effective filter step. 
    newlib=filter(lambda mol:not mol.GetBoolProp('failed'), lib)
    return newlib

def ExtendPool(pool, lib, newlib):
    return lib
