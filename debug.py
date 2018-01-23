#!/usr/bin/env python
#-*- coding: utf-8 -*-
import sys

def GetRefCounts():
    
    refDict = {}
    
    #collect all classes
    for module in sys.modules.values():
        for sym in dir(module):
            obj = getattr(module,sym)
            if type(obj) is type:
                refDict[obj] = sys.getrefcount(obj)

    #sort and reverse by refcount
    pairs = map(lambda x:(x[1],x[0]), refDict.items())
    pairs.sort()
    pairs.reverse()

    return pairs

def PrintRefCounts(num):
    #print the top N refcounts classes
    print '------ MOST ' + str(num) + ' REFERENCED CLASSES ------'
    
    for num, cls in GetRefCounts()[:num]:
        print '%10d %s' %(num, cls.__name__)

