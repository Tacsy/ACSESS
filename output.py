#!/usr/bin/env python
#-*- coding: utf-8 -*-

import time
import sys

startTime = time.time()

##########################
# Statistics function
##########################

statistics = {}
def PrintStat(nColumn, flush=True):
    global statistics
    if len(statistics) == 0:
        return
    PrintDict(statistics, nColumn)
    if flush:
        statistics = {}
    sys.stdout.flush()

##########################
# Timer function
##########################

timingHead = "\n    TIMINGS:\n"
totalTimingHead = """\n
--------------------------------------------------------
                        TOTAL TIMINGS                   
--------------------------------------------------------\n"""

timings = {}
timeRunning = {}
totalTimes = {}

def StartTimer(key):
    global timeRunning
    if debug:
        print "Starting timer: " + key
    if timeRunning.has_key(key):
        print "WARNING: " + key + "already running!"
    else:
        timeRunning[key] = time.time()

def EndTimer(key):
    global timeRunning, timings 
    if debug:
        print "Ending timer: " + key
    #get runtime 
    runTime = time.time() - timeRunning[key]
    timeRunning.pop(key)
    if timings.has_key(key):
        timings[key] += runTime
    else:
        timings[key] = runTime

def PrintTimings(nColumn, flush = True):
    global timings, totalTimes
    keys = timings.keys()
    if len(keys) == 0:
        return
    values = timings.values()
    print timingHead
    PrintDict(timings, nColumn, truncate = True, sort = 'val')
    if flush:
        for key, value in timings.iteritems():
            if totalTimes.has_key(key):
                totalTimes[key] += value
            else:
                totalTimes[key] = value
        timings = {}
    
    '''
    debug part should be considered
    '''

def PrintTotalTimings(nColumn):
    global totalTimes
    print totalTimingHead
    if len(totalTimes) == 0:
        return
    PrintDict(totalTimes, truncate = True, sort = 'val')


##########################
# Printing function
##########################

def PrintDcit(mydict, nColumn, truncate = False, sort = False):
    if sort != False:
        if sort == 'val':
            idx = 1
        elif sort  == 'key':
            idx = 0
        keys, values = zip(*sorted(mydict.items(), key = lambda x: -x[idx]))
        keys = list(keys)
        values = list(values)
    else:
        keys = mydict.keys()
        values = mydict.values()

    if truncate:
        values = ["%.2f"%v for v in values]
    
    #complation of keys/values in the fashion that fits the format
    while len(keys)%nColumn != 0:
        keys.append(' ')
        values.append(' ')

    #formating key and value in row
    keyformat = ''
    valformat = ''
    for i in xrange(nColumn):
        keyformat += '+--{:<15}--'
        valformat += '|  {:<15}}  '
    keyformat += '+'
    valformat += '|'

    for i in xrange(0, len(keys), nColumn):
        print keyformat.format(*keys[i:i+nColumn])
        print valformat.format(*values[i:i+nColumn])

    sys.stdout.flush()


