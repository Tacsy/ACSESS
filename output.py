#!/usr/bin/env python
#-*- coding: utf-8 -*-

import time
import sys
import mprms
import numpy as np

startTime = time.time()
debug = False

##########################
# Initialize output files
##########################

statcolumns = ('gen', 'diversity', 'nPool', 'nLib', 'nCand', 'nFilt', 'nDups', 'nExcp',
               'nUnFit', 'nAdd', 'nAddFail', 'nAddArRing', 'nAddArRingFail',
               'nAtomType', 'nAtomTypeFail', 'nBreak', 'nBreakFail', 'nFlip',
               'nFlipFail', 'nNewRing', 'nNewRingFail', 'nRemove',
               'nRemoveFail', 'nNoMutation')
fitnesscolumns = ('gen', 'NumIn', 'NumOut', 'nSwap', 'AvgFVal', 'MinFVal',
                  'MaxFVal', 'CutOff')


def Init():
    global statsFile, filterFile, fitnessFile
    if mprms.restart:
        statsFile = open('stats.dat', 'a')
        filterFile = open('filters.dat', 'a')
    else:
        statsFile = open('stats.dat', 'w')
        filterFile = open('filters.dat', 'w')
        statsFile.write(" ".join(statcolumns[:7]))
        statsFile.write(" \\\n    ")
        statsFile.write(" ".join(statcolumns[7:]))
        statsFile.write("\n")

    # objective output
    if mprms.optimize:
        if mprms.restart:
            fitnessFile = open('Fitness.dat', 'a')
        else:
            fitnessFile = open('Fitness.dat', 'w')
            fitnessFile.write(" ".join(fitnesscolumns))
            fitnessFile.write("\n")
    return


##########################
# Statistics function
##########################

statsHead = "\n\tSTATISTICS:"
from collections import defaultdict
_default = lambda: 0
_default.__name__ = 'lambda:0'
stats = defaultdict(_default)
obstats=defaultdict(_default)


def PrintStat(nColumn=4, flush=True):
    global stats
    if len(stats) == 0:
        return
    print statsHead
    keys = PrintDict(stats, nColumn, sort='key')
    for i, key in enumerate(statcolumns):
        try:
            statsFile.write(" {:4d} ".format(stats[key]))
        except ValueError:
            statsFile.write(" {:14.7f} ".format(stats[key]))
        if i == 6: statsFile.write("|")
    statsFile.write("\n")
    if flush:
        sys.stdout.flush()
        stats.clear()
    return


def PrintObjectiveStat(nColumn=4, flush=True):
    global obstats
    for key in fitnesscolumns[:4]:
        fitnessFile.write(" {:4d} ".format(obstats[key]))
    for function in (np.average, np.min, np.max):
        fitnessFile.write(" {:14.7f} ".format(function(obstats['fvals'])))
    fitnessFile.write("\n")
    if flush:
        fitnessFile.flush()


##########################
# Timer function
##########################

timingHead = "\n\tTIMINGS:"
totalTimingHead = """\n
--------------------------------------------------------
                        TOTAL TIMINGS:              
--------------------------------------------------------"""

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


from functools import wraps


def logtime():
    def decorate(f):
        @wraps(f)
        def wrapped(*args, **kwargs):
            name = f.__name__.replace('Drive', ' ')
            StartTimer(name)
            r = f(*args, **kwargs)
            EndTimer(name)
            return r

        return wrapped

    return decorate


def PrintTimings(nColumn=4, flush=True):
    global timings, totalTimes
    keys = timings.keys()
    if len(keys) == 0:
        return
    values = timings.values()
    print timingHead
    PrintDict(timings, nColumn, truncate=True, sort='val')
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


def PrintTotalTimings(nColumn=4):
    global totalTimes
    print totalTimingHead
    if len(totalTimes) == 0:
        return
    PrintDict(totalTimes, truncate=True, sort='val')


##########################
# Printing function
##########################


def PrintDict(mydict, nColumn=4, truncate=False, sort=False):
    if sort != False:
        if sort == 'val':
            idx = 1
            reverse = True
        elif sort == 'key':
            idx = 0
            reverse = False
        keys, values = zip(
            *sorted(mydict.items(), key=lambda x: x[idx], reverse=reverse))
        keys = list(keys)
        values = list(values)
    else:
        keys = mydict.keys()
        values = mydict.values()

    if truncate:
        values = ["%.2f" % v for v in values]

    #complation of keys/values in the fashion that fits the format
    while len(keys) % nColumn != 0:
        keys.append(' ')
        values.append(' ')

    #formating key and value in row
    keyformat = ''
    valformat = ''
    for i in xrange(nColumn):
        keyformat += '+--{:<15}--'
        valformat += '|  {:<15}  '
    keyformat += '+'
    valformat += '|'

    for i in xrange(0, len(keys), nColumn):
        print keyformat.format(*keys[i:i + nColumn])
        print valformat.format(*values[i:i + nColumn])

    sys.stdout.flush()
    return keys
