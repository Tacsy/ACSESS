#!/usr/bin/env python
#-*- coding: utf-8 -*-
import sys
import time

import numpy as np
from mpi4py import MPI

import output
mpi = False
'''
this module include various previous modules that coorespnds to parallelization
in order to make the code in a cleaner way, including:
    Parallel.py
    ParallelServer.py
'''
############################################################
#       Functions defined in Parallel.py
############################################################

############################################################
#       Class defined for Parallel Computing (ParallelServer.py)
############################################################

# A class for SIMD (single instruction multiple data) parallelization.
# Any "instruction" (i.e. python function) that will be called must
# be registered with the object on all nodes by calling
# MPITaskInstance.RegisterFunction.
#
# The function to be called is broadcast from the master using SetFunction.
# The master initiates a parallel calculation by calling MPIRun with
# the data array to be processed.
# The MPIRun function will return an array of the results for each datum.
QUIT = 42
SHUTDOWN = 'SHUTDOWN'
signaldone = [np.array([1], dtype='i'), 'i']


class MPITask():
    def __init__(self, myComm=MPI.COMM_WORLD):
        self.comm = myComm
        self.size = self.comm.Get_size()
        self.rank = self.comm.Get_rank()
        self.POLLTIME = 4
        self.verbose = False
        self.DistFunc = None
        self.RunMaster = 5
        self.RegisteredFunctions = {SHUTDOWN: None}
        self.IdleTime = 0
        self.BCastTime = 0
        if self.rank == 0:
            self.IsMaster = True
            self.stdout = sys.stdout
        else:
            self.IsMaster = False
            sys.stdout = open('node%i.out' % self.rank, 'a')
            self.stdout = sys.stdout
            sys.stderr = self.stdout

    #send shutdown message on deletion
    def Shutdown(self):
        if self.IsMaster:
            darpa = self.comm.bcast(SHUTDOWN, root=0)
            print 'Master node idle time: ', self.IdleTime
        else:
            sys.exit()

    #both master and children must call this simultaneously to register
    #the functions that will be called with MPI
    def RegisterFunction(self, funclist):
        for func in funclist:
            if not self.RegisteredFunctions.has_key(func.func_name):
                self.RegisteredFunctions[func.func_name] = func
                if self.verbose:
                    print 'MPI Task Registered: ' + func.func_name
            else:
                raise NameError('Function with name "' + func.func_name +
                                '" is already registered.')

    #set the function to distribute
    def SetFunction(self, *args):
        fname = ''
        if self.IsMaster:
            func = args[0]
            fname = func.func_name
        elif self.verbose:
            print 'Waiting for function broadcast ...'
            sys.stdout.flush()
        starttime = time.time()
        fname = self.comm.bcast(fname, root=0)
        self.BCastTime += time.time() - starttime
        if fname == SHUTDOWN:
            print 'Received shutdown message.'
            print 'Idle time waiting for instructions: ', self.BCastTime
            print 'Idle time during calculations: ', self.IdleTime
            sys.exit()

        self.DistFunc = self.RegisteredFunctions[fname]
        if self.verbose:
            print "Set to function: '" + fname + "'"
            sys.stdout.flush()

    #mpi run
    def RunMPI(self, distData):
        starttime = time.time()
        nData = self.comm.bcast([len(distData)], root=0)  #synchronize w/root
        self.IdleTime += time.time() - starttime
        if self.rank == 0:
            return self.ScatterDriver(distData)

        stats = MPI.Status()

        if self.verbose():
            print 'Starting calculation (' + str(nData) + ' total inputs)'

        while True:
            output.StartTimer('WAITING')
            if self.verbose:
                print 'Waiting for data ...'
            sys.stdout.flush()
            SLICE = self.comm.recv(source=0, tag=MPI.ANY_TAG, status=stats)
            output.EndTimer('WAITING')

            if stats.Get_tag() == QUIT:
                break
            if self.verbose:
                print 'Working on ' + str(len(SLICE)) + ' data slices'
                sys.stdout.flush()
            results = [self.DistFunc(item) for item in SLICE]

            output.StartTimer('WAITING')
            self.comm.Send(signaldone, dest=0)
            output.EndTimer('WAITING')

            self.comm.send(results, dest=0)

        self.stdout.flush()

    #main function for scatter calculation
    def ScatterDriver(self, distData):
        if self.verbose:
            print 'Scattering ' + str(len(distData)) + ' total inputs'
        assert self.rank == 0
        CurrentChunk = [
            -1
        ] * self.size  #CurrentChun[i] == -1 indicates node is idle
        request = [None] * self.size
        rflag = [np.array([0], dtype='i') for i in xrange(self.size)]

        chunk = max(int(len(distData) / (8.0 * self.size)), 1)
        Done = False
        length = len(distData)
        results = [None] * length
        lastitem = 0
        waitingfor = None
        while not Done:
            for i in xrange(1, self.size):
                #receive results from any workers that are done
                if request[i] is not None and request[i].Test():
                    resultbuff = self.comm.recv(source=i)
                    request[i] = None
                    results[CurrentChunk[i]:
                            CurrentChunk[i] + chunk] = resultbuff
                    sys.stdout.flush()
                    if self.verbose and lastitem < length and masterwork != lastitem:
                        print 'Master completed slices ', masterwork, 'to', lastitem
                        print 'RECV from', i, ': slices', CurrentChunk[
                            i], ' to ', CurrentChunk[i] + chunk - 1
                    CurrentChunk[i] = -1

                #instruct node to start on next segment
                if lastitem < length and request[i] is None:
                    self.comm.send(distData[lastitem:lastitem + chunk], dest=i)
                    CurrentChunk[i] = lastitem
                    lastitem += chunk
                    request[i] = self.comm.Irecv(rflag[i], source=i)
                    if self.verbose:
                        print 'SEND to ', i, ': slices', CurrentChunk[i], \
                              ' to ', CurrentChunk[i] + chunk - 1
                        masterwork = lastitem

            sys.stdout.flush()

            #see if everyone's finished
            if lastitem >= length:
                NotFinished = []
                for i, req in enumerate(request):
                    if req is not None:
                        NotFinished.append(i)
                NotFinished = tuple(NotFinished)

                if len(NotFinished) == 0:
                    if waitingfor is not None:
                        print ''
                    break
                else:
                    if waitingfor != NotFinished:
                        waittime = 0
                        if waitingfor is not None:
                            print ''
                        print 'Waiting for nodes', \
                              ', '.join(str(i) for i in NotFinished), \
                              ' to finish ...'
                        waitingfor = NotFinished
                        sleeptime = 1
                    else:
                        waittime += self.POLLTIME
                        print str(waittime) + 's'
                        sleeptime = self.POLLTIME
                    sys.stdout.flush()
                    time.sleep(sleeptime)
                    self.IdleTime += sleeptime

            #do computation while waiting for requests to be finished
            elif self.RunMaster > 0:
                numtodo = min(self.RunMaster, chunk / 2)
                if numtodo < 1:
                    numtodo = 1
                results[lastitem:lastitem+numtodo] = \
                    [self.DistFunc(item) for item in
                            DistData[lastitem:lastitem+numtodo]]
                lastitem += numtodo
            else:
                time.sleep(self.POLLTIME)
                self.IdleTime += self.POLLTIME

        #Once finished, send quit signal and return results
        if self.verbose:
            print '******MPI TASK COMPLETED******'
        for i in xrange(1, self.size):
            self.comm.isend(None, dest=i, tag=QUIT)

        return results
