#!/usr/bin/env python
#-*- coding: utf-8 -*-


class Counter(dict):
    #creat a counter of list in dictionary format

    def __init__(self, inlist, transformFunc=None, trunc=2):
        dict.__init__(self)
        self.transformFunc = transformFunc
        self.trunc = trunc
        self.AddList(inlist)

    def AddList(self, inlist):
        for item in inlist:
            if self.transformFunc is not None:
                self.Add(self.transformFunc(item))
            else:
                self.Add(item)

    def Add(self, item):
        if type(item) is float:
            item = round(item, self.trunc)

        if self.has_key(item):
            self[item] += 1
        else:
            self[item] = 1

    def Ranking(self, n=None):
        if n is None:
            n = len(self)
        items = self.items()
        items.sort(key=lambda x: -x[1])

        total = 0
        for i in items:
            total += i[1]

        maxKey = max(len(str(key)) for key in self.iterkeys())
        maxVal = max(len(str(val)) for val in self.itervalues())

        outformat = "{0:<" + str(maxKey) + "} {1:>" + str(maxVal) + "}    {2}"
        for key, value in items[0:n]:
            print outformat.format(key, value,
                                   ("%.2f" % (value * 100.0 / total)) + "%")

    def Combine(self, newDict):
        #combine two Counter objects together
        newKeys = newDict.keys()
        for key in newKeys:
            if type(newDict[key]) is not int:
                raise TypeError('passed object is not a counter')
        for key in newKeys:
            if not self.has_key(key):
                self[key] = 0
            self[key] += newDict[key]

    def Average(self):
        total = 0.0
        num = 0
        for key in self.keys():
            num += self[key]
            total += self[key] * key
        return total / (1.0 * num)

    def StdDev(self):
        avg = self.Average()
        num = 0
        stdDev = 0.0
        for key in self.keys():
            num += self[key]
            stdDev += self[key] * ((key - avg)**2)
        return stdDev / (1.0 * num)

    def Median(self):
        num = sum(self.values())
        nCount = 0
        for key in self.keys().sort():
            nCount += self[key]
            if nCount > num / 2:
                return key

    def Mode(self):
        return max(self.iteritems(), key=lambda x: x[1])[0]
