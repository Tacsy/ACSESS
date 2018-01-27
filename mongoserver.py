#!/usr/bin/env python
#-*- coding: utf-8 -*-
import sys
import os

'''
this is the script that copies from Chetan's PO-ACSESS version that calls
MongoDB to store some value
'''

hostfile = ''

def MongoServerInit():
    # establish connection to mongo server, save ACSESS database connection.
    global mydb, pymongo

    if not os.path.isfile(hostfile):
        raise IOError('Mongo server host file ' + hostfile + 'does not exist!')

    hostname = open(hostfile, 'r').readline().strip()

    import pymongo
    client = pymongo.MongoClient(hostname)
    mydb = client['ACSESS']

def UpdateDB(colname, values):
    # update database with all stored SMILES/value pairs

    mongoDocs = [{'_id':k, 'val':v} for k, v in values.iteritems()]

    col = mydb[colname]
    try:
        col.insert(mongoDocs, continue_on_error=True)
    except pymongo.errors.DuplicateKeyError as dke:
        print "WARNING: duplicate keys detected"
        print dke.args[0]

def LookupDB(colname, smiles):
    # return list with stored values for all SMILES in 'smiles'
    col = mydb[colname]
    query = col.find({"_id":{"$in":smiles}})

    memos = {q['_id']:q['val'] for q in query}
    results = [memos.get(s, None) for s in smiles]

    return results


