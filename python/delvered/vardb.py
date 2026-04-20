#!/usr/bin/env python                                                                                                                                                   
import os
import sys
import numpy as np
import glob
import warnings
from astropy.io import fits
from astropy.utils.exceptions import AstropyWarning
from astropy.table import Table, vstack, Column
from astropy.time import Time
import healpy as hp
from dlnpyutils import utils as dln, coords, db
import subprocess
import time
#from argparse import ArgumentParser
import socket
#from dustmaps.sfd import SFDQuery
from astropy.coordinates import SkyCoord
import sqlite3
import traceback

hostname = socket.gethostname()
host = hostname.split('.')[0]
if host=='thing' or host=='hulk' or host=='gp09':
    DBFILE = '/net/dl2/dnidever/delve/bricks/variables/delvemc_variables.db'
elif host=='tempest':
    DBFILE = '/home/group/davidnidever/nsc_variables/delvemc/delvemc_variables.db'
else:
    print('no dbfile')

def query(sql):
    """ Get rows from the database """
    t0 = time.time()
    sqlite3.register_adapter(np.int8, int)
    sqlite3.register_adapter(np.int16, int)
    sqlite3.register_adapter(np.int32, int)
    sqlite3.register_adapter(np.int64, int)
    sqlite3.register_adapter(np.float16, float)
    sqlite3.register_adapter(np.float32, float)
    sqlite3.register_adapter(np.float64, float)
    db = sqlite3.connect(DBFILE, detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES)
    cur = db.cursor()

    # Convert numpy data types to sqlite3 data types
    d2d = {"TEXT":(str,100), "INTEGER":int, "REAL":float}

    # Execute the select command
    cur.execute(sql)
    data = cur.fetchall()
    cols = [desc[0] for desc in cur.description]

    # No results
    if len(data)==0:
        return np.array([])

    # Get table column names and data types
    #cur.execute("select sql from sqlite_master where tbl_name = '"+table+"'")
    #dum = cur.fetchall()
    #db.close()
    #head = dum[0][0]
    ## 'CREATE TABLE exposure(expnum TEXT, nchips INTEGER, filter TEXT, exptime REAL, utdate TEXT, uttime TEXT, airmass REAL, wcstype TEXT)'
    #lo = head.find('(')
    #hi = head.find(')')
    #head = head[lo+1:hi]
    #cols = head.split(',')

    # infer types from first non-None value
    types = []
    for col_idx in range(len(cols)):
        col_values = [row[col_idx] for row in data if row[col_idx] is not None]
        inferred = type(col_values[0]) if col_values else None
        types.append(inferred)

    dt = []
    for i in range(len(cols)):
        dt.append( (cols[i], types[i]) )
    dtype = np.dtype(dt)

    # Convert to nump structured array
    tab = np.zeros(len(data),dtype=dtype)
    tab[...] = data
    del(data)

    return tab    


def getobject(objid,verbose=False):
    """ Get rows from the database """
    t0 = time.time()
    sqlite3.register_adapter(np.int8, int)
    sqlite3.register_adapter(np.int16, int)
    sqlite3.register_adapter(np.int32, int)
    sqlite3.register_adapter(np.int64, int)
    sqlite3.register_adapter(np.float16, float)
    sqlite3.register_adapter(np.float32, float)
    sqlite3.register_adapter(np.float64, float)
    db = sqlite3.connect(DBFILE, detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES)
    cur = db.cursor()

    # Convert numpy data types to sqlite3 data types
    d2d = {"TEXT":(str,100), "INTEGER":int, "REAL":float}

    sql = "SELECT * from object WHERE objid='"+str(objid)+"'"

    # Execute the select command
    #print(sql)
    cur.execute(sql)
    data = cur.fetchall()

    # No results
    if len(data)==0:
        return np.array([])

    # Get table column names and data types
    cur.execute("select sql from sqlite_master where tbl_name = 'object'")
    dum = cur.fetchall()
    db.close()
    head = dum[0][0]
    # 'CREATE TABLE exposure(expnum TEXT, nchips INTEGER, filter TEXT, exptime REAL, utdate TEXT, uttime TEXT, airmass REAL, wcstype TEXT)'
    lo = head.find('(')
    hi = head.find(')')
    head = head[lo+1:hi]
    cols = head.split(',')
    cols = dln.strip(cols)
    dt = []
    for c in cols:
        pair = c.split(' ')
        dt.append( (pair[0], d2d[pair[1]]) )
    dtype = np.dtype(dt)

    # Convert to nump structured array
    tab = np.zeros(len(data),dtype=dtype)
    tab[...] = data
    del(data)

    if verbose: print('got data in '+str(time.time()-t0)+' sec.')

    return tab


def getlc(objid,verbose=False):
    """ Get rows from the database """
    t0 = time.time()
    sqlite3.register_adapter(np.int8, int)
    sqlite3.register_adapter(np.int16, int)
    sqlite3.register_adapter(np.int32, int)
    sqlite3.register_adapter(np.int64, int)
    sqlite3.register_adapter(np.float16, float)
    sqlite3.register_adapter(np.float32, float)
    sqlite3.register_adapter(np.float64, float)
    db = sqlite3.connect(DBFILE, detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES)
    cur = db.cursor()

    # Convert numpy data types to sqlite3 data types
    d2d = {"TEXT":(str,100), "INTEGER":int, "REAL":float}

    sql = "SELECT * from meas WHERE objid='"+str(objid)+"'"

    # Execute the select command
    #print(sql)
    cur.execute(sql)
    data = cur.fetchall()

    # No results
    if len(data)==0:
        return np.array([])

    # Get table column names and data types
    cur.execute("select sql from sqlite_master where tbl_name = 'meas'")
    dum = cur.fetchall()
    db.close()
    head = dum[0][0]
    # 'CREATE TABLE exposure(expnum TEXT, nchips INTEGER, filter TEXT, exptime REAL, utdate TEXT, uttime TEXT, airmass REAL, wcstype TEXT)'
    lo = head.find('(')
    hi = head.find(')')
    head = head[lo+1:hi]
    cols = head.split(',')
    cols = dln.strip(cols)
    dt = []
    for c in cols:
        pair = c.split(' ')
        dt.append( (pair[0], d2d[pair[1]]) )
    dtype = np.dtype(dt)

    # Convert to nump structured array
    tab = np.zeros(len(data),dtype=dtype)
    tab[...] = data
    del(data)

    if verbose: print('got data in '+str(time.time()-t0)+' sec.')

    return tab
