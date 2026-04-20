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
    db = sqlite3.connect(dbfile, detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES)
    cur = db.cursor()

    # Convert numpy data types to sqlite3 data types
    d2d = {"TEXT":(str,100), "INTEGER":int, "REAL":float}

    # Execute the select command
    cur.execute(sql)
    data = cur.fetchall()

    # No results
    if len(data)==0:
        return np.array([])

    # Get table column names and data types
    cur.execute("select sql from sqlite_master where tbl_name = '"+table+"'")
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

def getlc(dbfile,table='meas',cols='*',rar=None,decr=None,where=None,verbose=False):
    """ Get rows from the database """
    t0 = time.time()
    sqlite3.register_adapter(np.int8, int)
    sqlite3.register_adapter(np.int16, int)
    sqlite3.register_adapter(np.int32, int)
    sqlite3.register_adapter(np.int64, int)
    sqlite3.register_adapter(np.float16, float)
    sqlite3.register_adapter(np.float32, float)
    sqlite3.register_adapter(np.float64, float)
    db = sqlite3.connect(dbfile, detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES)
    cur = db.cursor()

    # Convert numpy data types to sqlite3 data types
    d2d = {"TEXT":(str,100), "INTEGER":int, "REAL":float}

    # Start the SELECT statement
    cmd = 'SELECT '+cols+' FROM '+table
    # RA constraints
    if rar is not None:
        if cmd.find('WHERE') == -1:
            cmd += ' WHERE '
        else:
            cmd += ' AND '
        cmd += 'ra>='+str(rar[0])+' AND ra<'+str(rar[1])
    # DEC constraints
    if decr is not None:
        if cmd.find('WHERE') == -1:
            cmd += ' WHERE '
        else:
            cmd += ' AND '
        cmd += 'dec>='+str(decr[0])+' AND dec<'+str(decr[1])
    # WHERE
    if where is not None:
        if cmd.lower().find('where')>-1:
            cmd += ' AND '+where
        else:
            cmd += ' WHERE '+where

    # Execute the select command
    #print('CMD = '+cmd)
    cur.execute(cmd)
    data = cur.fetchall()

    # No results
    if len(data)==0:
        return np.array([])

    # Get table column names and data types
    cur.execute("select sql from sqlite_master where tbl_name = '"+table+"'")
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
