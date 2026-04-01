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

def writecat2db(cat,dbfile,table='meas'):
    """ Write a catalog to the database """
    ncat = dln.size(cat)
    sqlite3.register_adapter(np.int8, int)
    sqlite3.register_adapter(np.int16, int)
    sqlite3.register_adapter(np.int32, int)
    sqlite3.register_adapter(np.int64, int)
    sqlite3.register_adapter(np.float16, float)
    sqlite3.register_adapter(np.float32, float)
    sqlite3.register_adapter(np.float64, float)
    db = sqlite3.connect(dbfile, detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES)
    c = db.cursor()

    # Convert numpy data types to sqlite3 data types
    d2d = {"S":"TEXT", "s":"TEXT", "U":"TEXT", "u":"TEXT", "i":"INTEGER", "f":"REAL"}

    # Get the column names
    cnames = cat.dtype.names
    cdict = dict(cat.dtype.fields)
    # Create the table
    #   the primary key ROWID is automatically generated
    if len(c.execute('SELECT name from sqlite_master where type= "table" and name="'+table+'"').fetchall()) < 1:
        columns = cnames[0].lower()+' '+d2d[cdict[cnames[0]][0].kind]
        for n in cnames[1:]: columns+=', '+n.lower()+' '+d2d[cdict[n][0].kind]
        c.execute('CREATE TABLE '+table+'('+columns+')')
    # Insert statement
    columns = []
    for n in cnames: columns.append(n.lower())
    qmarks = np.repeat('?',dln.size(cnames))
    data = list(cat)
    # Replace nan with 'nan'  
    data = [
        tuple('NAN' if isinstance(i, np.floating) and np.isnan(i) else i for i in t)
        for t in list(cat)
    ]
    # Replace inf with 'inf'  
    data = [
        tuple(str(i).upper() if isinstance(i, np.floating) and np.isinf(i) else i for i in t)
        for t in list(data)
    ]
    c.executemany('INSERT INTO '+table+'('+','.join(columns)+') VALUES('+','.join(qmarks)+')', data)
    #c.executemany('INSERT INTO '+table+'('+','.join(columns)+') VALUES('+','.join(qmarks)+')', list(cat))
    db.commit()
    db.close()

def createindexdb(dbfile,col='measid',table='meas',unique=True,verbose=False):
    """ Index a column in the database """
    t0 = time.time()
    db = sqlite3.connect(dbfile, detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES)
    c = db.cursor()
    index_name = 'idx_'+col+'_'+table
    # Check if the index exists first
    c.execute('select name from sqlite_master')
    d = c.fetchall()
    for nn in d:
        if nn[0]==index_name:
            print(index_name+' already exists')
            return
    # Create the index
    if verbose: print('Indexing '+col)
    if unique:
        c.execute('CREATE UNIQUE INDEX '+index_name+' ON '+table+'('+col+')')
    else:
        c.execute('CREATE INDEX '+index_name+' ON '+table+'('+col+')')
    data = c.fetchall()
    db.close()
    if verbose: print('indexing done after '+str(time.time()-t0)+' sec')

def getdatadb(dbfile,table='meas',cols='*',rar=None,decr=None,where=None,verbose=False):
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
    cat = np.zeros(len(data),dtype=dtype)
    cat[...] = data
    del(data)

    if verbose: print('got data in '+str(time.time()-t0)+' sec.')

    return cat

def createsumtable(dbfile=None,delvedir='/net/dl2/dnidever/delve/'):
    """ Create the DELVE-MC exposures and chips database"""

    expdir = delvedir+'exposures/'
    if dbfile is None:
        dbfile = delvedir+'bricks/db/delvered_summary.db'

    # Loop over list of nightly summary files

    # Get all of the nights
    dirs = glob.glob(expdir+'20??????')
    dirs = np.array(dirs)
    ndirs = dln.size(dirs)
    nights = np.zeros(dirs.size,(str,10))
    for i,d in enumerate(dirs): nights[i]=os.path.basename(d)
    si = np.argsort(nights)
    nights = nights[si]
    dirs = dirs[si]
    # Nightly summary files
    allnightsumfiles = dln.strjoin(dirs,'/',nights)
    allnightsumfiles = dln.strjoin(allnightsumfiles,'_summary.fits')
    gdnightsumfiles,nnightsumfiles = dln.where(dln.exists(allnightsumfiles))
    print(str(nnightsumfiles)+' nightly summary files found')
    nightsumfiles = allnightsumfiles[gdnightsumfiles]
    # Load in the data
    expcount = 0
    chcount = 0
    lexpstr = []  # lists
    lchstr = []
    for i in range(nnightsumfiles):
        print('Loading '+nightsumfiles[i])
        #st = os.stat(nightsumfiles[i])
        try:
            expstr1 = fits.getdata(nightsumfiles[i],1,memmap=False)
            nexpstr1 = dln.size(expstr1)
            lexpstr.append(expstr1)
            chstr1 = fits.getdata(nightsumfiles[i],2,memmap=False)
            nchstr1 = dln.size(chstr1)
            lchstr.append(chstr1)
        except:
            print('Problem loading '+nightsumfiles[i])
        
       # # Load the database
       # writecat2db(expstr1,dbfile,'exposure')
       # writecat2db(chstr1,dbfile,'chip')

    # Concatenate them
    expstr = dln.concatenate(lexpstr)
    chstr = dln.concatenate(lchstr)

    # Load the database
    writecat2db(expstr,dbfile,'exposure')
    writecat2db(chstr,dbfile,'chip')


    # Create the indices
    print('Indexing')
    createindexdb(dbfile,'ra','exposure',unique=False)
    createindexdb(dbfile,'dec','exposure',unique=False)
    db.analyzetable(dbfile,'exposure')
    createindexdb(dbfile,'ra','chip',unique=False)
    createindexdb(dbfile,'dec','chip',unique=False)
    db.analyzetable(dbfile,'chip')


def createbricksumtable(dbfile=None,delvedir='/net/dl2/dnidever/delve/'):
    """ Create the DELVE-MC brick metadata summary database """

    basedir = delvedir+'bricks/'
    if dbfile is None:
        dbfile = delvedir+'bricks/db/delvered_brick_summary.db'


    bricks = Table.read('/home/dnidever/projects/delvered/data/delvemc_bricks_0.25deg.fits.gz')
    for c in bricks.colnames: bricks[c].name = c.lower()
    si = np.argsort(np.random.rand(len(bricks)))
    bricks = bricks[si]
    # Load in the data

    
    dt = [('field', 'U10'), ('file', 'U200'), ('expnum', 'U8'), ('chip', '>i8'),
          ('base', 'U50'), ('filter', 'U3'), ('exptime', '>f8'), ('utdate', 'U15'),
          ('uttime', 'U20'), ('airmass', '>f8'), ('gain', '>f8'), ('rdnoise', '>f8'),
          ('nx', '>i8'), ('ny', '>i8'), ('wcstype', 'U10'), ('pixscale', '>f8'),
          ('ra', '>f8'), ('dec', '>f8'), ('wcsrms', '>f8'), ('fwhm', '>f8'),
          ('skymode', '>f8'), ('skysig', '>f8'), ('dao_nsources', '>i8'),
          ('dao_depth', '>f8'), ('dao_npsfstars', '>i8'), ('dao_psftype', 'U10'),
          ('dao_psfboxsize', '>i8'), ('dao_psfvarorder', '>i8'), ('dao_psfchi', '>f8'),
          ('alf_nsources', '>i8'), ('alf_depth', '>f8'), ('calib_depth', '>f8'),
          ('calib_color', 'U10'), ('calib_zpterm', '>f8'), ('calib_zptermsig', '>f8'),
          ('calib_amterm', '>f8'), ('calib_amtermsig', '>f8'),
          ('calib_colorterm', '>f8'), ('calib_colortermsig', '>f8'),
          ('calib_magname', 'U10'), ('apcor', '>f8'), ('ebv', '>f8'),
          ('fieldname', 'U50'), ('depth', '>f4'), ('eta', '>f4'),
          ('background', '>f4'), ('tau', '>f4'), ('teff', '>f4'),
          ('fracoverlap', '>f4'), ('mnx', '>f4'), ('mny', '>f4'),
          ('brickname', 'U10'), ('vx1', '>f4'), ('vx2', '>f4'), ('vx3', '>f4'),
          ('vx4', '>f4'), ('vy1', '>f4'), ('vy2', '>f4'), ('vy3', '>f4'), ('vy4', '>f4')]
    data = []
    for i in range(len(bricks)):
        name = bricks['brickname'][i]
        print(i,name)
        metafile = basedir+name[:4]+'/'+name+'/'+name+'_meta.fits'
        if os.path.exists(metafile)==False:
            continue
        try:
            meta = Table.read(metafile)
            for c in meta.colnames:meta[c].name=c.lower()
            meta['brickname'] = name
            if 'vx' not in meta.colnames:
                # should make the vx/vy arrays
                meta['vx1'] = 0.0
                meta['vx2'] = 0.0
                meta['vx3'] = 0.0
                meta['vx4'] = 0.0
                meta['vy1'] = 0.0
                meta['vy2'] = 0.0
                meta['vy3'] = 0.0
                meta['vy4'] = 0.0
            else:
                meta['vx1'] = meta['vx'][:,0]
                meta['vx2'] = meta['vx'][:,1]
                meta['vx3'] = meta['vx'][:,2]
                meta['vx4'] = meta['vx'][:,3]
                meta['vy1'] = meta['vy'][:,0]
                meta['vy2'] = meta['vy'][:,1]
                meta['vy3'] = meta['vy'][:,2]
                meta['vy4'] = meta['vy'][:,3]
                del meta[['vx','vy']]

            newmeta = np.zeros(len(meta),dtype=np.dtype(dt))
            for c in newmeta.dtype.names:
                if newmeta[c].dtype.kind=='f':
                    newmeta[c] = np.nan            
            for c in meta.colnames: newmeta[c] = meta[c]
            data.append(newmeta)
        except KeyboardInterrupt:
            raise
        except:
            traceback.print_exc()
            print('Problem loading '+metafile)

    # Concatenate them
    #data = vstack(data)
    data = np.concatenate(data)

    # Write the database
    print('Writing the database')
    writecat2db(data,dbfile,'chip')

    # Create the indices
    print('Indexing')
    createindexdb(dbfile,'brickname','chip',unique=False)
    createindexdb(dbfile,'ra','chip',unique=False)
    createindexdb(dbfile,'dec','chip',unique=False)
    createindexdb(dbfile,'file','chip',unique=False)
    db.analyzetable(dbfile,'chip')


    import pdb; pdb.set_trace()
