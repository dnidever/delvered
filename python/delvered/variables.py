import os
import numpy as np
from astropy.table import Table,vstack
from dlnpyutils import utils as dln,db as dlndb
from glob import glob
from delvered import delvered_db as db
import traceback

def getvars(brick,ndetthresh=20,save=False,overwrite=False,loop=True):
    """ Get variable star data for a brick. """
    bdir = '/net/dl2/dnidever/delve/bricks/'+brick[:4]+'/'+brick
    metafile = bdir+'/'+brick+'_joint_meta.fits'
    if os.path.exists(metafile)==False:
        print(metafile,'NOT FOUND')
        return
    jobjfile = bdir+'/'+brick+'_joint_object.fits.gz'
    if os.path.exists(jobjfile)==False:
        print(jobjfile,'NOT FOUND')
        return
    jobj = Table.read(jobjfile)
    for c in jobj.colnames:jobj[c].name=c.lower()
    gd, = np.where((jobj['variable10sig']==1) & (jobj['ndet']>=ndetthresh))
    if len(gd)==0:
        print('  No variables for',brick,'pass the thresholds')
        return
    print('  ',len(gd),'variables')
    varid = jobj['objid'][gd].data.astype(str)
    varid = np.array([str(f).strip() for f in varid])  # strip extra spaces
    varobj = jobj[gd]
    
    # Load measurements
    measfile = bdir+'/'+brick+'_joint_meas.fits.gz'
    if os.path.exists(measfile) and loop==False:
        print('Using single meas file',measfile)
        meas = Table.read(measfile)
        for c in meas.colnames: meas[c].name = c.lower()
        meas['objid'] = meas['objid'].data.astype(str)
        index = dln.create_index(meas['objid'])
        objid = [str(f).strip() for f in index['value'].astype(str)]
        _,ind1,ind2 = np.intersect1d(objid,varid,return_indices=True)
        allind = []
        for i in range(len(ind1)):
            ind = index['index'][index['lo'][ind1[i]]:index['hi'][ind1[i]]+1]
            allind.append(ind)
        allind = np.concatenate(allind)
        allind = np.sort(allind)
        meas = meas[allind]

    # Split in multiple meas files
    else:
        #measfiles = glob(bdir+'/'+brick+'_joint_meas_????????.fits')
        #measfiles = glob(bdir+'/'+brick+'_joint_meas_*.fits')
        #measfiles.sort()
        meta = Table.read(metafile,1)
        for c in meta.colnames:meta[c].name=c.lower()
        measfiles = Table.read(metafile,2)
        for c in measfiles.colnames:measfiles[c].name=c.lower()
        meas = []
        for i in range(len(measfiles)):
            filename = str(measfiles['filename'][i]).strip()
            meas1 = Table.read(filename)
            for c in meas1.colnames:meas1[c].name=c.lower()
            expnum1 = os.path.basename(filename)[:-5].split('_')[-1]
            meas1['expnum'] = 8*' '
            meas1['expnum'][:] = expnum1
            objid = [str(f).strip() for f in meas1['objid'].data.astype(str)]
            _,ind1,ind2 = np.intersect1d(objid,varid,return_indices=True)
            print('  {:5d}  {:<40s}  {:10s}  {:8d}  {:5d}'.format(i+1,os.path.basename(filename),expnum1,len(meas1),len(ind1)))
            if len(ind1)>0:
                meas.append(meas1[ind1])
        meas = vstack(meas)
        
    print(' ',len(meas),'variable measurements')

    if save:
        print('  Saving')
        varobj.write(bdir+'/'+brick+'_variables_object.fits',overwrite=overwrite)
        meas.write(bdir+'/'+brick+'_variables_meas.fits',overwrite=overwrite)

    return varobj,meas

def allbricks():
    """ get all bricks """

    bricks = Table.read('/home/dnidever/projects/delvered/data/delvemc_bricks_0.25deg.fits.gz')
    for c in bricks.colnames: bricks[c].name = c.lower()
    #for i in range(len(bricks)):
    for i in range(4332,len(bricks)):
        print(i+1,bricks['brickname'][i])
        try:
            out = getvars(bricks['brickname'][i],save=True,overwrite=True)
        except KeyboardInterrupt:
            raise
        except:
            traceback.print_exc()

def load_database(dbfile=None):
    """" Load all of the brick variable star data into an sqlite3 database """

    if dbfile is None:
        dbfile = '/net/dl2/dnidever/delve/bricks/variables/delvemc_variables.db'
    if os.path.exists(dbfile):
        print(dbfile,'already exists')
        return

    bricks = Table.read('/home/dnidever/projects/delvered/data/delvemc_bricks_0.25deg.fits.gz')
    for c in bricks.colnames: bricks[c].name = c.lower()

    for i in range(len(bricks)):
        name = bricks['brickname'][i]
        bdir = '/net/dl2/dnidever/delve/bricks/'
        bdir = os.path.join(bdir,name[:4],name)
        measfile = bdir+'/'+name+'_variables_meas.fits'
        objfile = bdir+'/'+name+'_variables_object.fits'
        if os.path.exists(measfile) and os.path.exists(objfile):
            meas = Table.read(measfile)
            obj = Table.read(objfile)
            obj['brick'] = name
            print(i+1,name,len(obj),len(meas))
            db.writecat2db(meas,dbfile,'meas')
            db.writecat2db(obj,dbfile,'object')
        else:
            print(i+1,name,'variable files not found')

    # make indices
    print('Indexing')
    db.createindexdb(dbfile,'ra','object',unique=False)
    db.createindexdb(dbfile,'dec','object',unique=False)
    db.createindexdb(dbfile,'objid','object')
    dlndb.analyzetable(dbfile,'object')
    db.createindexdb(dbfile,'objid','meas',unique=False)
    dlndb.analyzetable(dbfile,'meas')

