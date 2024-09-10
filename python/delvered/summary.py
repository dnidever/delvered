import os
import numpy as np
from astropy.table import Table,vstack
from astropy.io import fits
from dlnpyutils import utils as dln,coords
from astropy.coordinates import SkyCoord
import astropy.units as u
from gala.coordinates import MagellanicStreamNidever08

def bricks():
    """ Brick-level summary."""

    bricks = Table.read('/home/dnidever/projects/delvered/data/delvemc_bricks_0.25deg.fits.gz')
    for c in bricks.colnames: bricks[c].name = c.lower()

    # Get Magellanic Stream coordinates
    coo = SkyCoord(ra=bricks['ra']*u.degree, dec=bricks['dec']*u.degree,frame='icrs')
    mcoo = coo.transform_to(MagellanicStreamNidever08)
    bricks['mlon'] = mcoo.L.deg
    bricks['mlat'] = mcoo.B.deg

    bricks['nchips'] = -1
    bricks['nexp'] = -1
    bricks['ngexp'] = -1
    bricks['nrexp'] = -1
    bricks['niexp'] = -1
    bricks['nzexp'] = -1
    bricks['nyexp'] = -1
    bricks['exptime'] = 0.0
    bricks['gexptime'] = 0.0
    bricks['rexptime'] = 0.0
    bricks['iexptime'] = 0.0
    bricks['zexptime'] = 0.0
    bricks['yexptime'] = 0.0
    bricks['nobject'] = -1
    bricks['njmeas'] = -1
    bricks['njobject'] = -1
    basedir = '/net/dl2/dnidever/delve/bricks/'
    filters = ['g','r','i','z','Y']
    for i in range(len(bricks)):
        name = bricks['brickname'][i]
        print(i,name)
        metafile = basedir+name[:4]+'/'+name+'/'+name+'_meta.fits'
        #metafile = basedir+name[:4]+'/'+name+'/'+name+'_joint_meta.fits'
        if os.path.exists(metafile):
            meta = Table.read(metafile,1)
            bricks['nchips'][i] = len(meta)
            bricks['nexp'][i] = len(np.unique(meta['EXPNUM']))
            _,ui = np.unique(meta['EXPNUM'],return_index=True)
            bricks['exptime'][i] = np.sum(meta['EXPTIME'][ui])
            for f in filters:
                ind, = np.where(meta['FILTER']==f)
                if len(ind)>0:
                    fmeta = meta[ind]
                    _,ui = np.unique(fmeta['EXPNUM'],return_index=True)
                    bricks['n'+f.lower()+'exp'][i] = len(ui)
                    bricks[f.lower()+'exptime'][i] = np.sum(fmeta['EXPTIME'][ui])
        else:
            print(metafile,'not found')
        objfile = basedir+name[:4]+'/'+name+'/'+name+'_object.fits.gz'
        jmeasfile = basedir+name[:4]+'/'+name+'/'+name+'_joint_meas.fits.gz'
        jobjfile = basedir+name[:4]+'/'+name+'/'+name+'_joint_object.fits.gz'
        if os.path.exists(objfile):
            ohead = fits.getheader(objfile,1)
            bricks['nobject'][i] = ohead['naxis2']
        if os.path.exists(jmeasfile):
            mhead = fits.getheader(jmeasfile,1)
            bricks['njmeas'][i] = mhead['naxis2']
        if os.path.exists(jobjfile):
            ojhead = fits.getheader(jobjfile,1)
            bricks['njobject'][i] = ojhead['naxis2']

    return bricks
