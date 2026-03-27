import os
import numpy as np
from astropy.table import Table,vstack
from dlnpyutils import utils as dln
from glob import glob
import traceback

def getvars(brick,ndetthresh=20,save=False,overwrite=False):
    """ Get variable star data for a brick. """
    bdir = '/net/dl2/dnidever/delve/bricks/'+brick[:4]+'/'+brick
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
    varobj = jobj[gd]
    
    # Load measurements
    measfile = bdir+'/'+brick+'_joint_meas.fits.gz'
    if os.path.exists(measfile):
        meas = Table.read(measfile)
        meas['objid'] = meas['objid'].data.astype(str)
        index = dln.create_index(meas['objid'])
        _,ind1,ind2 = np.intersect1d(index['value'],varid,return_indices=True)
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
        measfiles = glob(bdir+'/'+brick+'_joint_meas_*.fits')
        measfiles.sort()
        meas = []
        for i in range(len(measfiles)):
            meas1 = Table.read(measfiles[i])
            expnum1 = os.path.basename(measfiles[i])[:-5].split('_')[-1]
            meas1['expnum'] = 8*' '
            meas1['expnum'][:] = expnum1
            for c in meas1.colnames:meas1[c].name=c.lower()
            meas1['objid'] = meas1['objid'].data.astype(str)
            _,ind1,ind2 = np.intersect1d(meas1['objid'],varid,return_indices=True)
            print('  {:5d}  {:<40s}  {:10s}  {:8d}  {:5d}'.format(i+1,os.path.basename(measfiles[i]),expnum1,len(meas1),len(ind1)))
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
    for i in range(1528,len(bricks)):
        print(i+1,bricks['brickname'][i])
        try:
            out = getvars(bricks['brickname'][i],save=True,overwrite=True)
        except KeyboardInterupt:
            raise
        except:
            traceback.print_exc()
