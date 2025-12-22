import os
import numpy as np
from glob import glob
from astropy.table import Table,vstack
from astropy.io import fits
from dlnpyutils import utils as dln,coords,job_daemon as jd
from astropy.coordinates import SkyCoord
import astropy.units as u
from gala.coordinates import MagellanicStreamNidever08
from scipy.stats import binned_statistic,binned_statistic_2d

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
            try:
                ohead = fits.getheader(objfile,1)
                bricks['nobject'][i] = ohead['naxis2']
            except:
                print('problem loading',objfile)
        if os.path.exists(jmeasfile):
            try:
                mhead = fits.getheader(jmeasfile,1)
                bricks['njmeas'][i] = mhead['naxis2']
            except:
                print('problem loading',jmeasfile)
        if os.path.exists(jobjfile):
            try:
                ojhead = fits.getheader(jobjfile,1)
                bricks['njobject'][i] = ojhead['naxis2']
            except:
                print('problem loading',jobjfile)

    return bricks

def allepochs(clobber=False,nmulti=10):
    """ Run through all of the bricks and getting the epoch data"""
    btab = Table.read('/home/dnidever/projects/delvered/data/delvemc_bricks_0.25deg.fits.gz')
    for c in btab.colnames: btab[c].name = c.lower()
    btab['brickname'] = [b.strip() for b in btab['brickname']]
    dt = [('name',str,50),('dir',str,300),('cmd',str,100),('outfile',str,50),
          ('logfile',str,300),('errfile',str,300)]
    tasks = np.zeros(len(btab),dtype=np.dtype(dt))
    count = 0
    for i in range(len(btab)):
        brickname = btab['brickname'][i]
        bdir = '/net/dl2/dnidever/delve/bricks/'
        bdir = os.path.join(bdir,brickname[:4],brickname)
        print(i+1,brickname)
        if os.path.exists(bdir):
            outfile = os.path.join(bdir,brickname+'_epochdata.fits')
            if os.path.exists(outfile) and clobber==False:
                print(outfile,'exists and clobber not set')
                continue
            if nmulti<=1:
                epochs(brickname,clobber=clobber)
            else:
                cmd = 'python -c "from delvered import summary; summary.epochs('
                cmd += "'"+brickname+"'"
                if clobber:
                    cmd += ',clobber)"'
                else:
                    cmd += ')"'
                dir1 = '/net/dl2/dnidever/delve/bricks/summary/epochs/'
                logfile = os.path.join(dir1,brickname+'_epochs.log')
                tasks['name'][count] = brickname
                tasks['dir'][count] = dir1
                tasks['cmd'][count] = cmd  
                tasks['outfile'][count] = brickname+'_epochs'
                tasks['logfile'][count] = logfile
                tasks['errfile'][count] = logfile.replace('.log','.err')
                count += 1
        else:
            print(' no directory for',brickname)
    # trim tasks
    if count > 0 and count < len(tasks):
        tasks = tasks[:count]

    # use job_daemon
    if nmulti > 1:
        jobs = jd.job_daemon(tasks['cmd'],tasks['dir'],
                             prefix='epochs',hyperthread=True,idle=False,
                             nmulti=nmulti,nobreak=True)

    import pdb; pdb.set_trace()

def epochs(brickname,clobber=False):
    """ Get information for the temporal/epoch information """

    # Load each brick and 
    # nepochs per object, and time baseline per object

    bdir = '/net/dl2/dnidever/delve/bricks/'
    bdir = os.path.join(bdir,brickname[:4],brickname)

    outfile = os.path.join(bdir,brickname+'_epochdata.fits')
    if os.path.exists(outfile) and clobber==False:
        print(outfile,'already exists and clobber not set')
        return

    jmfile = os.path.join(bdir,brickname+'_joint_meas.fits.gz')
    mfile = os.path.join(bdir,brickname+'_expforced.fits.gz')
    if os.path.exists(jmfile):
        measfile = jmfile
    elif os.path.exists(jmfile)==False and os.path.exists(mfile):
        measfile = mfile
    else:
        #files = glob(bdir+'/*')
        #print(files)
        #import pdb; pdb.set_trace()
        return
    print('Loading',measfile)
    meas = Table.read(measfile)
    for c in meas.colnames: meas[c].name = c.lower()
    meas['objid'] = [o.strip() for o in meas['objid']]
    index = dln.create_index(meas['objid'])
    nobj = len(index['value'])
    print(nobj,'objects')
    dt = [('objid',str,20),('ndet',np.int32),('dt',np.float32)]
    res = np.zeros(nobj,dtype=np.dtype(dt))
    for i in range(nobj):
        ind = index['index'][index['lo'][i]:index['hi'][i]+1]
        nind = len(ind)
        res['objid'][i] = index['value'][i]
        res['ndet'][i] = nind
        res['dt'][i] = np.max(meas['mjd'][ind])-np.min(meas['mjd'][ind])

    # save
    print('Writing to',outfile)
    Table(res).write(outfile,overwrite=True)

def allepochscollate():
    """ Collate all of the epoch summary information. """

    btab = Table.read('/home/dnidever/projects/delvered/data/delvemc_bricks_0.25deg.fits.gz')
    for c in btab.colnames: btab[c].name = c.lower()
    btab['brickname'] = [b.strip() for b in btab['brickname']]
    tab = []
    hist_ndet = np.zeros(5000,int)
    hist_dt = np.zeros(5000,int)
    for i in range(len(btab)):
    #for i in range(1000):
        brickname = btab['brickname'][i]
        bdir = '/net/dl2/dnidever/delve/bricks/'
        bdir = os.path.join(bdir,brickname[:4],brickname)
        print(i+1,brickname)
        if os.path.exists(bdir):
            outfile = os.path.join(bdir,brickname+'_epochdata.fits')
            if os.path.exists(outfile)==False:
                print(outfile,'not found')
                continue
            tab1 = Table.read(outfile)
            print(' ',len(tab1))
            hist_ndet1,x1 = np.histogram(tab1['ndet'],bins=np.arange(5001))
            hist_dt1,x2 = np.histogram(tab1['dt'],bins=np.arange(5001))
            hist_ndet += hist_ndet1
            hist_dt += hist_dt1
            #import pdb; pdb.set_trace()
            #tab.append(tab1)
    import pdb; pdb.set_trace()
    #tab = vstack(tab)
    #tab.write('/net/dl2/dnidever/delve/bricks/summary/bricks_epochs.fits',overwrite=True)
    np.save('/net/dl2/dnidever/delve/bricks/summary/hist_ndet.npy',hist_ndet)
    np.save('/net/dl2/dnidever/delve/bricks/summary/hist_dt.npy',hist_dt)
    return hist_ndet,hist_dt


def metacollate():
    """ Collate all of the brick meta information. """

    btab = Table.read('/home/dnidever/projects/delvered/data/delvemc_bricks_0.25deg.fits.gz')
    for c in btab.colnames: btab[c].name = c.lower()
    btab['brickname'] = [b.strip() for b in btab['brickname']]
    meta = []
    for i in range(len(btab)):
    #for i in range(300):
        brickname = btab['brickname'][i]
        bdir = '/net/dl2/dnidever/delve/bricks/'
        bdir = os.path.join(bdir,brickname[:4],brickname)
        print(i+1,brickname)
        if os.path.exists(bdir):
            metafile = os.path.join(bdir,brickname+'_joint_meta.fits')
            jmetafile = os.path.join(bdir,brickname+'_joint_meta.fits')
            if os.path.exists(metafile)==False and os.path.exists(jmetafile)==False:
                print('no meta files found')
                continue
            if os.path.exists(jmetafile):
                meta1 = Table.read(jmetafile)
            else:
                meta1 = Table.read(metafile)
            print(' ',len(meta1))
            _,ui = np.unique(meta1['EXPNUM'],return_index=True)
            meta1 = meta1[ui]
            #import pdb; pdb.set_trace()
            meta.append(meta1)
    #import pdb; pdb.set_trace()
    meta = vstack(meta)
    print('Writing to /net/dl2/dnidever/delve/bricks/summary/allmeta.fits')
    meta.write('/net/dl2/dnidever/delve/bricks/summary/allmeta.fits',overwrite=True)


def allphoterrsummary(clobber=False,nmulti=5,group=10):
    """ Run through all of the bricks and getting the cmd summary data """
    btab = Table.read('/home/dnidever/projects/delvered/data/delvemc_bricks_0.25deg.fits.gz')
    for c in btab.colnames: btab[c].name = c.lower()
    btab['brickname'] = [b.strip() for b in btab['brickname']]
    dt = [('name',str,50),('dir',str,300),('cmd',str,200),('outfile',str,50),
          ('logfile',str,300),('errfile',str,300)]
    count = 0
    njobs = len(btab)//group + 1
    tasks = np.zeros(njobs,dtype=np.dtype(dt))
    for i in range(njobs):
        blist = []
        for j in range(group):
            blist.append(btab['brickname'][count])
            count += 1
            if count > len(btab)-1:
                break
            
        cmd = 'python -c "from delvered import summary; summary.photerrsummarymulti('
        cmd += "['" + "','".join(blist) + "']"
        if clobber:
            cmd += ',clobber=True)"'
        else:
            cmd += ')"'
        dir1 = '/net/dl2/dnidever/delve/bricks/summary/photerr/'
        #logfile = os.path.join(dir1,brickname+'_epochs.log')
        tasks['name'][i] = 'group'+str(i+1)
        tasks['dir'][i] = dir1
        tasks['cmd'][i] = cmd  
        #tasks['outfile'][i] =
        #tasks['logfile'][i] = logfile
        #tasks['errfile'][i] = logfile.replace('.log','.err')

    #import pdb; pdb.set_trace()

    # use job_daemon
    jobs = jd.job_daemon(tasks['cmd'],tasks['dir'],
                         prefix='photerrsum',hyperthread=True,idle=False,
                         nmulti=nmulti,nobreak=True)

    import pdb; pdb.set_trace()


def photerrsummarymulti(bricknames,clobber=False):
    """ Thin wrapper to run multiple bricks at a time """
    for i in range(len(bricknames)):
        print(i+1,bricknames[i])
        photerrsummary(bricknames[i],clobber=clobber)

def photerrsummary(brickname,clobber=False):
    """ Photometric error summary information """

    # Load each brick and 

    bdir = '/net/dl2/dnidever/delve/bricks/'
    bdir = os.path.join(bdir,brickname[:4],brickname)

    outfile = os.path.join(bdir,brickname+'_photerrsummary.fits')
    if os.path.exists(outfile) and clobber==False:
        print(outfile,'already exists and clobber not set')
        return

    jmfile = os.path.join(bdir,brickname+'_joint_object.fits.gz')
    mfile = os.path.join(bdir,brickname+'_object.fits.gz')
    if os.path.exists(jmfile):
        objfile = jmfile
    elif os.path.exists(jmfile)==False and os.path.exists(mfile):
        objfile = mfile
    else:
        #files = glob(bdir+'/*')
        #print(files)
        #import pdb; pdb.set_trace()
        return
    print('Loading',objfile)
    obj = Table.read(objfile)
    for c in obj.colnames: obj[c].name = c.lower()
    print(len(obj),'objects')

    bands = ['u','g','r','i','z','y']
    dt = [('band',str,20),('count',np.int32,(169,119))]
    res = np.zeros(len(bands),dtype=np.dtype(dt))
    #xbins = np.arange(10,27,0.1)
    xbins = np.arange(10,30,0.1)
    ybins = np.arange(-4,2,0.05)
    for i in range(len(bands)):
        band = bands[i]
        res['band'][i] = band
        if band+'mag' not in obj.colnames: continue
        mag = obj[band+'mag']
        err = obj[band+'err']
        gd, = np.where((mag > 5) & (mag < 50) & np.isfinite(mag))
        #               (np.abs(obj['sharp'])<2) & (obj['chi']<5) & (obj['prob']>0.2))
        if len(gd)>0:
            out = binned_statistic_2d(mag[gd],np.log10(err[gd]),mag[gd],
                                      statistic='count',bins=(xbins,ybins))
            result,xedge,yedge,binnumber = out
            res['count'][i] = result

    # save
    print('Writing to',outfile)
    Table(res).write(outfile,overwrite=True)

def photerrcollate():
    """ Collate all of the brick photerr information. """

    btab = Table.read('/home/dnidever/projects/delvered/data/delvemc_bricks_0.25deg.fits.gz')
    for c in btab.colnames: btab[c].name = c.lower()
    btab['brickname'] = [b.strip() for b in btab['brickname']]
    res = []
    for i in range(len(btab)):
    #for i in range(300):
        brickname = btab['brickname'][i]
        bdir = '/net/dl2/dnidever/delve/bricks/'
        bdir = os.path.join(bdir,brickname[:4],brickname)
        print(i+1,brickname)
        if os.path.exists(bdir)==False:
            continue
        jmfile = os.path.join(bdir,brickname+'_joint_object.fits.gz')
        mfile = os.path.join(bdir,brickname+'_object.fits.gz')
        if os.path.exists(jmfile)==False and os.path.exists(mfile)==False:
            print('no object files found')
            continue
        outfile = os.path.join(bdir,brickname+'_photerrsummary.fits')
        if os.path.exists(outfile)==False:
            print(outfile,'not found')
            return
        res1 = Table.read(outfile)
        res.append(res1)
    #import pdb; pdb.set_trace()
    res = vstack(res)
    print('Writing to /net/dl2/dnidever/delve/bricks/summary/allphoterr.fits')
    res.write('/net/dl2/dnidever/delve/bricks/summary/allphoterr.fits',overwrite=True)

def allcmdsummary(clobber=False,nmulti=5,group=10):
    """ Run through all of the bricks and getting the cmd summary data """
    btab = Table.read('/home/dnidever/projects/delvered/data/delvemc_bricks_0.25deg.fits.gz')
    for c in btab.colnames: btab[c].name = c.lower()
    btab['brickname'] = [b.strip() for b in btab['brickname']]
    dt = [('name',str,50),('dir',str,300),('cmd',str,200),('outfile',str,50),
          ('logfile',str,300),('errfile',str,300)]
    count = 0
    njobs = len(btab)//group + 1
    tasks = np.zeros(njobs,dtype=np.dtype(dt))
    for i in range(njobs):
        blist = []
        for j in range(group):
            blist.append(btab['brickname'][count])
            count += 1
            if count > len(btab)-1:
                break
            
        cmd = 'python -c "from delvered import summary; summary.cmdsummarymulti('
        cmd += "['" + "','".join(blist) + "']"
        if clobber:
            cmd += ',clobber=True)"'
        else:
            cmd += ')"'
        dir1 = '/net/dl2/dnidever/delve/bricks/summary/cmdsummary/'
        #logfile = os.path.join(dir1,brickname+'_epochs.log')
        tasks['name'][i] = 'group'+str(i+1)
        tasks['dir'][i] = dir1
        tasks['cmd'][i] = cmd  
        #tasks['outfile'][i] =
        #tasks['logfile'][i] = logfile
        #tasks['errfile'][i] = logfile.replace('.log','.err')

    #import pdb; pdb.set_trace()

    # use job_daemon
    jobs = jd.job_daemon(tasks['cmd'],tasks['dir'],
                         prefix='cmdsum',hyperthread=True,idle=False,
                         nmulti=nmulti,nobreak=True)

    import pdb; pdb.set_trace()

def cmdsummarymulti(bricknames,clobber=False):
    """ Thin wrapper to run multiple bricks at a time """
    for i in range(len(bricknames)):
        print(i+1,bricknames[i])
        cmdsummary(bricknames[i],clobber=clobber)

def cmdsummary(brickname,clobber=False):
    """ Get CMD information for a brick """

    # Load each brick and 

    bdir = '/net/dl2/dnidever/delve/bricks/'
    bdir = os.path.join(bdir,brickname[:4],brickname)

    outfile = os.path.join(bdir,brickname+'_cmdforcedsummary.fits')
    #outfile = os.path.join(bdir,brickname+'_cmdsummary.fits')
    if os.path.exists(outfile) and clobber==False:
        print(outfile,'already exists and clobber not set')
        return

    jmfile = os.path.join(bdir,brickname+'_joint_object.fits.gz')
    mfile = os.path.join(bdir,brickname+'_object.fits.gz')
    if os.path.exists(jmfile):
        objfile = jmfile
    elif os.path.exists(jmfile)==False and os.path.exists(mfile):
        objfile = mfile
    else:
        #files = glob(bdir+'/*')
        #print(files)
        #import pdb; pdb.set_trace()
        return
    objfile = mfile
    print('Loading',objfile)
    obj = Table.read(objfile)
    for c in obj.colnames: obj[c].name = c.lower()
    print(len(obj),'objects')

    xbins = np.arange(-4,6,0.05)
    ybins = np.arange(10,30,0.1)
    shape = (len(xbins)-1,len(ybins)-1)
    dt = [('count',np.int32,shape),('stars',np.int32,shape),
          ('variables',np.int32,shape),('variables10',np.int32,shape)]
    res = np.zeros(1,dtype=np.dtype(dt))

    gd, = np.where((obj['gmag']>5) & (obj['gmag']<50) &
                   (obj['imag']>5) & (obj['imag']<50))
    if len(gd)>0:
        out = binned_statistic_2d(obj['gmag'][gd]-obj['imag'][gd],obj['gmag'][gd],obj['gmag'][gd],
                                  statistic='count',bins=(xbins,ybins))
        result,xedge,yedge,binnumber = out
        res['count'][0] = result

    gds, = np.where((obj['gmag']>5) & (obj['gmag']<50) &
                    (obj['imag']>5) & (obj['imag']<50) &
                    (np.abs(obj['sharp'])<1) & (obj['chi']<3) & (obj['prob']>0.5))
    if len(gds)>0:
        out = binned_statistic_2d(obj['gmag'][gds]-obj['imag'][gds],obj['gmag'][gds],obj['gmag'][gds],
                                  statistic='count',bins=(xbins,ybins))
        result,xedge,yedge,binnumber = out
        res['stars'][0] = result

    # Do variables as well
    if 'variable10sig' in obj.colnames:
        gdvar, = np.where((obj['gmag']>5) & (obj['gmag']<50) &
                          (obj['imag']>5) & (obj['imag']<50) &
                          (obj['variable10sig']==1))
        if len(gdvar)>0:
            out = binned_statistic_2d(obj['gmag'][gdvar]-obj['imag'][gdvar],obj['gmag'][gdvar],
                                      obj['gmag'][gdvar],statistic='count',bins=(xbins,ybins))
            result,xedge,yedge,binnumber = out
            res['variables'][0] = result

        gdvar10, = np.where((obj['gmag']>5) & (obj['gmag']<50) &
                            (obj['imag']>5) & (obj['imag']<50) &
                            (obj['variable10sig']==1) & (obj['ndet']>=10))
        if len(gdvar10)>0:
            out = binned_statistic_2d(obj['gmag'][gdvar10]-obj['imag'][gdvar10],obj['gmag'][gdvar10],
                                      obj['gmag'][gdvar10],statistic='count',bins=(xbins,ybins))
            result,xedge,yedge,binnumber = out
            res['variables10'][0] = result

    # save
    print('Writing to',outfile)
    Table(res).write(outfile,overwrite=True)

def cmdcollate():
    """ Collate all of the brick cmd information. """

    btab = Table.read('/home/dnidever/projects/delvered/data/delvemc_bricks_0.25deg.fits.gz')
    for c in btab.colnames: btab[c].name = c.lower()
    btab['brickname'] = [b.strip() for b in btab['brickname']]
    res = []
    for i in range(len(btab)):
    #for i in range(300):
        brickname = btab['brickname'][i]
        bdir = '/net/dl2/dnidever/delve/bricks/'
        bdir = os.path.join(bdir,brickname[:4],brickname)
        print(i+1,brickname)
        if os.path.exists(bdir)==False:
            continue
        jmfile = os.path.join(bdir,brickname+'_joint_object.fits.gz')
        mfile = os.path.join(bdir,brickname+'_object.fits.gz')
        if os.path.exists(jmfile)==False and os.path.exists(mfile)==False:
            print('no object files found')
            continue
        outfile = os.path.join(bdir,brickname+'_cmdsummary.fits')
        if os.path.exists(outfile)==False:
            print(outfile,'not found')
            continue
        try:
            res1 = Table.read(outfile)
            res.append(res1)
        except:
            print('Problem load',outfile)
    #import pdb; pdb.set_trace()
    res = vstack(res)
    print('Writing to /net/dl2/dnidever/delve/bricks/summary/allcmdsummary.fits')
    res.write('/net/dl2/dnidever/delve/bricks/summary/allcmdsummary.fits',overwrite=True)
