#!/usr/bin/env python

import os
import time
import numpy as np
import socket
import logging
import tempfile
import shutil
from glob import glob
from datetime import datetime
from dlnpyutils import utils as dln, coords
import photred
from photred import io,daomatch
import photred.allframe as alf
import dill as pickle
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table
from astropy.coordinates import SkyCoord
from dustmaps.sfd import SFDQuery
import subprocess

def make_brick_wcs(brickstr):
    """
    This creates a brick WCS given the brick structure
    
    Each brick covers 0.25 x 0.25 deg at a pixel scale of 0.262"/pix
    with 3600x3600 pixels that gives us an overlap of ~82 pixels on
    each side.  The *UNIQUE* brick area (0.25x0.25 deg) is defined by
    the BRAMIN/MAX and BDECMIN/MAX keywords.
    """

    # Make the tiling file
    #---------------------
    nx = 3600
    ny = 3600
    step = 0.262 / 3600    # 0.262" per pixel, DECam pixel scale
    xref = nx/2
    yref = ny/2

    #  Make the header as well
    tilehead = fits.Header()
    tilehead['NAXIS1'] = nx
    tilehead['CDELT1'] = step
    tilehead['CRPIX1'] = xref+1
    tilehead['CRVAL1'] = brickstr['ra']
    tilehead['CTYPE1'] = 'RA---TAN'
    tilehead['NAXIS2'] = ny
    tilehead['CDELT2'] = step
    tilehead['CRPIX2'] = yref+1
    tilehead['CRVAL2'] = brickstr['dec']
    tilehead['CTYPE2'] = 'DEC--TAN'
    tilehead['BRAMIN'] = brickstr['ra1'],'RA min of unique brick area'
    tilehead['BRAMAX'] = brickstr['ra2'],'RA max of unique brick area'
    tilehead['BDECMIN'] = brickstr['dec1'],'DEC min of unique brick area'
    tilehead['BDECMAX'] = brickstr['dec2'],'DEC max of unique brick area'
    wcs = WCS(tilehead)
    
    # Create the TILE structure
    tilestr = {'type':'WCS','naxis':np.array([nx,ny]),'cdelt':np.array([step,step]).astype(float),
               'crpix':np.array([xref+1,yref+1]).astype(float),
               'crval':np.array([brickstr['ra'],brickstr['dec']]).astype(float),'ctype':['RA--TAN','DEC--TAN'],
               'head':tilehead,'wcs':wcs,'xrange':[0,nx],'yrange':[0,ny],'nx':nx,'ny':ny}

    return tilestr

def forcebrick(brick,scriptsdir=None,irafdir=None,workdir=None,redo=False,
               update=False,logfile=None,delvedir=None,delvereddir=None):
    """
    Process a single DELVE brick and perform ALLFRAME FORCED photometry 
 
    Parameters
    ----------
    brick : str
      The DELVE brick name, e.g. 1234m045 
 
    By D. Nidever  August 2019 
    Translated to Python by D. Nidever,  April 2022
    """
     
    # Limit the number of threads 
    #CPU,tpool_nthreads=4 
     
    # This bricks pre-processing script gets DELVE and community MC data ready 
    # to run PHOTRED ALLFRAME on it. 
     
    t0 = time.time()
    curdir = os.path.abspath(os.curdir) 
     
    # Defaults
    if delvedir is not None:
        if delvedir.endswith('/')==False:
            delvedir += '/'        
    else: 
        delvedir = '/net/dl2/dnidever/delve/' 
    if os.path.exists(delvedir)==False:
        os.makedirs(delvedir)
    if delvereddir is not None:
        if delvereddir.endswith('/')==False:
            delvereddir += '/'
    else: 
        delvereddir = '/'.join(os.path.abspath(__file__).split('/')[:-3])+'/'
    if scriptsdir is not None:
        if scriptsdir.endswith('/')==False:
            scriptsdir += '/'
    else: 
        scriptsdir = '/'.join(photred.__file__.split('/')[:-3])+'/scripts/'
    if irafdir is not None:
        if irafdir.endswith('/')==False:
            irafdir += '/'
    else: 
        irafdir = os.path.expanduser('~/iraf/')
    tempdir = '/tmp/'
    if workdir is None:
        host = socket.gethostname()
        hostname = host.split('.')[0]
        workdir = '/data0/dnidever/delve/' 
        if os.path.exists(workdir)==False:
            os.makedirs(workdir)
    # Exposures directory
    expdir = delvedir+'exposures/'
    # Bricks directory
    brickdir = delvedir+'bricks/' 
    if os.path.exists(brickdir)==False:
        os.makedirs(brickdir)
    logsdir = brickdir+'logs/' 
    if os.path.exists(logsdir)==False:
        os.makedirs(logsdir)
     
    # Start the logfile 
    #------------------ 
    # format is delvered_forcebrick.brick.host.DATETIME.log
    host = socket.gethostname()
    hostname = host.split('.')[0]
    logtime = datetime.now().strftime("%Y%m%d%H%M%S") 
    # Set up logging to screen and logfile
    logFormatter = logging.Formatter("%(asctime)s [%(levelname)-5.5s]  %(message)s")
    logger = logging.getLogger() 
    while logger.hasHandlers(): # some existing loggers, remove them   
        logger.removeHandler(logger.handlers[0]) 
    logger = logging.getLogger()
    logtime = datetime.now().strftime("%Y%m%d%H%M%S")
    if logfile is None:
        logfile = logsdir+'delvered_brick.'+brick+'.'+hostname+'.'+logtime+'.log' 
    if os.path.exists(logfile): os.remove(logfile)
    fileHandler = logging.FileHandler(logfile)
    fileHandler.setFormatter(logFormatter)
    logger.addHandler(fileHandler)
    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    logger.addHandler(consoleHandler)
    logger.setLevel(logging.NOTSET)
     
    # Chip name and CCDNUM relations
    decam = Table.read(delvereddir+'data/decam.txt',format='ascii')
    for n in decam.colnames: decam[n].name = n.lower()  # lowercase column names
     
    # Load the brick information
    brickstr = Table.read(delvereddir+'data/delvemc_bricks_0.25deg.fits.gz')
    for n in brickstr.colnames: brickstr[n].name = n.lower()  # lowercase column names

    # Get the brick information 
    bind, = np.where(brickstr['brickname'] == brick) 
    if len(bind) == 0: 
        logger.info(brick+' not in DELVE-MC brick list')
        return 
    brickstr1 = brickstr[bind[0]]
     
    # Subdirectory is the first 4 digits of the brickname, e.g., 0952 of 0952m462, the RA portion of the name 
    subdir = brickdir+brick[0:4]+'/' 
    if os.path.exists(subdir)==False:
        os.makedirs(subdir)
    # Create brick directory if it doesn't exist yet 
    bdir = subdir+brick+'/' 
    if os.path.exists(bdir)==False:
        os.makedirs(bdir)
    logfile = bdir+brick+'.'+logtime+'.log' 
     
    # Check the output file 
    photfile = bdir+brick+'.fits'
    if os.path.exists(photfile+'.gz') and redo==False and update==False:
        logger.info(photfile+'.gz EXISTS and redo or update NOT set')
        return 
     
    # DECam imager 
    thisimager = {'telescope':'BLANCO','instrument':'DECam','namps':62,'separator':'_'} 
     
    # Print information 
    logger.info('--------------------------------------------------------------------')
    logger.info(' Run ALLFRAME forced photometry on DELVE Brick = '+brick)
    logger.info('--------------------------------------------------------------------')
    logger.info('DELVEDIR = '+delvedir)
    logger.info('SCRIPTSDIR = '+scriptsdir)
    logger.info('IRAFDIR = '+irafdir)
    logger.info('EXPOSUREDIR = '+expdir)
    logger.info('WORKINGDIR = '+workdir)
    logger.info('HOST = '+host)
    logger.info(datetime.now().strftime("%a %b %d %H:%M:%S %Y")) 
     
    # The photed.setup file 
    setup = ['##### REQUIRED #####',
             'scriptsdir  '+scriptsdir,
             'irafdir     '+irafdir,
             'telescope   Blanco',
             'instrument  DECAM',
             'observatory CTIO',
             'nmulti      10',
             'nmulti_wcs       40',
             'nmulti_daophot   30',
             'nmulti_allframe  10',
             'filtref     g,i,r,z,u',
             'trans       delve.trans',
             '##### OPTIONAL #####',
             'sepfielddir  1',
             'sepchipdir   1',
             'keepmef      0',
             'catformat    FITS',
             'workdir      '+workdir,
             'clean        1',
             'skipcheck    1',
             'redo         0',
             'wcsrefname   GAIADR2',
             'searchdist   20',
             '#wcsrmslim   1.0',
             'hyperthread   1',
             'daopsfva      1',
             'daofitradfwhm 1.0',
             'psfcomsrc     0',
             'psfcomglobal  0',
             'psfcomgauss   0',
             '#mchmaxshift  50.0',
             'finditer      1',       # 2->1 on 7/22/20 
             'alfdetprog  sextractor',
             '#alfnocmbimscale 0',
             'alftrimcomb   0',
             '#ddo51radoffset  1',
             'cmbforce      1',
             'keepinstr     1',
             'avgmag        1',
             'avgonlymag    0',
             'todered       u,g,r,i,z,g-i',
             '##### STAGES #####',
             '#rename',
             '#split',
             '#wcs',
             '#daophot',
             '#zeropoint',
             '#match',
             ' allframe',
             '#apcor',
             '#astrom',
             '#calib',
             '#combine',
             '#deredden',
             '#save',
             '#html']
    dln.writelines(bdir+'photred.setup',setup)
    setup = io.readsetup(setupdir=bdir)
     
    # Print out some information about the brick 
    logger.info('RA = %.5f' % brickstr1['ra'])
    logger.info('DEC = %.3f' % brickstr1['dec'])
    logger.info('RA range  = [ %.5f,%.5f ]' % (brickstr1['ra1'],brickstr1['ra2']))
    logger.info('DEC range = [ %.5f,%.5f ]' % (brickstr1['dec1'],brickstr1['dec2']))                
     
    # Get the Brick WCS information 
    tilestr = make_brick_wcs(brickstr1)
     
    # Step 1: Get the list of exposures/chips that overlap this brick 
    #---------------------------------------------------------------- 
    logger.info('Step 1: Get list of chip files that overlap this brick')
    #logger.info('Getting chips that overlap this brick' 
    cenra = brickstr1['ra']
    cendec = brickstr1['dec']
    tid,tmpfile = tempfile.mkstemp(prefix="tmp",suffix='.fits',dir=tempdir)
    dln.touch(tmpfile)
    for f in [tmpfile]:
        if os.path.exists(f): os.remove(f)
    # /noshell causes problems on gp09 because it gives python2 instead of python3
    out = subprocess.check_output(delvereddir+'bin/query_delvered_summary_table '+str(cenra)+' '+str(cendec)+' '+tmpfile+' --lim 0.5',shell=True)
    exists = os.path.exists(tmpfile)
    size = os.path.getsize(tmpfile)
    if size == 0: 
        logger.info('No overlapping chips found')
        if os.path.exists(tmpfile): os.remove(tmpfile)
        return
    chstr = Table.read(tmpfile)
    if os.path.exists(tmpfile): os.remove(tmpfile)
    nchstr = len(chstr) 
    logger.info('Found '+str(nchstr)+' overlapping chips within 0.5 deg of brick center')
    logger.info(datetime.now().strftime("%a %b %d %H:%M:%S %Y")) 
     
    # Make sure the chips are unique,  some were duplicated on SMASH nights 
    chid = np.char.array(chstr['expnum'].astype(str))+'-'+np.char.array(chstr['chip'].astype(str))
    uchstr,ui = np.unique(chid,return_index=True)
    chstr = chstr[ui] 
    nchstr = len(chstr) 
     
    # Do more rigorous overlap checking 
    #  the brick region with overlap 
    print('Performing more rigorous overlap checking')

    bcoo = tilestr['wcs'].pixel_to_world([0,tilestr['nx']-1,tilestr['nx']-1,0],[0,0,tilestr['ny']-1,tilestr['ny']-1])
    bvra = bcoo.ra.deg
    bvdec = bcoo.dec.deg
    #HEAD_XYAD,tilestr.head,[0,tilestr.nx-1,tilestr.nx-1,0],[0,0,tilestr.ny-1,tilestr.ny-1],bvra,bvdec,/deg
    olap = np.zeros(nchstr,bool)
    vxarr = np.zeros((nchstr,4),float)
    vyarr = np.zeros((nchstr,4),float)    
    for i in range(nchstr): 
        if (i % 100 == 0) and (i > 0): 
            print(i)
        hd1 = io.readfile(chstr['file'][i],header=True)
        nx = hd1['naxis1']
        ny = hd1['naxis2']
        w1 = WCS(hd1)
        vcoo = w1.pixel_to_world([0,nx-1,nx-1,0],[0,0,ny-1,ny-1])
        vra = vcoo.ra.deg
        vdec = vcoo.dec.deg
        olap[i] = coords.doPolygonsOverlap(bvra,bvdec,vra,vdec) 
        vx,vy = tilestr['wcs'].world_to_pixel(vcoo)
        vxarr[i,:] = vx 
        vyarr[i,:] = vy 
    # Require at least a 2 pixel overlap in X and Y 
    g, = np.where((olap==True) & (np.max(vxarr,axis=1) >= 2) & (np.max(vyarr,axis=1) >= 2) &
                  (np.min(vxarr,axis=1) <= tilestr['nx']-3) & (np.min(vyarr,axis=1) <= tilestr['ny']-3))
    if len(g) == 0: 
        logger.info('No chips overlap this brick')
        return 
    logger.info(str(len(g))+' chips overlap this brick')
    chstr = chstr[g] 
    nchstr = len(chstr)
     
    # APPLY cuts on the exposures 
    #----------------------------- 
    logger.info('Applying quality, zero-point, filter and exptime cuts')
     
    # Zero-point structure, from NSC
    zpstr = np.zeros(7,dtype=np.dtype([('instrument',np.str,10),('filter',np.str,10),('amcoef',(float,2)),('thresh',float)]))
    #zpstr = replicate({instrument:'',filter:'',amcoef:fltarr(2),thresh:0.5},7)
    zpstr['instrument'] = 'c4d'
    zpstr['filter'] = ['u','g','r','i','z','Y','VR'] 
    zpstr['amcoef'][0] = [-1.60273, -0.375253]  # c4d-u 
    zpstr['amcoef'][1] = [0.277124, -0.198037]  # c4d-g 
    zpstr['amcoef'][2] = [0.516382, -0.115443]  # c4d-r 
    zpstr['amcoef'][3] = [0.380338, -0.067439]  # c4d-i 
    zpstr['amcoef'][4] = [0.074517, -0.067031]  # c4d-z 
    zpstr['amcoef'][5] = [-1.07800, -0.060014]  # c4d-Y 
    zpstr['amcoef'][6] = [1.111859, -0.083630]  # c4d-VR 
     
    # Convert to additive zero-point as used in NSC 
    zpterm = -chstr['calib_zpterm']
    # Fix early DES exposures that used different units/gain 
    gdes, = np.where(chstr['gain'] < 2)
    if len(gdes) > 0: 
        zpterm[gdes] -= 1.55 
    # global zpterm and airmass correction 
    for i in range(len(zpstr)): 
        ind, = np.where(chstr['filter']==zpstr['filter'][i])
        if len(ind)>0: 
            zpterm[ind] -= np.polyval(np.flip(zpstr['amcoef'][i]),chstr['airmass'][ind])
                     
    fwhmthresh = 2.0  # seeing 2.0" threshold 
    filt = np.char.array(chstr['filter'].astype(str)).ljust(1)  # first characater
    gdch, = np.where((chstr['fwhm']*chstr['pixscale'] <= fwhmthresh) & (chstr['exptime'] >= 90.) & (zpterm >= -0.5) &
                     np.isfinite(zpterm) & np.isfinite(chstr['apcor']) &
                     ((filt== 'u') | (filt=='g') | (filt=='r') | (filt=='i') | (filt=='z') | (filt=='Y')))
    if len(gdch) == 0: 
        logger.info('No chips passed the cuts')
        return 
    logger.info(str(len(gdch))+' chips passed the cuts')
    chstr = chstr[gdch] 
    nchstr = len(chstr)
    chstr['file'] = np.char.array(chstr['file'].astype(str)).strip()
    chstr['base'] = np.char.array(chstr['base'].astype(str)).strip()
     
    # Check if we need to update
    if update and redo==False:
        # Load previous meta 
        metafile = bdir+brick+'_meta.fits'
        if os.path.exists(metafile):
            meta0 = Table.read(metafile,1)
            meta0['base'] = str(meta0['base'])
            ind1,ind2 = dln.match(meta0['base'],chstr['base'])
            if nchstr == len(meta0) and nchstr == nmatch: 
                logger.info('Nothing to UPDATE') 
                return 
            # Moving previous catalogs to a backup
            oldfiles = glob(bdir+brick+'*fits*')
            oldfiles += glob(bdir+brick+'allframe.opt')
            oldfiles += glob(bdir+brick+'default.*')
            if len(oldfiles)>0: 
                bakdir = bdir+'bak'+logtime
                logger.info('Backing up old files to '+bakdir)
                os.makedirs(bakdir)
                for f in oldfiles:
                    shutil.move(f,bakdir)
            else:
                logger.info('No old files to backup')
            logger.info('There are new exposures to include.  UPDATING')
         
                  
    # Create temporary local directory to perform the work/processing 
    #  copy everything to it 
    if os.path.exists(workdir)==False:
        os.makedirs(workdir)
                     
    procdir = tempfile.mkdtemp(prefix="brk",dir=workdir)
    procdir += '/'
    os.chmod(procdir,0o755)
    logger.info('Working in temporary directory '+procdir)
    logger.info(datetime.now().strftime("%a %b %d %H:%M:%S %Y")) 
                    
    # Copy the files 
    #---------------- 
    logger.info('Copying over the necessary files to '+procdir)
    # Chip loop 
    for i in range(nchstr): 
        logger.info('Copying over files for '+chstr['file'][i])
        # New filename 
        #chstr[i].newbase = chstr[i].base 
        odir = os.path.dirname(chstr['file'][i])+'/'
        obase = photred.utils.fitsext(chstr['file'][i],basename=True)
        obase = obase[0]
        # Need symlinks to .psf, .als 
        if os.path.exists(odir+obase+'.psf') == False: 
            logger.info(odir+obase+'.psf NOT FOUND.  Skipping this chip') 
            continue
        if os.path.exists(odir+obase+'.als') == False: 
            logger.info(odir+obase+'.als NOT FOUND.  Skipping this chip')
            continue
        for e in ['.psf','.als','.ap','.opt','.als.opt','.log']:
            if os.path.exists(procdir+chstr[i].base+e): os.remove(procdir+chstr[i].base+e)
        for e in ['.psf','.als','.ap','.opt','.als.opt','.log']:
            if os.path.exists(procdir+'/'+obase+e): os.remove(procdir+'/'+obase+e)
            shutil.copyfile(odir+obase+e,procdir)
            os.chmod(procdir+chstr['base'][i]+e,0o755)  # make sure they are writable 
        # Copy the fits, fits resource file and header files locally 
        if os.path.exists(odir+obase+'.fits') == False: 
            logger.info(odir+obase+'.fits NOT FOUND.  Skipping this chip')
            continue
        if os.path.exists(odir+obase+'.fits'): os.remove(odir+obase+'.fits')
        shutil.copyfile(odir+obase+'.fits',procdir)
        if os.path.exists(odir+'.'+obase+'.fits'):
            if os.path.exists(procdir+'.'+obase+'.fits'): os.remove(procdir+'.'+obase+'.fits')
            shutil.copyfile(odir+'.'+obase+'.fits',procdir)
        if os.path.exists(odir+obase+'.fits.head'):
            if os.path.exists(procdir+obase+'.fits.head'): os.remove(procdir+obase+'.fits.head')
            shutil.copyfile(odir+obase+'.fits.head',procdir)
        os.chmod(procdir+chstr['base'][i]+'.fits',0o755)
         
    # Step 2: Run DAOMATCH_TILE.PRO on the files 
    #-------------------------------------------- 
    logger.info('Step 2: Matching up objects with DAOMATCH_TILE')
    logger.info(datetime.now().strftime("%a %b %d %H:%M:%S %Y"))
    os.chdir(procdir)
    groupstr = {'x0':0,'y0':0} 
    daomatch_tile(chstr['base']+'.als',tilestr,groupstr)
                  
    # Step 3: Run ALLFRAME 
    #---------------------- 
    logger.info('Step 3: Run ALLFRAME')
    logger.info(datetime.now().strftime("%a %b %d %H:%M:%S %Y"))                     
    # DO I NEED TO HAVE IT TRIM THE COMBINED IMAGE??? 
    mchbase = procdir+chstr[0].base 
    mchfile = mchbase+'.mch'# allframe needs absolute path 
    alf.allframe(mchfile,tile=tilestr,setupdir=bdir,scriptsdir=scriptsdir,irafdir=irafdir,
                 logfile=logfile,catformat='FITS',imager=thisimager,geocoef=0)
    magfile = chstr['base'][0]+'.mag' 
    if os.path.exists(magfile) == False: 
        logger.info(magfile+' NOT FOUND')
        return 
         
    # Step 4: Calculate coordinates 
    #------------------------------- 
    logger.info('Step 4: Adding coordinates')
    logger.info(datetime.now().strftime("%a %b %d %H:%M:%S %Y"))                     
    # Load the MCH file
    alsfiles,trans,magoff = io.readmch(mchfile)
    nalsfiles = len(alsfiles) 
    # Load the photometry file 
    instphot = io.readfile(magfile) 
    ninstphot = len(instphot) 
    logger.info('Nstars = '+str(ninstphot)) 
    # Converting to IDL X/Y convention, starting at (0,0) 
    # DAOPHOT has X/Y start at (1,1)
    ra,dec = tilestr['wcs'].pixel_to_world(instphot.x-1.0,instphot.y-1.0)
    #HEAD_XYAD,tilestr.head,instphot.x-1.0,instphot.y-1.0,ra,dec,/degree 
              
    # Step 5: Calibrating photometry with zero-points 
    #------------------------------------------------- 
    logger.info('Step 5: Calibrating photometry with zero-points')
    logger.info(datetime.now().strftime("%a %b %d %H:%M:%S %Y"))
    cmag = np.zeros((ninstphot,nchstr),float)+99.99
    cerr = np.zeros((ninstphot,nchstr),float)+9.99        
    # Chip loop 
    for i in range(nchstr): 
        # id, x, y, ra, dec, unsolved magnitudes, chi, sharp 
        imag = instphot['mag'+str(i+1)] 
        ierr = instphot['err'+str(i+1)]
        gdmag, = np.where(imag < 50)
        if len(gdmag) > 0: 
            # exptime, aperture correction, zero-point 
            # aperture correction is SUBTRACTIVE, makes it brighter 
            # ZPTERM is a SUBTRACTIVE constant offset 
            cmag[gdmag,i] = imag[gdmag] + 2.5*np.log10(chstr['exptime'][i]) - chstr['apcor'][i] - chstr['calib_zpterm'][i]
            # Add zero-point error in quadrature 
            cerr[gdmag,i] = np.sqrt(ierr[gdmag]**2+chstr['calib_zptermsig'][i]**2) 
    # Calculate average photometry per filter
    ufilt = np.unique(chstr['filter'])
    nufilt = len(ufilt)
    avgmag = np.zeros((ninstphot,nufilt),float)
    avgerr = np.zeros((ninstphot,nufilt),float)        
    ndet = np.zeros((ninstphot,nufilt),int)
    for i in range(nufilt): 
        gdf, = np.where(chstr['filter'] == ufilt[i]) 
        # Single exposure 
        if len(gdf) == 1: 
            avgmag[:,i] = cmag[:,gdf[0]] 
            avgerr[:,i] = cerr[:,gdf[0]] 
            ndet[:,i] = np.sum(cmag[:,gdf[0]] < 50) 
            # Multiple exposures 
        else: 
            # Loop through all of the exposures and add up the flux, totalwt, etc.
            totalwt = np.zeros(ninstphot,float)
            totalfluxwt = np.zeros(ninstphot,float)
            for k in range(ngdf): 
                gdmag, = np.where(cmag[:,gdf[k]] < 50) 
                if len(gdmag) > 0: 
                    totalwt[gdmag] += 1.0/cerr[gdmag,gdf[k]]**2 
                    totalfluxwt[gdmag] += 2.5118864**cmag[gdmag,gdf[k]] * (1.0/cerr[gdmag,gdf[k]]**2) 
                    ndet[gdmag,i] += 1
            newflux = totalfluxwt/totalwt 
            newmag = 2.5*np.log10(newflux) 
            newerr = np.sqrt(1.0/totalwt) 
            bdmag, = np.where(~np.isfinite(newmag)) 
            if len(bdmag) > 0: 
                newmag[bdmag] = 99.99 
                newerr[bdmag] = 9.99 
            avgmag[:,i] = newmag 
            avgerr[:,i] = newerr 
    # Measure scatter
    scatter = np.zeros((ninstphot,nufilt),float)+99.99
    for i in range(nufilt): 
        gdf, = np.where(chstr['filter'] == ufilt[i]) 
        if len(gdf) > 1:
            totaldiff = np.zeros(ninstphot,float)
            for k in range(ngdf): 
                gdmag, = np.where(cmag[:,gdf[k]] < 50) 
                if len(gdmag) > 0: 
                    totaldiff[gdmag] += (avgmag[gdmag,i]-cmag[gdmag,gdf[k]])**2 
            scatter[:,i] = np.sqrt( totaldiff/np.maximum(ndet[:,i],1) ) 
            bd, = np.where(ndet[:,i] <= 1) 
            if len(bd) > 0 : 
                scatter[bd,i]=99.99 
    # Create final catalog schema
    photdt = [('objid',(np.str,100)),('x',float),('y',float),('ra',float),('dec',float)]

    # Add columns for calibrated single-epoch photometry columns 
    cmagnames = np.zeros(nchstr,(np.str,100))
    cerrnames = np.zeros(nchstr,(np.str,100))        
    for i in range(nufilt): 
        ind, = np.where(chstr['filter'] == ufilt[i]) 
        cmagnames[ind] = [ufilt[i].upper()+str(j+1)+'MAG' for j in np.arange(nind)]
        cerrnames[ind] = [ufilt[i].upper()+str(j+1)+'ERR' for j in np.arange(nind)]            
    for i in range(nchstr):
        photdt += [(cmagnames[i],np.float32),(cerrnames[i],np.float32)]
    # Add columns for average photometry per filter 
    for i in range(nufilt):
        photdt += [(ufilt[i]+'MAG',np.float32),(ufilt[i]+'ERR',np.float32),(ufilt[i]+'SCATTER',np.float32),('NDET'+ufilt[i],int)]
    # Extra columns
    photdt += [('chi',np.float32),('sharp',np.float32),('prob',np.float32),('ebv',np.float32)]
    # other SE columns
    photdt += [('mag_auto',np.float32),('magerr_auto',np.float32),('asemi',np.float32),('bsemi',np.float32),
               ('theta',np.float32),('ellipticity',np.float32),('fwhm',np.float32)]
    # in unique brick area
    dt += [('brickuniq',boolean)]
    # Create final catalog
    phot = np.zeros(ninstphot,dtype=np.dtype(photdt))
    for n in instphot:
        phot[n] = instphot[n]
    phtags = phot.colnames
    # object IDs 
    phot['objid'] = brick+'.'+np.char.array(instphot['id'])
    # Stuff in the coordinates calculated above 
    phot['ra'] = ra 
    phot['dec'] = dec 
    # Stuff in the calibrated single-epoch photometry columns 
    for i in range(nchstr): 
        phot[cmagnames[i]] = cmag[:,i] 
        phot[cerrnames[i]] = cerr[:,i] 
    # Stuff in the average photometry per filter 
    for i in range(nufilt): 
        phot[ufilt[i].upper()+'MAG'] = avgmag[:,i] 
        phot[ufilt[i].upper()+'ERR'] = avgerr[:,i] 
        phot[ufilt[i].upper()+'SCATTER'] = scatter[:,i] 
        phot['NDET'+ufilt[i].upper()] = ndet[:,i] 
         
    # Calculate SFD E(B-V) 
    #GLACTC,phot.ra,phot.dec,2000.0,glon,glat,1,/deg
    coo = SkyCoord(phot['ra'], phot['dec'], unit='deg', frame='icrs')
    coo_gal = coo.transform_to('galactic')
    sfd = SFDQuery()
    phot['ebv'] = sfd(coords_gal)
         
    # THIS IS NOW BEING DONE IN DELVERED_FINALCAT.PRO THAT COMBINES ALL CATALOGS 
    # Only include objects that are INSIDE the UNIQUE brick area 
    #   INCLUSIVE at the lower RA and DEC limit 
    # Getting objects that are in the UNIQUE brick area 
    if brickstr1.dec == -90: 
        # the brick right at the pole does not have any RA limits 
        ginside, = np.where(phot['dec'] < brickstr1['dec2']) 
    else: 
        ginside, = np.where((phot['ra'] >= brickstr1['ra1']) & (phot['ra'] < brickstr1['ra2']) &
                            (phot['dec'] >= brickstr1['dec1']) & (phot['dec'] < brickstr1['dec2']))
    if len(ginside)>0: 
        phot['brickuniq'][ginside] = True
         
    # Get some meta-data 
    for i in range(nchstr): 
        alffile = chstr['base'][i]+'.alf' 
        if os.path.exists(alffile):
            chstr['alf_nsources'][i] = dln.numlines(alffile)-3 
         
         
    # Make the exposure-level forced photometry catalog 
    #-------------------------------------------------- 
    # Load the individual ALF files to get chi, sharp 
    # Load TFR file to conver ALF IDs to final object ID 
    tfrfile = mchbase+'_comb.tfr' 
    alffiles,tfrstr = io.readtfr(tfrfile)
    expdt = [('id',(np.str,100)),('objid',(np.str,100)),('exposure',(np.str,100)),('ccdnum',int),('filter',(np.str,10)),
             ('mjd',float),('x',float),('y',float),('ra',float),('dec',float),('imag',float),('ierr',float),
             ('mag',float),('err',float),('sky',float),('chi',float),('sharp',float)]
    expcat = np.zeros(int(np.sum(cmag < 50))+10000,dtype=np.dtype(expdt))
    cnt = 0
    for i in range(nchstr): 
        base1 = chstr['base'][i]
        fitsfile = base1+'.fits' 
        if os.path.exists(alffiles[i]) and os.path.exists(fitsfile):
            alf = io.readals(alffiles[i])
            if len(alf)==0: 
                continue
            # Sometimes the rows are duplicated in the ALF file 
            ui = np.uniq(alf.id,np.argsort(alf.id)) 
            if len(ui) < nalf: 
                alf = alf[ui] 
                nalf = len(alf) 
            head = io.readfile(fitsfile,header=True) 
                 
            # Calibrate the photometry 
            # exptime, aperture correction, zero-point 
            # aperture correction is SUBTRACTIVE, makes it brighter 
            # ZPTERM is a SUBTRACTIVE constant offset 
            cmag1 = alf['mag'] + 2.5*np.log10(chstr['exptime'][i]) - chstr['apcor'][i] - chstr['calib_zpterm'][i]
            # Add zero-point error in quadrature 
            cerr1 = np.sqrt(alf['err']**2+chstr['calib_zptermsig'][i]**2) 
                 
            # Coordinates
            ra1,dec1 = head.pixel_to_world(alf['x']-1,alf['y']-1)
            #HEAD_XYAD,head,alf.x-1,alf.y-1,ra1,dec1,/deg 
                 
            # MATCH up ALF IDs to TFR INDEX 
            # One row per unique object 
            # the INDEX values are 1-based indices into the ALF files
            ind1,ind2 = dln.match(tfrstr.index[i],lindgen(nalf)+1)
            #MATCH,tfrstr.index[i],lindgen(nalf)+1,ind1,ind2,/sort,count=nmatch 
            objid = brick+'.'+np.char.array(tfrstr['id'][ind1])
                 
            # Create the new catalog
            newcat = np.zeros(nalf,dtype=np.dtype(dt))
            newcat['objid'] = objid 
            newcat['id'] = chstr['expnum'][i]+'_'+chstr['chip'][i]+'.'+np.char.array(alf['id'])
            newcat['exposure'] = chstr['base'][i]
            newcat['ccdnum'] = chstr['chip'][i]
            newcat['filter'] = chstr['filter'][i]
            dateobs = chstr['utdate'][i]+'T'+chstr['uttime'][i]
            t = Time(dateobs)
            newcat['mjd'] = t.mjd
            newcat['x'] = alf['x']
            newcat['y'] = alf['y']
            newcat['ra'] = ra1 
            newcat['dec'] = dec1 
            newcat['imag'] = alf['mag']
            newcat['ierr'] = alf['err']
            newcat['mag'] = cmag1 
            newcat['err'] = cerr1 
            newcat['sky'] = alf['sky']
            newcat['chi'] = alf['chi']
            newcat['sharp'] = alf['sharp']
                 
            # Add more elements 
            if cnt+nalf > len(expcat): 
                expcat = add_elements(expcat,100000>nalf) 
                 
            # Add to global catalog 
            expcat[cnt:cnt+nalf-1] = newcat 
            cnt += nalf
                
        else:
            logger.info(alffile+' NOT FOUND')
        expcat = expcat[0:cnt]   # this should not be needed 
             
        # Object catalog 
        #---------------
        # do not want all of the individual epoch photometry
        #  these are between the DEC and the first unique mean magnitude (e.g. GMAG, and GERR)
        lo = where(phtags=='DEC')[0][0]
        hi = where(np.char.array(phtags).upper() == ufilt[0].upper()+'MAG')[0][0]
        objdt = photdt[0:lo]
        objdt += photdt[hi:]
        obj = np.zeros(ninstphot,dtype=np.dtype(objdt))
        for n in obj.colnames:
            obj[n] = phot[h]
            
        # Saving final catalog 
        logger.info(datetime.now().strftime("%a %b %d %H:%M:%S %Y"))                         
        photfile = bdir+brick+'.fits'
        if len(phot.colnames) <= 999:
            logger.info('Writing photometry to '+photfile+'.gz')
            phot.writeto(photfile,overwrite=True)
            out = subprocess.check_output(['gzip','-f',photfile],shell=False)
        else: 
            logger.info('Too many columns for FITS.  Saving as pickle file instead. '+bdir+brick+'.pkl')
            with open(bdir+brick+'.pkl','wb') as f:
                pickle.dump(phot,f)
             
        # Saving object photometry catalog 
        objfile = bdir+brick+'_object.fits'
        obj.writeto(objfile,overwrite=True)
        out = subprocess.check_output(['gzip','-f',objfile],shell=False)                            
             
        # Saving exposure-level forced photometry catalog 
        expfile = bdir+brick+'_expforced.fits'
        expcat.writeto(expfile,overwrite=True)
        out = subprocess.check_output(['gzip','-f',expfile],shell=False)
             
        # Save metadata 
        metafile = bdir+brick+'_meta.fits' 
        logger.info('Writing meta-data to '+metafile)
        chstr.writeto(metafile,overwrite=True)
             
             
        # Delete files we don't want to keep 
        #----------------------------------- 
        # Individual fits files 
        alsbase = os.path.basename(alsfiles,'.als') 
        # if fits files have resources files then replace the fits file by a 
        #  dummy file 
        for i in range(len(alsbase)):
            exists = os.path.exists(alsbase[i]+'.fits')
            size = os.path.getsize(alsbase[i]+'.fits')
            if exists and size > 1: 
                os.remove(alsbase[i]+'.fits')
                if os.path.exists('.'+alsbase[i]+'.fits'):
                    dln.writelines(alsbase[i]+'.fits','')
        # Combined files 
        #   _comb  lst, lst1, lst2, lst1.chi, grp, nst, lst2.chi, plst.chi, psfini.ap 
        #   nei, als.inp, a.fits, cmn.log, cmn.coo, cmn.ap, cmn.lst, 
        #   _sub.fits, _sub.cat, _sub.als, _all.coo, makemag 
        base = os.path.basename(mchfile,'.mch')
        ext = ['.lst','.lst1','.lst2','.lst1.chi','.lst2.chi','.grp','.nst','.plst.chi','.nei',
               '.als.inp','.cmn.log','.cmn.coo','.cmn.ap','.cmn.lst','a.fits','a.fits.fz',
               '_sub.fits','_sub.cat','_sub.als','_all.coo','.makemag']
        for e in ext:
            if os.path.exists(base+e): os.remove(base+e)
        if os.path.exists('check.fits'): os.remove('check.fits')
             
        # fpack _comb.fits and _combs.fits
        for f in [base+'_comb.fits.fz',base+'_combs.fits.fz']:
            if os.path.exists(f): os.remove(f)
        out = subprocess.check_output(['fpack','-D','-Y',base+'_comb.fits'],shell=False)
        out = subprocess.check_output(['fpack','-D','-Y',base+'_combs.fits'],shell=False)
        # gzip _comb.mask.fits and _comb.bpm.fits 
        out = subprocess.check_output(['gzip','-f',base+'_comb.mask.fits'],shell=False)
        out = subprocess.check_output(['gzip','-f',base+'_comb.bpm.fits'],shell=False)
                      
        # Copy everything left back to the original directory 
        logger.info('Copying files back to '+bdir)
        allfiles = glob(procdir+'*')
        allfiles += glob(procdir+'.*.fits')
        nallfiles = len(allfiles)
        for f in allfiles:
            base = os.path.basename(f)
            if os.path.exists(bdir+'/'+base): os.remove(bdir+'/'+base)
            shutil.move(f,bdir)
        os.rmdir(procdir)  # delete temporary processing directory 
             
        logger.info(datetime.now().strftime("%a %b %d %H:%M:%S %Y"))                       
        logger.info('DELVERED_FORCEBRICK:ne after '+str(time.time()-t0,2)+' sec.')

        os.chdir(curdir)  # back to original directory 
             
        # Create JOINT catalogs 
        logger.info('')
        logger.info('CREATE JOINT CATALOGS')
        logger.info('')
        logger.info(datetime.now().strftime("%a %b %d %H:%M:%S %Y"))                       
             
        if redo or update:
            jntredo = True 
        else: 
            jntredo = False
        delvered_jointbrickcats(brick,logfile=logfile,redo=jntredo)
             
 
        
