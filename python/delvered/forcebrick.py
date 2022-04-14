#!/usr/bin/env python

import os
import time
import numpy as np
import socket
import logging
import tempfile
from datetime import datetime
from dlnpyutils import utils as dln

 
def forcebrick(brick,scriptsdir=None,irafdir=None,workdir=None,redo=False,
               update=False,logfile=None,delvedir=None):
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
        delvereddir = '/home/dnidever/projects/delvered/'
    if scriptsdir is not None:
        if scriptsdir.endswith('/')==False:
            scriptsdir += '/'
    else: 
        scriptsdir = '/home/dnidever/projects/PHOTRED/scripts/'
    if irafdir is not None:
        if irafdir.endswith('/')==False:
            irafdir += '/'
    else: 
        irafdir='/home/dnidever/iraf/' 
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
    brickdir = trailingslash(delvedir)+'bricks/' 
    if os.path.exists(brickdir)==False:
        os.makedirs(brickdir)
    logsdir = brickdir+'logs/' 
    if os.path.exists(logsdir)==False:
        os.makedirs(logsdir)
    if logfile is None:
        logfile = -1 
     
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
    rootLogger.addHandler(fileHandler)
    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    logger.addHandler(consoleHandler)
    logger.setLevel(logging.NOTSET)
     
    # Chip name and CCDNUM relations 
    decam = IMPORTASCII(delvereddir+'data/decam.txt',/header,/silent) 
     
    # Load the brick information 
    brickstr = MRDFITS(delvereddir+'data/delvemc_bricks_0.25deg.fits.gz',1,/silent) 
     
    # Get the brick information 
    bind , = np.where(brickstr.brickname == brick,nbind) 
    if nbind == 0: 
        logger.info(ibrick+' not in DELVE-MC brick list' 
        return 
    brickstr1 = brickstr[bind[0]] 
     
    # Subdirectory is the first 4 digits of the brickname, e.g., 0952 of 0952m462, the RA portion of the name 
    subdir = brickdir+strmid(brick,0,4)+'/' 
    if os.path.exists(subdir,/directory) == 0 : 
        file_mkdir,subdir 
    # Create brick directory if it doesn't exist yet 
    bdir = subdir+brick+'/' 
    if os.path.exists(bdir,/directory) == 0 : 
        file_mkdir,bdir 
    logfile = bdir+brick+'.'+logtime+'.log' 
     
    # Check the output file 
    photfile = bdir+brick+'.fits'
    if os.path.exists(photfile+'.gz') and redo==False and update==False:
        logger.info(photfile+'.gz EXISTS and redo or update NOT set')
        return 
     
    # DECam imager 
    thisimager = {telescope:'BLANCO',instrument:'DECam',namps:62,separator:'_'} 
     
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
    setup = photred.utils.loadsetup(setupdir=bdir)
     
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
    tid,tmpfile = tempfile.mkstemp(prefix="tmp",dir=tempdir)
    dln.touchzero(tmpfile+'.fits')
    for f in [tmpfile,tmpfile+'.fits']:
        if os.path.exists(f): os.remove(f)
    tmpfile += '.fits' 
    # /noshell causes problems on gp09 because it gives python2 instead of python3 
    spawn,delvereddir+'bin/query_delvered_summary_table '+str(cenra,2)+' '+str(cendec,2)+' '+tmpfile+' --lim 0.5',out,errout 
    info = file_info(tmpfile) 
    if info.size == 0: 
        logger.info('No overlapping chips found')
        return
    chstr = Table.read(tmpfile)
    if os.path.exists(tmpfile): os.remove(tmpfile)
    nchstr = len(chstr) 
    logger.info('Found '+str(nchstr)+' overlapping chips within 0.5 deg of brick center')
    logger.info(datetime.now().strftime("%a %b %d %H:%M:%S %Y")) 
     
    # Make sure the chips are unique,  some were duplicated on SMASH nights 
    chid = chstr['expnum']+'-'+str(chstr['chip'])
    uchstr,ui = np.unique(chid,return_index=True)
    chstr = chstr[ui] 
    nchstr = len(chstr) 
     
    # Do more rigorous overlap checking 
    #  the brick region with overlap 
    print('Performing more rigorous overlap checking')
    HEAD_XYAD,tilestr.head,[0,tilestr.nx-1,tilestr.nx-1,0],[0,0,tilestr.ny-1,tilestr.ny-1],bvra,bvdec,/deg 
    olap = intarr(nchstr) 
    vxarr = fltarr(nchstr,4) 
    vyarr = fltarr(nchstr,4) 
    for i in range(nchstr): 
        if (i mod 100 == 0) and (i > 0): 
            print(i)
        hd1 = PHOTRED_READFILE(chstr[i].file,/header) 
        nx = sxpar(hd1,'naxis1') 
        ny = sxpar(hd1,'naxis2') 
        head_xyad,hd1,[0,nx-1,nx-1,0],[0,0,ny-1,ny-1],vra,vdec,/degree 
        olap[i] =:polygonsoverlap(bvra,bvdec,vra,vdec) 
        head_adxy,tilestr.head,vra,vdec,vx,vy,/deg 
        vxarr[i,:] = vx 
        vyarr[i,:] = vy 
    # Require at least a 2 pixel overlap in X and Y 
    g , = np.where(olap == 1 and max(vxarr,dim=2) >= 2 and max(vyarr,dim=2) >= 2 and min(vxarr,dim=2) <= tilestr.nx-3 and min(vyarr,dim=2) <= tilestr.ny-3,ng) 
    if ng == 0: 
        logger.info('No chips overlap this brick')
        return 
    logger.info(str(ng)+' chips overlap this brick')
    chstr = chstr[g] 
    nchstr = ng 
     
    # APPLY cuts on the exposures 
    #----------------------------- 
    logger.info('Applying quality, zero-point, filter and exptime cuts')
     
    # Zero-point structure, from NSC 
    zpstr = replicate({instrument:'',filter:'',amcoef:fltarr(2),thresh:0.5},7) 
    zpstr.instrument = 'c4d' 
    zpstr.filter = ['u','g','r','i','z','Y','VR'] 
    zpstr[0].amcoef = [-1.60273, -0.375253]  # c4d-u 
    zpstr[1].amcoef = [0.277124, -0.198037]  # c4d-g 
    zpstr[2].amcoef = [0.516382, -0.115443]  # c4d-r 
    zpstr[3].amcoef = [0.380338, -0.067439]  # c4d-i 
    zpstr[4].amcoef = [0.074517, -0.067031]  # c4d-z 
    zpstr[5].amcoef = [-1.07800, -0.060014]  # c4d-Y 
    zpstr[6].amcoef = [1.111859, -0.083630]  # c4d-VR 
     
    # Convert to additive zero-point as used in NSC 
    zpterm = -chstr.calib_zpterm 
    # Fix early DES exposures that used different units/gain 
    gdes, = np.where(chstr.gain < 2,ngdes) 
    if len(gdes) > 0: 
        zpterm[gdes] -= 1.55 
    # global zpterm and airmass correction 
    for i in range(len(zpstr)): 
        ind , = np.where(chstr.filter == zpstr[i].filter,nind) 
        if nind > 0 : 
            zpterm[ind] -= poly(chstr[ind].airmass,zpstr[i].amcoef) 
     
    fwhmthresh = 2.0# seeing 2.0" threshold 
    filt = strmid(chstr.filter,0,1) 
    gdch , = np.where(chstr.fwhm*chstr.pixscale <= fwhmthresh and chstr.exptime >= 90. and zpterm >= -0.5 and              finite(zpterm) == 1 and finite(chstr.apcor) == 1 and              (filt == 'u' or filt == 'g' or filt == 'r' or filt == 'i' or filt == 'z' or filt == 'Y'),ngdch) 
     
    if ngdch == 0: 
        logger.info('No chips passed the cuts')
        return 
    logger.info(str(ngdch)+' chips passed the cuts')
    chstr = chstr[gdch] 
    nchstr = ngdch 
    chstr.file = str(chstr.file,2) 
    chstr.base = str(chstr.base,2) 
     
     
    # Check if we need to update
    if update and redo==False:
        # Load previous meta 
        metafile = bdir+brick+'_meta.fits' 
        if os.path.exists(metafile): 
            meta0 = mrdfits(metafile,1,/silent) 
            meta0.base = str(meta0.base,2) 
            MATCH,meta0.base,chstr.base,ind1,ind2,/sort,count=nmatch 
            if nchstr == len(meta0) and nchstr == nmatch: 
                logger.info('Nothing to UPDATE') 
                return 
            # Moving previous catalogs to a backup 
            oldfiles = file_search(bdir+[brick+'*fits*','allframe.opt','default.*'],count=noldfiles) 
            if noldfiles > 0: 
                bakdir = bdir+'bak'+logtime
                logger.info('Backing up old files to '+bakdir)
                FILE_MKDIR,bakdir 
                FILE_MOVE,oldfiles,bakdir 
             else logger.info('No old files to backup')
            logger.info('There are new exposures to include.  UPDATING')
         
         
         
        # Create temporary local directory to perform the work/processing 
        #  copy everything to it 
        if os.path.exists(workdir,/directory) == 0 : 
            FILE_MKDIR,workdir 
        procdir = first_el(MKTEMP('brk',outdir=workdir,/directory))+'/' 
        procdir = repstr(procdir,'//','/') 
        FILE_CHMOD,procdir,/a_execute 
        logger.info('Working in temporary directory '+procdir)
        logger.info(datetime.now().strftime("%a %b %d %H:%M:%S %Y")) 
                    
        # Copy the files 
        #---------------- 
        logger.info('Copying over the necessary files to '+procdir)
        # Chip loop 
        for i in range(nchstr): 
            logger.info('Copying over files for ',chstr[i].file 
            # New filename 
            #chstr[i].newbase = chstr[i].base 
            odir = os.path.dirname(chstr[i].file)+'/' 
            obase = PHOTRED_GETFITSEXT(chstr[i].file,/basename) 
            obase = first_el(obase) 
            # Need symlinks to .psf, .als 
            if os.path.exists(odir+obase+'.psf') == 0: 
                logger.info(odir+obase+'.psf NOT FOUND.  Skipping this chip' 
                goto,BOMB1 
            if os.path.exists(odir+obase+'.als') == 0: 
                logger.info(odir+obase+'.als NOT FOUND.  Skipping this chip' 
                goto,BOMB1 
            os.remove(procdir+chstr[i].base+['.psf','.als','.ap','.opt','.als.opt','.log'],/allow 
            FILE_COPY,odir+obase+['.psf','.als','.ap','.opt','.als.opt','.log'],procdir,/allow_same,/overwrite 
            FILE_CHMOD,procdir+chstr[i].base+['.psf','.als','.ap','.opt','.als.opt','.log'],'755'o# make sure they are writable 
            # Copy the fits, fits resource file and header files locally 
            if os.path.exists(odir+obase+'.fits') == 0: 
                logger.info(odir+obase+'.fits NOT FOUND.  Skipping this chip' 
                goto,BOMB1 
            FILE_COPY,odir+obase+'.fits',procdir,/over 
            if os.path.exists(odir+'.'+obase+'.fits') == 1 : 
                FILE_COPY,odir+'.'+obase+'.fits',procdir,/over 
            if os.path.exists(odir+obase+'.fits.head') == 1 : 
                FILE_COPY,odir+obase+'.fits.head',procdir,/over 
            FILE_CHMOD,procdir+chstr[i].base+'.fits','755'o 
            BOMB1: 
         
         
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
        ALLFRAME,mchfile,tile=tilestr,setupdir=bdir,scriptsdir=scriptsdir,irafdir=irafdir,         logfile=logfile,catformat='FITS',imager=thisimager,geocoef=0 
        magfile = chstr[0].base+'.mag' 
        if os.path.exists(magfile) == 0: 
            logger.info(magfile+' NOT FOUND')
            return 
         
        # Step 4: Calculate coordinates 
        #------------------------------- 
        logger.info('Step 4: Adding coordinates')
        logger.info(datetime.now().strftime("%a %b %d %H:%M:%S %Y"))                     
        # Load the MCH file
        alsfiles,trans,magoff = photred.utils.loadmch(mchfile)
        nalsfiles = len(alsfiles) 
        # Load the photometry file 
        instphot = PHOTRED_READFILE(magfile) 
        ninstphot = len(instphot) 
        logger.info('Nstars = '+str(ninstphot)) 
        # Converting to IDL X/Y convention, starting at (0,0) 
        # DAOPHOT has X/Y start at (1,1) 
        HEAD_XYAD,tilestr.head,instphot.x-1.0,instphot.y-1.0,ra,dec,/degree 
         
         
        # Step 5: Calibrating photometry with zero-points 
        #------------------------------------------------- 
        logger.info('Step 5: Calibrating photometry with zero-points')
        logger.info(datetime.now().strftime("%a %b %d %H:%M:%S %Y"))                     
        cmag = fltarr(ninstphot,nchstr)+99.99 
        cerr = fltarr(ninstphot,nchstr)+9.99 
        # Chip loop 
        for i in range(nchstr): 
            # id, x, y, ra, dec, unsolved magnitudes, chi, sharp 
            imag = instphot.(2*i+3) 
            ierr = instphot.(2*i+4) 
            gdmag , = np.where(imag < 50,ngdmag) 
            if ngdmag > 0: 
                # exptime, aperture correction, zero-point 
                # aperture correction is SUBTRACTIVE, makes it brighter 
                # ZPTERM is a SUBTRACTIVE constant offset 
                cmag[gdmag,i] = imag[gdmag] + 2.5*alog10(chstr[i].exptime) - chstr[i].apcor - chstr[i].calib_zpterm 
                # Add zero-point error in quadrature 
                cerr[gdmag,i] = sqrt(ierr[gdmag]**2+chstr[i].calib_zptermsig**2) 
        # Calculate average photometry per filter 
        ufilt = chstr[np.uniq(chstr.filter,np.argsort(chstr.filter))].filter 
        nufilt = len(ufilt) 
        avgmag = fltarr(ninstphot,nufilt) 
        avgerr = fltarr(ninstphot,nufilt) 
        ndet = lonarr(ninstphot,nufilt) 
        for i in range(nufilt): 
            gdf , = np.where(chstr.filter == ufilt[i],ngdf) 
            # Single exposure 
            if ngdf == 1: 
                avgmag[:,i] = cmag[:,gdf[0]] 
                avgerr[:,i] = cerr[:,gdf[0]] 
                ndet[:,i] = (cmag[:,gdf[0]] < 50) 
                # Multiple exposures 
            else: 
                # Loop through all of the exposures and add up the flux, totalwt, etc. 
                totalwt = dblarr(ninstphot) 
                totalfluxwt = dblarr(ninstphot) 
                for k in range(ngdf): 
                    gdmag , = np.where(cmag[:,gdf[k]] < 50,ngdmag) 
                    if ngdmag > 0: 
                        totalwt[gdmag] += 1.0d0/cerr[gdmag,gdf[k]]**2 
                        totalfluxwt[gdmag] += 2.5118864d**cmag[gdmag,gdf[k]] * (1.0d0/cerr[gdmag,gdf[k]]**2) 
                        ndet[gdmag,i]++ 
                newflux = totalfluxwt/totalwt 
                newmag = 2.50*alog10(newflux) 
                newerr = sqrt(1.0/totalwt) 
                bdmag , = np.where(finite(newmag) == 0,nbdmag) 
                if nbdmag > 0: 
                    newmag[bdmag] = 99.99 
                    newerr[bdmag] = 9.99 
                avgmag[:,i] = newmag 
                avgerr[:,i] = newerr 
        # Measure scatter 
        scatter = fltarr(ninstphot,nufilt)+99.99 
        for i in range(nufilt): 
            gdf , = np.where(chstr.filter == ufilt[i],ngdf) 
            if ngdf > 1: 
                totaldiff = fltarr(ninstphot) 
                for k in range(ngdf): 
                    gdmag , = np.where(cmag[:,gdf[k]] < 50,ngdmag) 
                    if ngdmag > 0 : 
                        totaldiff[gdmag] += (avgmag[gdmag,i]-cmag[gdmag,gdf[k]])**2 
                scatter[:,i] = sqrt( totaldiff/(ndet[:,i]>1) ) 
                bd , = np.where(ndet[:,i] <= 1,nbd) 
                if nbd > 0 : 
                    scatter[bd,i]=99.99 
        # Create final catalog schema 
        newschema = {objid:'',x:0.0,y:0.0,ra:0.0d0,dec:0.0d0} 
        # Add columns for calibrated single-epoch photometry columns 
        cmagnames = strarr(nchstr) 
        cerrnames = strarr(nchstr) 
        for i in range(nufilt): 
            ind , = np.where(chstr.filter == ufilt[i],nind) 
            cmagnames[ind] = strupcase(ufilt[i])+str(lindgen(nind)+1,2)+'MAG' 
            cerrnames[ind] = strupcase(ufilt[i])+str(lindgen(nind)+1,2)+'ERR' 
        for i in range(nchstr): 
            newschema = create_struct(newschema,cmagnames[i],0.0,cerrnames[i],0.0) 
        # Add columns for average photometry per filter 
        for i in range(nufilt): 
            newschema = create_struct(newschema,ufilt[i]+'MAG',0.0,ufilt[i]+'ERR',0.0,ufilt[i]+'SCATTER',0.0,'NDET'+ufilt[i],0L) 
        # Extra columns 
        newschema = create_struct(newschema,'chi',0.0,'sharp',0.0,'prob',0.0,'ebv',0.0) 
        # other SE columns 
        newschema = create_struct(newschema,'mag_auto',0.0,'magerr_auto',0.0,'asemi',0.0,'bsemi',0.0,'theta',0.0,'ellipticity',0.0,'fwhm',0.0) 
        # in unique brick area 
        newschema = create_struct(newschema,'brickuniq',0B) 
        # Create final catalog 
        phot = replicate(newschema,ninstphot) 
        struct_assign,instphot,phot,/nozero 
        phtags = tag_names(phot) 
        # object IDs 
        phot.objid = brick+'.'+str(instphot.id,2) 
        # Stuff in the coordinates calculated above 
        phot.ra = ra 
        phot.dec = dec 
        # Stuff in the calibrated single-epoch photometry columns 
        for i in range(nchstr): 
            magind , = np.where(strupcase(phtags) == cmagnames[i],nmagind) 
            errind , = np.where(strupcase(phtags) == cerrnames[i],nerrind) 
            phot.(magind) = cmag[:,i] 
            phot.(errind) = cerr[:,i] 
        # Stuff in the average photometry per filter 
        for i in range(nufilt): 
            magind , = np.where(strupcase(phtags) == strupcase(ufilt[i])+'MAG',nmagind) 
            errind , = np.where(strupcase(phtags) == strupcase(ufilt[i])+'ERR',nerrind) 
            scatind , = np.where(strupcase(phtags) == strupcase(ufilt[i])+'SCATTER',nscatind) 
            detind , = np.where(strupcase(phtags) == 'NDET'+strupcase(ufilt[i]),ndetind) 
            phot.(magind) = avgmag[:,i] 
            phot.(errind) = avgerr[:,i] 
            phot.(scatind) = scatter[:,i] 
            phot.(detind) = ndet[:,i] 
         
        # Calculate SFD E(B-V) 
        GLACTC,phot.ra,phot.dec,2000.0,glon,glat,1,/deg 
        phot.ebv = dust_getval(glon,glat,/noloop,/interp) 
         
        # THIS IS NOW BEING DONE IN DELVERED_FINALCAT.PRO THAT COMBINES ALL CATALOGS 
        # Only include objects that are INSIDE the UNIQUE brick area 
        #   INCLUSIVE at the lower RA and DEC limit 
        # Getting objects that are in the UNIQUE brick area 
        if brickstr1.dec == -90: 
            # the brick right at the pole does not have any RA limits 
            ginside , = np.where(phot.dec < brickstr1.dec2,ninside) 
        else: 
            ginside , = np.where(phot.ra >= brickstr1.ra1 and phot.ra < brickstr1.ra2 and                   phot.dec >= brickstr1.dec1 and phot.dec < brickstr1.dec2,ninside) 
        if ninside > 0 : 
            phot[ginside].brickuniq=1B 
         
        # Get some meta-data 
        for i in range(nchstr): 
            alffile = chstr[i].base+'.alf' 
            if os.path.exists(alffile) == 1 : 
                chstr[i].alf_nsources=file_lines(alffile)-3 
         
         
        # Make the exposure-level forced photometry catalog 
        #-------------------------------------------------- 
        # Load the individual ALF files to get chi, sharp 
        # Load TFR file to conver ALF IDs to final object ID 
        tfrfile = mchbase+'_comb.tfr' 
        LOADTFR,tfrfile,alffiles,tfrstr 
        schema = {id:'',objid:'',exposure:'',ccdnum:0,filter:'',mjd:0.0d0,x:0.0,y:0.0,ra:0.0d0,dec:0.0d0,          imag:0.0,ierr:0.0,mag:0.0,err:0.0,sky:0.0,chi:0.0,sharp:0.0} 
        expcat = replicate(schema,int(np.sum(cmag < 50))+10000L) 
        cnt = 0LL 
        for i in range(nchstr): 
            base1 = chstr[i].base 
            fitsfile = base1+'.fits' 
            if os.path.exists(alffiles[i]) == 1 and os.path.exists(fitsfile) == 1: 
                LOADALS,alffiles[i],alf,count=nalf 
                if nalf == 0 : 
                    goto,BOMB2 
                # Sometimes the rows are duplicated in the ALF file 
                ui = np.uniq(alf.id,np.argsort(alf.id)) 
                if len(ui) < nalf: 
                    alf = alf[ui] 
                    nalf = len(alf) 
                head = photred_readfile(fitsfile,/header) 
                 
                # Calibrate the photometry 
                # exptime, aperture correction, zero-point 
                # aperture correction is SUBTRACTIVE, makes it brighter 
                # ZPTERM is a SUBTRACTIVE constant offset 
                cmag1 = alf.mag + 2.5*alog10(chstr[i].exptime) - chstr[i].apcor - chstr[i].calib_zpterm 
                # Add zero-point error in quadrature 
                cerr1 = sqrt(alf.err**2+chstr[i].calib_zptermsig**2) 
                 
                # Coordinates 
                HEAD_XYAD,head,alf.x-1,alf.y-1,ra1,dec1,/deg 
                 
                # MATCH up ALF IDs to TFR INDEX 
                # One row per unique object 
                # the INDEX values are 1-based indices into the ALF files 
                MATCH,tfrstr.index[i],lindgen(nalf)+1,ind1,ind2,/sort,count=nmatch 
                objid = brick+'.'+str(tfrstr[ind1].id,2) 
                 
                # Create the new catalog 
                newcat = replicate(schema,nalf) 
                newcat['objid'] = objid 
                newcat['id'] = chstr['expnum'][i]+'_'+str(chstr['chip'][i])+'['+str(alf['id'])
                newcat['exposure'] = chstr['base'][i]
                newcat['ccdnum'] = chstr['chip'][i]
                newcat['filter'] = chstr['filter'][i] 
                newcat['mjd'] = date2jd(chstr['utdate'][i]+'T'+chstr['uttime'][i],mjd=True)
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
                if cnt+nalf > len(expcat) : 
                    expcat=add_elements(expcat,100000L>nalf) 
                 
                # Add to global catalog 
                expcat[cnt:cnt+nalf-1] = newcat 
                cnt += nalf 
                 
                BOMB2: 
             else:
                    logger.info(alffile+' NOT FOUND')
            expcat = expcat[0:cnt-1]# this should not be needed 
             
            # Object catalog 
            #--------------- 
            lo = first_el(where(phtags == 'DEC',nlo)) 
            hi = first_el(where(strupcase(phtags) == strupcase(ufilt[0])+'MAG',nhi)) 
            obj_schema = create_struct(phtags[0],fix('',type=size(phot.(0),/type))) 
            for i in np.arange(1,lo+1): 
                obj_schema = create_struct(obj_schema,phtags[i],fix('',type=size(phot.(i),/type))) 
            for i in np.arange(hi,len(phtags)-2+1): 
                obj_schema = create_struct(obj_schema,phtags[i],fix('',type=size(phot.(i),/type))) 
            obj_schema = create_struct(obj_schema,'brickuniq',0B) 
            obj = replicate(obj_schema,ninstphot) 
            STRUCT_ASSIGN,phot,obj,/nozero 
             
            # Saving final catalog 
            logger.info(datetime.now().strftime("%a %b %d %H:%M:%S %Y"))                         
            photfile = bdir+brick+'.fits' 
            if n_tags(phot) <= 999: 
                logger.info('Writing photometry to '+photfile+'.gz')
                MWRFITS,phot,photfile,/create 
                spawn,['gzip','-f',photfile],/noshell 
            else: 
                logger.info('Too many columns for FITS.  Saving as IDL SAVE file instead. '+bdir+brick+'.dat') 
                SAVE,phot,file=bdir+brick+'.dat' 
             
            # Saving object photometry catalog 
            objfile = bdir+brick+'_object.fits'
            obj.writeto(objfile,overwrite=True)
            out = subprocess.check_output(['gzip','-f',objfile],shell=False)                            
             
            # Saving exposure-level forced photometry catalog 
            expfile = bdir+brick+'_expforced.fits' 
            MWRFITS,expcat,expfile,/create 
            out = subprocess.check_output(['gzip','-f',expfile],shell=False)
             
            # Save metadata 
            metafile = bdir+brick+'_meta.fits' 
            logger.info('Writing meta-data to '+metafile)
            MWRFITS,chstr,metafile,/create 
             
             
            # Delete files we don't want to keep 
            #----------------------------------- 
            # Individual fits files 
            alsbase = os.path.basename(alsfiles,'.als') 
            # if fits files have resources files then replace the fits file by a 
            #  dummy file 
            for i in range(len(alsbase)):
                exists = os.path.exists(alsbase[i]+'.fits')
                size = os.size(alsbase[i]+'.fits')
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
            if os.path.exists('check.fits): os.remove('check.fits')
             
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
            allfiles = file_search(procdir+['*','.*.fits'],count=nallfiles) 
            FILE_MOVE,allfiles,bdir,/allow,/overwrite 
            os.remove(procdir)  # delete temporary processing directory 
             
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
             
 
        
