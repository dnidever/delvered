#!/usr/bin/env python

import os
import time
import numpy as np
import socket
from dlnpyutils import utils as dln
from . import forcebrick

def bricks(input,nmulti=10,redo=False,update=False,delvedir=None,delvereddir=None):
    """
    This is a script to run PHOTRED ALLFRAME on DELVE MC exposures to create 
    forced-photometry catalogs. 
 
    Parameters
    ----------
    input : str or list
      What bricks to run PHOTRED/ALLFRAME on.  Either an array or 
        a range such as 100-110. 
    nmulti : int, optional
       Number of simultaneous jobs to run.  Default is 10.
    redo : boolean, optional
       Redo bricks that were previously done.  Default is False.
    update : boolean, optional
       Check for updates and rerun if there are.  Default is False.
    delvedir : str, optional
       The main delve directory.  Default is /net/dl2/dnidever/delve/ 
    delvereddir : str, optional
       The main delvered software directory.  Default is /home/dnidever/projects/delvered/'. 
 
    Returns
    -------
    PHOTRED_ALLFRAME will be run on each brick and a final catalog and 
    summary file created for each brick. 
 
    Example
    -------

    bricks(100)
    
    By D. Nidever  Aug 2019 
    Translated to Python by D. Nidever,  April 2022
    """

    t0 = time.time() 
     
    # Defaults
    if delvedir is not None:
        if delvedir.endswith('/')==False:
            delvedir += '/'
    else: 
        delvedir = '/net/dl2/dnidever/delve/'
    if delvereddir is not None:
        if delvereddir.endswith('/')==False:
            delvereddir += '/'
    else: 
        delvereddir = '/home/dnidever/projects/delvered/' 
    # Exposures directory 
    expdir = delvedir+'exposures/' 
    # Bricks directory 
    brickdir = delvedir+'bricks/' 
    # Logs directory 
    logsdir = expdir+'logs/' 
    if os.path.exists(logsdir)==False:
        os.makedirs(logsdir)
    tempdir = '/tmp/'
    if workdir is None:
        host = socket.gethostname()
        hostname = host.split('.')[0]
        workdir = '/data0/dnidever/delve/' 
        if workdir is not None:
            if os.path.exists(workdir)==False:
                os.makedirs(workdir)
     
    # Start the logfile 
    #------------------ 
    # format is delvered_bricks.host.DATETIME.log
    # Set up logging to screen and logfile
    logFormatter = logging.Formatter("%(asctime)s [%(levelname)-5.5s]  %(message)s")
    logger = logging.getLogger() 
    while logger.hasHandlers(): # some existing loggers, remove them   
        logger.removeHandler(logger.handlers[0]) 
    logger = logging.getLogger()
    logtime = datetime.now().strftime("%Y%m%d%H%M%S")
    host = socket.gethostname()
    hostname = host.split('.')[0]
    logfile = logsdir+'delvered_bricks.'+hostname+'.'+logtime+'.log'     
    if os.path.exists(logfile): os.remove(logfile)
    fileHandler = logging.FileHandler(logfile)
    fileHandler.setFormatter(logFormatter)
    rootLogger.addHandler(fileHandler)
    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    logger.addHandler(consoleHandler)
    logger.setLevel(logging.NOTSET)

    # Load the brick information 
    brickstr = Table.read(delvereddir+'data/delvemc_bricks_0.25deg.fits.gz') 
     
    # Parse the input nights
    bricknames = []
    for i in range(len(input)): 
        input1 = input[i] 
        # See if there is a - 
        # Use RA-DEC coordinate as a BOX selection 
        #   could also think of them as an ordered list and take all bricks 
        #   between the two (inclusive) 
        if input1.find('-') != -1: 
            brickrange = input1.split('-')
            br = 2*[None]  # RA and DEC components 
            for j in range(2):
                if brickrange[j].find('p')>-1:
                    br1 = brickrange[j].split('p')
                else:
                    br1 = brickrange[j].split('m')
                br[j] = np.array(br1).astype(float) / 10.0
            ind_bricks , = np.where((brickstr['ra'] >= br[0][0]) & (brickstr['ra'] <= br[1][0]) &
                                    (brickstr['dec'] >= br[0][1]) & (brickstr['dec'] <= br[1][1]))
            nind_bricks = len(ind_bricks)
            if nind_bricks > 0:
                bricknames += list(brickstr[ind_bricks])
            else: 
                logger.info('No brick found matching '+input1)
        else: 
            if input1[0] == '*':
                bricknames += list(brickstr['brickname'])
            else:
                ind1,ind2 = dln.match(brickstr['brickanme'],input1)
                if len(ind1)>0:
                    bricknames += list(brickstr['bricknames'][ind2])
                else: 
                    logger.info(input1+' brick not found')
                    nbricks = len(bricknames) 
    if nbricks == 0: 
        print('No bricks to process')
        return 

    bricknames = [b.strip() for b in bricknames]
    bricknames = np.unique(bricknames)
    nbricks = len(bricknames) 
 
    # Print info 
    #----------- 
    logger.info('')
    logger.info('############################################')
    logger.info('Starting DELVERED_BRICKS   ',systime(0))
    logger.info('Running on '+host)
    logger.info('############################################')
    logger.info('')
 
    # Check if bricks were previously done
    if redo==False and update==False:
        logger.info('Checking for any bricks were previously done')
        outfiles = brickdir+np.char.array(bricknames)+'/'+np.char.array(bricknames)+'.fits.gz'
        exists = [os.path.exists(f) for f in outfiles]
        bd, = np.where(np.array(exists)==False)
        gd, = np.where(np.array(exists)==True)        
        if len(gd) == 0: 
            logger.info('All bricks were previously processed.  Nothing to do')
            return 
        if len(bd) > 0: 
            logger.info('REDO not set and '+str(nbd)+' bricks were previously processed. Removing them from current list to process.')
            bricknames = bricknames[gd] 
            nbricks = ngd 
 
    # Check if we need to update the exposures database
    dbfile = '/net/dl2/dnidever/delve/bricks/db/delvered_summary.db' 
    lockfile = dbfile+'.lock' 
    logger.info('Checking if exposures database needs to be updated')
    # Check for lockfile 
    while (os.path.exists(lockfile) == 1): 
        logger.info('Lock file found. Waiting 60 seconds.')
        time.sleep(60)
    
    dbmtime = os.getmtime(dbfile)
    nightdirs = glob(delvedir+'exposures/20??????')
    sumfiles = np.array([ndir+'/'+os.path.basename(ndir)+'_summary.fits' for ndir in nightdirs])
    sumexists = np.array([os.path.exists(f) for f in sumfiles])    
    sumsize = np.array([os.path.getsize(f) for f in sumfiles])
    summtime = np.array([os.path.getmtime(f) for f in sumfiles])    
    gsum, = np.where((sumexists==True) & (sumsize>0))
    if len(gsum) > 0:
        if np.max(summtime[gsum]) > dbmtime: 
            logger.info('Need to update the database')
            dln.touch(lockfile)
            logger.info('Waiting 10 sec to let current queries finish')
            time.sleep(10)
            # This takes about 2.5 min to run
            out = subprocess.check_output(delvereddir+'bin/make_delvered_summary_table',shell=False)
            if os.path.exists(lockfile): os.remove(lockfile)

        
    ########################################## 
    ##  STARTING THE PROCESSING 
    ##########################################
    logger.info('Processing '+str(nbricks)+' brick(s)')
    cmd = "delvered_forcebrick,'"+np.char.array(bricknames)+"',delvedir='"+delvedir+"'" 
    if redo:
        cmd += ',/redo' 
    if update:
        cmd += ',/update' 
    cmddirs = np.char.array(np.zeros(nbricks,(np.str,200)))+workdir 
    jobs = job_daemon(cmd,cmddirs,prefix='dlvbrcks',hyperthread=True,nmulti=nmulti,wait=5)
 
    logger.info('DELVERED_BRICKS:ne after ',str(time.time()-t0,2),' sec')
 
 
