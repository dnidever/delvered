#!/usr/bin/env python

# ;+
# ;
# ; DELVERED_JOINTBRICKCATS
# ;
# ; Process a single DELVE brick and perform ALLFRAME FORCED photometry
# ;
# ;-

import numpy as np
from astropy.io import fits
import os
import time
import socket
import logging
from astropy.wcs import WCS


def delvered_jointbrickcats(brick, scriptsdir=None, irafdir=None, workdir=None, redo=False, logfile=None):
    # Start timing
    t0 = time.time()
    
    # Validate inputs
    if not brick:
        print("Syntax - delvered_jointbrickcats(brick, scriptsdir=scriptsdirs, irafdir=irafdir, workdir=workdir, redo=redo, logfile=logfile)")
        return
    
    # Defaults
    delvedir = '/net/dl2/dnidever/delve/'
    delvereddir = '/home/dnidever/projects/delvered/'
    scriptsdir = scriptsdir or '/home/dnidever/projects/PHOTRED/scripts/'
    irafdir = irafdir or '/home/dnidever/iraf/'
    tempdir = '/tmp/'
    
    # Workdir setup
    if not workdir:
        host = socket.gethostname()
        hostname = host.split('.')[0]
        workdir = '/data0/dnidever/delve/'
        os.makedirs(workdir, exist_ok=True)

    # Directories
    expdir = os.path.join(delvedir, 'exposures/')
    brickdir = os.path.join(delvedir, 'bricks/')
    logsdir = os.path.join(brickdir, 'logs/')

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
    decam = np.genfromtxt(delvereddir+'data/decam.txt')
    
    # Load brick information
    brick_info_path = os.path.join(delvereddir, 'data/delvemc_bricks_0.25deg.fits.gz')
    with fits.open(brick_info_path) as hdul:
        brick_data = hdul[1].data

    # Find brick in data
    brick_index = np.where(brick_data['BRICKNAME'] == brick)[0]
    if len(brick_index) == 0:
        log_message(logfile, f"{brick} not in DELVE-MC brick list")
        return
    brick_info = brick_data[brick_index[0]]

    # Subdirectory
    subdir = os.path.join(brickdir, brick[:4])
    bdir = os.path.join(subdir, brick)
    os.makedirs(bdir, exist_ok=True)

    # Check output file
    joint_object_path = os.path.join(bdir, f"{brick}_joint_object.fits.gz")
    if os.path.exists(joint_object_path) and not redo:
        print(f"{joint_object_path} EXISTS and /redo NOT set")
        return

    # Imager settings
    this_imager = {"telescope": "BLANCO", "instrument": "DECam", "namps": 62, "separator": "_"}

    # Logging information
    log_info = {
        "BRICK": brick,
        "DELVEDIR": delvedir,
        "SCRIPTSDIR": scriptsdir,
        "IRAFDIR": irafdir,
        "EXPOSUREDIR": expdir,
        "WORKINGDIR": workdir,
        "HOST": os.getenv('HOST', 'Unknown Host'),
    }
    for key, value in log_info.items():
        log_message(logfile, f"{key} = {value}")


    # Print out some information about the brick
    log_info = {
        "RA" : brick_info['RA'],
        "DEC" : brick_info['DEC'],
        "RA range" : f"[{brick_info['RA1']}, {brick_info['RA2']}]",
        "DEC range": f"[{brick_info['DEC1']}, {brick_info['DEC2']}]",
    }
    for key, value in log_info.items():
        log_message(logfile, f"{key} = {value}")

    # Get the Brick WCS information 
    tile = make_brick_wcs(brickstr1)
    
    # Process DELVE and community MC data
    process_brick(brick_info, bdir, expdir, scriptsdir, redo, logfile)

    # Finalize and print elapsed time
    elapsed_time = time.time() - t0
    log_message(logfile, f"Completed in {elapsed_time:.2f} seconds")


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
    tile = {'type':'WCS','naxis':np.array([nx,ny]),'cdelt':np.array([step,step]).astype(float),
            'crpix':np.array([xref+1,yref+1]).astype(float),
            'crval':np.array([brickstr['ra'],brickstr['dec']]).astype(float),'ctype':['RA--TAN','DEC--TAN'],
            'head':tilehead,'wcs':wcs,'xrange':[0,nx],'yrange':[0,ny],'nx':nx,'ny':ny}

    return tile


def log_message(logfile, message):
    """Log a message to the console and optionally to a file."""
    print(message)
    if logfile and logfile != -1:
        with open(logfile, 'a') as log:
            log.write(message + '\n')


def process_brick(brick_info, bdir, expdir, scriptsdir, redo, logfile):
    """
    Process the DELVE brick data.
    """

    


# Example usage
if __name__ == "__main__":
    delvered_jointbrickcats(
        brick="1234m045",
        scriptsdir="/path/to/scripts/",
        irafdir="/path/to/iraf/",
        workdir="/path/to/workdir/",
        redo=True,
        logfile="/path/to/logfile.log"
    )
