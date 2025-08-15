import os
import numpy as np
from dlnpyutils import utils as dln,coords
from astropy.table import Table
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from glob import glob
import traceback

def loadheader(headfile):
    """ Load header file """

    # FITS file
    if dln.file_isfits(headfile):
        hdu = fits.open(headfile)
        headdict = {}
        # Loop over the extendions
        for i in range(len(hdu)):
            head = hdu[i].header
            if i==0:
                if headfile.endswith('.fits.fz'):
                    headdict['main'] = hdu[0].header
                else:
                    # the main header is blank
                    # use the next header which should have filter and other info
                    headdict['main'] = hdu[1].header
            else:
                ccdnum = head['CCDNUM']
                headdict[ccdnum] = head
        hdu.close()
    # ASCII file
    else:
        headlines = dln.readlines(headfile)
        lo = dln.grep(headlines,'^SIMPLE  =',index=True)
        begind = dln.grep(headlines,'^XTENSION',index=True)
        begind = lo+begind
        endind = dln.grep(headlines,'^END',index=True)
        if len(endind) != len(begind):
            endind = np.array(begind)[1:]-1
            endind = np.concatenate((endind,[len(headlines)-1]))
        headdict = {}
        # Loop over the extendions                                                                                          
        for i in range(len(begind)):
            lines = headlines[begind[i]:endind[i]+1]
            head = fits.Header.fromstring('\n'.join(lines),sep='\n')
            if i==0:
                headdict['main'] = head
            else:
                ccdnum = head['CCDNUM']
                headdict[ccdnum] = head
    return headdict

def checkchip(chipfile):
    """ Check chip for weird bad image issue """
    # input filename for chip fits file
    fitsheadfile = chipfile+'.head'
    head = fits.Header.fromtextfile(fitsheadfile)
    resourcefile = os.path.join(os.path.dirname(chipfile),'.'+os.path.basename(chipfile))
    rlines = dln.readlines(resourcefile)
    mssfile = rlines[0].split('=')[1].split('[')[0].strip()
    extname = rlines[0].split('=')[1].split('[')[1][:-1]
    try:
        fhead = fits.getheader(mssfile,extname=extname)
    except:
        traceback.print_exc()
        return False,0.0
    radiff = head['crval1']-fhead['crval1']
    decdiff = head['crval2']-fhead['crval2']
    dist = coords.sphdist(head['crval1'],head['crval2'],fhead['crval1'],fhead['crval2'])

    if dist > 0.01:
        bad = True
    else:
        bad = False

    return bad,dist

def checkexposure(fdir,expnum,verbose=False):

    files = glob(fdir+'/chip??/F*-'+expnum+'_??.fits')
    files.sort()
    if verbose:
        print(len(files),'chip files')
    bad = len(files)*[None]
    dist = len(files)*[None]
    for i in range(len(files)):
        chipfile = files[i]
        bad1,dist1 = checkchip(chipfile)
        bad[i] = bad1
        dist[i] = dist1
        if verbose:
            print(i+1,chipfile,bad1,dist1)

    if np.sum(bad) > 10:
        bad = True
    else:
        bad = False

    return bad,np.sum(bad)

def checkexposure2(fdir,expnum,verbose=False):

    files = glob(fdir+'/chip??/F*-'+expnum+'_??.fits')
    files.sort()
    if verbose:
        print(len(files),'chip files')

    # Load the first resource file to get the mss filename
    fitsheadfile = files[0]+'.head'
    head = fits.Header.fromtextfile(fitsheadfile)
    resourcefile = os.path.join(os.path.dirname(files[0]),'.'+os.path.basename(files[0]))
    rlines = dln.readlines(resourcefile)
    mssfile = rlines[0].split('=')[1].split('[')[0].strip()
    extname = rlines[0].split('=')[1].split('[')[1][:-1]
    # Load the full image header file
    headerfile = mssfile.replace('.fits.fz','.hdr')
    #headlines = dln.readlines(headerfile)
    try:
        headdict = loadheader(headerfile)
    except:
        traceback.print_exc()
        return False,-1

    bad = len(files)*[False]
    dist = len(files)*[0.0]
    for i in range(len(files)):
        chipfile = files[i]
        base = os.path.basename(files[i])
        ccdnum = int(base.split('_')[-1][:-5])
        fitsheadfile = chipfile+'.head'
        head = fits.Header.fromtextfile(fitsheadfile)
        fhead = headdict.get(ccdnum)
        if fhead is None:
            bad[i] = False
            dist[i] = 0.0
            continue
        fcrval1 = fhead.get('crval1')
        fcrval2 = fhead.get('crval2')
        if fcrval1 is None or fcrval2 is None:
            bad[i] = False
            dist[i] = 0.0
            continue
        radiff = head['crval1']-fcrval1
        decdiff = head['crval2']-fcrval2
        dist1 = coords.sphdist(head['crval1'],head['crval2'],fhead['crval1'],fhead['crval2'])

        if dist1 > 0.01:
            bad1 = True
        else:
            bad1 = False

        #chipfile = files[i]
        #bad1,dist1 = checkchip(chipfile)
        bad[i] = bad1
        dist[i] = dist1
        if verbose:
            print(i+1,chipfile,bad1,dist1)

    if np.sum(bad) > 10:
        bad = True
    else:
        bad = False

    return bad,np.sum(bad)


def checknight(nightdir):
    night = os.path.basename(os.path.abspath(nightdir))[:8]
    outfile = nightdir+'/'+night+'_badimagecheck.fits'
    if os.path.exists(outfile):
        print(outfile,'exists already')
        return None

    expfile = os.path.join(nightdir,night+'_exposures.fits')
    exptab = Table.read(expfile)
    exptab['bad'] = False
    fieldlines = dln.readlines(os.path.join(nightdir,'fields'))

    # Get the exposures in each field directory
    fdirs1 = glob(os.path.join(nightdir,'F?'))
    fdirs1.sort()
    fdirs2 = glob(os.path.join(nightdir,'F??'))
    fdirs2.sort()
    fdirs = fdirs1+fdirs2
    data = []
    for i in range(len(fdirs)):
        fitsfiles = glob(fdirs[i]+'/chip01/F*-*_01.fits')
        expnum = [os.path.basename(f).split('-')[1].split('_')[0] for f in fitsfiles]
        data1 = [(e,fdirs[i]) for e in expnum]
        print(data1)
        data += data1
    dt = [('expnum',str,10),('fdir',str,100)]
    data = np.array(data,dtype=np.dtype(dt))  # convert to numpy structured array

    # Add that information to the exptab table
    _,ind1,ind2 = np.intersect1d(exptab['EXPNUM'],data['expnum'],return_indices=True)
    exptab['fdir'] = 100*' '
    exptab['fdir'][ind1] = data['fdir'][ind2]
    exptab.sort('EXPNUM')

    # Now check the xposures
    exptab['bad'] = False
    exptab['nbad'] = 0
    for i in range(len(exptab)):
        fdir = exptab['fdir'][i]
        expnum = exptab['EXPNUM'][i]
        #bad,nbad = checkexposure(fdir,expnum)
        bad,nbad = checkexposure2(fdir,expnum)
        exptab['bad'][i] = bad
        exptab['nbad'][i] = nbad
        print(i+1,expnum,exptab['bad'][i],exptab['nbad'][i])

    exptab.write(nightdir+'/'+night+'_badimagecheck.fits')

    return exptab

def checkallnights(dirs=None):
    if dirs is None:
        dirs = glob('/net/dl2/dnidever/delve/exposures/20??????')
        dirs.sort()
    print(len(dirs),'nights')
    for i in range(len(dirs)):
        print(i+1,'checking',dirs[i])
        dum = checknight(dirs[i])
