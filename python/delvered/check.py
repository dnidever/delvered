
# Check exposures to see if they are corrupted

import os
import numpy as np
from glob import glob
from astropy.table import Table
from astropy.io import fits
from dlnpyutils import utils as dln
import subprocess

# get rid of annoying warnings
import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', AstropyWarning)

def fitscheck(filename):
    """ Check that all of the data looks okay."""
    okay = True
    error = []
    if os.path.exists(filename)==False:
        return False,[str(filename)+' NOT FOUND']

    try:
        hdu = fits.open(filename)
        nhdu = len(hdu)
    except Exception as e:
        okay = False
        error.append('Error - '+str(e))
        try:
            hdu.close()
        except:
            pass
        return okay,error
    
    for i in range(nhdu):
        try:
            d = hdu[i].data
        except Exception as e:
            okay = False
            error.append('HDU '+str(i)+' error - '+str(e))
            break
    hdu.close()
    return okay,error

def collatefiles():
    """
    Get the list of all mass archive exposure files that we are using
    """

    dirs = glob('/net/dl2/dnidever/delve/exposures/20??????/F*')
    dirs = [d for d in dirs if os.path.isdir(d)]
    print(len(dirs),' directories to check')
    info = []
    badresourcefiles = []   # resources files with missing filenames
    for i,d in enumerate(dirs):
        field = os.path.basename(d)
        chfiles = glob(d+'/chip01/.'+field+'-????????_01.fits')
        print(i,d,len(chfiles))
        for f in chfiles:
            expnum = os.path.basename(f).split('-')[1].split('_')[0]
            #try:
            lines = dln.readlines(f)
            fluxfile = lines[0].split('=')[1].split('[')[0].strip()
            wtfile = lines[1].split('=')[1].split('[')[0].strip()
            maskfile = lines[2].split('=')[1].split('[')[0].strip()
            if fluxfile=='' or wtfile=='' or maskfile=='':
                print('blank filename')
                badresourcefiles.append(f)
            info.append({'expnum':expnum,'filename':fluxfile,'filetype':'flux'})
            info.append({'expnum':expnum,'filename':wtfile,'filetype':'wt'})
            info.append({'expnum':expnum,'filename':maskfile,'filetype':'mask'})
            #except:
            #    pass

    # Turn it into a table
    ninfo = len(info)
    dt = [('expnum',str,10),('filename',str,100),('filetype',str,10)]
    data = np.zeros(ninfo,dtype=np.dtype(dt))
    data['expnum'] = [d['expnum'] for d in info]
    data['filename'] = [d['filename'] for d in info]
    data['filetype'] = [d['filetype'] for d in info]
    # Get unique files
    _,ui = np.unique(data['filename'],return_index=True)
    data = data[ui]
    print(len(data),' unique exposure files')

    return data,badresourcefiles

def checkexposures(files):
    """
    Check all DELVE-MC exposures to see if they are corrupted
    """

    dt = [('filename',str,100),('okay',bool)]
    results = np.zeros(len(files),dtype=np.dtype(dt))
    for i in range(len(files)):
        results['filename'][i] = files[i]
        res = subprocess.run(['checkfits',files[i]],capture_output=True,shell=False)
        returncode = res.returncode
        out = res.stdout.decode().strip()
        err = res.stderr.decode().strip()
        okay1 = True
        if returncode==0:
            if out != '':
                if out.split()[0]=='BAD':
                    okay1 = False
                error1 = out
            else:
                okay1 = False
                error1 = ['Unknown problem']
        else:
            okay1 = False
            error1 = 'checkfits failed'
        results['okay'][i] = okay1
        if okay1:
            cmt = 'OK'
            ecmt = ''
            print(i,cmt,"'"+files[i]+"'",ecmt)
        else:
            cmt = 'BAD'
            ecmt = error1
            print(i,ecmt)

    return results
