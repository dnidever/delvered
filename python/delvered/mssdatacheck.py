import os
from astropy.io import fits
import subprocess
import tempfile

def check(filename):
    tid,tfile = tempfile.mkstemp(prefix='chk')
    os.remove(tfile)
    #res = subprocess.check_output(['funpack','-E','60','-O',tfile,filename],shell=False,stderr=subprocess.STDOUT)
    res = subprocess.run(['funpack','-E','60','-O',tfile,filename],shell=False,stderr=subprocess.DEVNULL,stdout=subprocess.DEVNULL)
    okay = os.path.exists(tfile)
    if okay: os.remove(tfile)
    return okay
