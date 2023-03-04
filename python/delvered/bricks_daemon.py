#!/usr/bin/env python
#
# DBJOB_DAEMON.PY - Job manager running off of a database table
#

from __future__ import print_function

__authors__ = 'David Nidever <dnidever@montana.edu>'
__version__ = '20230301'  # yyyymmdd

# Modifcation of job_daemon.py

import os
import sys
import numpy as np
import warnings
import socket
import time
import subprocess
import tempfile
from datetime import datetime
from dlnpyutils import utils as dln
from astropy.io import fits
import psycopg2 as pg
from psycopg2.extras import execute_values

from psycopg2.extensions import register_adapter, AsIs
def addapt_np_float16(np_float16):
    return AsIs(np_float16)
def addapt_np_float32(np_float32):
    return AsIs(np_float32)
def addapt_np_float64(np_float64):
    return AsIs(np_float64)
def addapt_np_int8(np_int8):
    return AsIs(np_int8)
def addapt_np_int16(np_int16):
    return AsIs(np_int16)
def addapt_np_int32(np_int32):
    return AsIs(np_int32)
def addapt_np_int64(np_int64):
    return AsIs(np_int64)
def addapt_np_uint64(np_uint64):
    return AsIs(np_uint64)
def addapt_np_bool(np_bool):
    return AsIs(np_bool)
register_adapter(np.float16, addapt_np_float16)
register_adapter(np.float32, addapt_np_float32)
register_adapter(np.float64, addapt_np_float64)
register_adapter(np.int8, addapt_np_int8)
register_adapter(np.int16, addapt_np_int16)
register_adapter(np.int32, addapt_np_int32)
register_adapter(np.int64, addapt_np_int64)
register_adapter(np.uint64, addapt_np_uint64)
register_adapter(np.bool, addapt_np_bool)
register_adapter(np.bool_, addapt_np_bool)

from psycopg2.extensions import register_type
def cast_date(value, cursor):
    return value
oids = (1082, 1114, 1184) 
new_type = pg.extensions.new_type(oids, "DATE", cast_date)
register_type(new_type) 
# Cast None's for bool
oids_bool = (16,)
def cast_bool_none(value, cursor):
    if value is None:
        return False
    return np.bool(value)
new_type = pg.extensions.new_type(oids_bool, "BOOL", cast_bool_none)
register_type(new_type)
# Cast None's for text/char
oids_text = (18,25)
def cast_text_none(value, cursor):
    if value is None:
        return ''
    return value
new_type = pg.extensions.new_type(oids_text, "TEXT", cast_text_none)
register_type(new_type)
# Cast None's for integers
oids_int = (20,21,23)
def cast_int_none(value, cursor):
    if value is None:
        return -9999
    return np.int(value)
new_type = pg.extensions.new_type(oids_int, "INT", cast_int_none)
register_type(new_type)
# Cast None's for floats
oids_float = (700,701)
def cast_float_none(value, cursor):
    if value is None:
        return -999999.
    return np.float(value)
new_type = pg.extensions.new_type(oids_float, "FLOAT", cast_float_none)
register_type(new_type)


def register_date_typecasters(connection):
    """
    Casts date and timestamp values to string, resolves issues with out of
    range dates (e.g. BC) which psycopg2 can't handle
    """

    def cast_date(value, cursor):
        return value

    cursor = connection.cursor()
    cursor.execute("SELECT NULL::date")
    date_oid = cursor.description[0][1]
    cursor.execute("SELECT NULL::timestamp")
    timestamp_oid = cursor.description[0][1]
    cursor.execute("SELECT NULL::timestamp with time zone")
    timestamptz_oid = cursor.description[0][1]
    oids = (date_oid, timestamp_oid, timestamptz_oid)
    #oids = (1082, 1114, 1184)
    new_type = psycopg2.extensions.new_type(oids, "DATE", cast_date)
    pg.extensions.register_type(new_type) 


# Ignore these warnings, it's a bug
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

class DBSession(object):

    def __init__(self):
        """ Initialize the database session object. The connection is opened."""
        self.open()
        self.username = os.getlogin()
        self.pid = os.getpid()
        self.hostname = socket.gethostname()
        self.host = self.hostname.split('.')[0]
        self.runcounter = 1
        
    def open(self):
        """ Open the database connection."""
        connection = pg.connect(user="datalab",host="gp09.datalab.noirlab.edu",
                                password="",port = "5432",database = "tapdb")
        connection.autocommit = True
        self.connection = connection

    def close(self):
        """ Close the database connection."""
        self.connection.close()

    def restart(self):
        """ Restart the connection."""
        self.close()
        self.open()

    def nextbrick(self):
        """ Get the next brick."""
        cur = self.connection.cursor()
        # Loop until we have one that is TODO or REDO
        done = False
        niter = 0
        while (done==False):
            # Get the next vaue in the sequence
            cur.execute("select nextval('delvered_bricks_seq')")
            nextval = cur.fetchone()[0]
            # Now select that single row
            cur.execute('select * from delvered_processing.bricks where priority='+str(nextval))
            row = cur.fetchall()[0]
            brickname = row[0]
            brickid = row[1]
            status = row[6]
            if status=='TODO' or status=='REDO':
                done = True
            else:
                continue
            # RunID
            runid = self.username+'-'+self.host+'-'+str(self.pid)+'-'+str(self.runcounter)
            self.runcounter += 1
            # Now update the table to put our RUNID in there
            starttime = datetime.now().isoformat()
            cmd = "update delvered_processing.bricks set status='R',runid='"+runid+"',"
            cmd += "runstart='"+starttime+"',runhost='"+self.hostname+"',runuser='"+self.username+"'"
            cmd += " where priority="+str(nextval)
            cur.execute(cmd)
            niter += 1

        cur.close()
        return brickname,brickid,runid,status
            
    def query(self,sql,verbose=False,fmt='numpy',raw=False):
        """ Query the database."""
        cur = self.connection.cursor()
        # Execute the command
        if verbose:
            print('CMD = '+sql)
        cur.execute(sql)
        data = cur.fetchall()
        ndata = len(data)

        if len(data)==0:
            cur.close()
            return np.array([])

        # Return the raw results
        if fmt=='raw':
            cur.close()
            return data
    
        # Return fmt="list" format
        if fmt=='list':
            colnames = [desc[0] for desc in cur.description]
            data = [tuple(colnames)]+data
            cur.close()
            return data

        # Get table column names and data types
        colnames = [desc[0] for desc in cur.description]
        colnames = np.array(colnames)
        # Fix duplicate column names
        cindex = dln.create_index(colnames)
        bd,nbd = dln.where(cindex['num']>1)
        for i in range(nbd):
            ind = cindex['index'][cindex['lo'][bd[i]]:cindex['hi'][bd[i]]+1]
            ind.sort()
            nind = len(ind)
            for j in np.arange(1,nind):
                colnames[ind[j]] += str(j+1)

        # Use the data returned to get the type
        dt = []
        for i,c in enumerate(colnames):
            type1 = type(data[0][i])
            if type1 is str:
                dt.append( (c, type(data[0][i]), 300) )
            elif type1 is list:  # convert list to array
                nlist = len(data[0][i])
                dtype1 = type(data[0][i][0])
                dt.append( (c, dtype1, nlist) )
            else:
                dtype1 = type(data[0][i])
                # Check for None/null, need a "real" value to get the type
                if (data[0][i] is None):
                    cnt = 0
                    while (data[cnt][i] is None) and (cnt<(ndata-1)) and (cnt<100): cnt += 1
                    if data[cnt][i] is not None:
                        dtype1 = type(data[cnt][i])
                    else:  # still None, use float
                        dtype1 = float
                dt.append( (c, dtype1) )
        dtype = np.dtype(dt)

        # Convert to numpy structured array
        tab = np.zeros(len(data),dtype=dtype)
        tab[...] = data
        del(data)
        
        # For string columns change size to maximum length of that column
        dt2 = []
        names = dtype.names
        nplen = np.vectorize(len)
        needcopy = False
        for i in range(len(dtype)):
            type1 = type(tab[names[i]][0])
            if type1 is str or type1 is np.str_:
                maxlen = np.max(nplen(tab[names[i]]))
                dt2.append( (names[i], str, maxlen+10) )
                needcopy = True
            else:
                dt2.append(dt[i])  # reuse dt value

        # We need to copy
        if needcopy==True:
            dtype2 = np.dtype(dt2)
            tab2 = np.zeros(len(tab),dtype=dtype2)
            for n in names:
                tab2[n] = tab[n]
            tab = tab2
            del tab2

        # Convert to astropy table
        if fmt=='table':
            tab = Table(tab)
            
        return tab

    def setstatus(self,brickid,status):
        """ Set the status of a list of bricks."""
        if type(status) is not str:
            raise ValueError('Status must be a string')
        cur = self.connection.cursor()
        cmd = "update delvered_processing.bricks set status='"+status+"'"
        if dln.size(brickid)==1:
            cmd += " where brickid="+str(brickid)
        else:
            cmd += " where brickid IN ("+','.join(brickid)+")"
        cur.execute(cmd)
        cur.close()


# Set up the database connection
db = DBSession()


def mkstatstr(n=None):
    """ This returns the stat structure schema or an instance of the stat structure."""
    dtype = np.dtype([('jobid',str,20),('name',str,100),('user',str,100),
                      ('timeuse',str,100),('status',str,10),('queue',str,20)])
    if n is None:
        return dtype
    else:
        statstr = np.zeros(n,dtype=dtype)    
        return statstr


def mkjobstr(n=None):
    """ This returns the job structure schema or an instance of the job structure."""
    dtype = np.dtype([('host',str,20),('jobid',(str,100)),('input',str,1000),('dir',str,500),
                      ('name',str,100),('scriptname',str,200),('logfile',str),
                      ('submitted',bool),('done',bool),('begtime',float),('endtime',float),
                      ('duration',float),('brickname',str,100),('brickid',int),('success',bool)])
    if n is None:
        return dtype
    else:
        jobstr = np.zeros(n,dtype=dtype)
        jobstr['submitted'] = False
        jobstr['done'] = False
        return jobstr


def mkrunbatch():
    curdir = os.getcwd()
    batchfile = os.path.join(curdir,'runbatch')
    if os.path.exists(batchfile) is False:
        lines = []
        lines.append("#!/bin/sh\n")
        lines.append("if test $# -eq 0\n")
        lines.append("then\n")
        lines.append("  echo 'Syntax - runbatch program'\n")
        lines.append("else\n")
        lines.append("  echo 'Log file: '$1'.log'\n")
        lines.append("  ( nohup  $1 > $1.log 2>&1 ) &\n")
        lines.append("  echo 'JOBID='$!\n")
        lines.append("fi\n")
        dln.writelines(batchfile,lines,overwrite=True,raw=True)
        os.chmod(batchfile,0o755)
    return batchfile


def mkidlbatch():
    curdir = os.getcwd()
    batchfile = os.path.join(curdir,'idlbatch')
    if os.path.exists(batchfile) is False:
        try:
            out = subprocess.check_output(['which','idl'],stderr=subprocess.STDOUT,shell=False)
        except subprocess.CalledProcessError:
            raise Exception("IDL program not available")
        idlprog = out.decode().strip()
        if os.path.exists(idlprog) is False:
            raise Exception("IDL program "+idlprog+" not found")
        lines = []
        lines.append("#!/bin/sh\n")
        lines.append("if test $# -eq 0\n")
        lines.append("then\n")
        lines.append("  echo 'Syntax - idlbatch idl.batch'\n")
        lines.append("else\n")
        lines.append("  echo 'Log file: '$1'.log'\n")
        lines.append("  ( nohup "+idlprog+" < $1 > $1.log 2>&1 ) &\n")
        lines.append("  echo 'JOBID='$!\n")
        lines.append("fi\n")
        dln.writelines(batchfile,lines,overwrite=True,raw=True)
        os.chmod(batchfile,0o755)
    return batchfile


def check_diskspace(indir=None,updatestatus=False):
    """ This checks if there's enough free disk space. """
    if indir is None: raise ValueError("Must give directory")
    statvfs = os.statvfs(indir)
    available = statvfs.f_frsize * statvfs.f_bavail / 1e6  # in MB
    if updatestatus: print('Disk Space: '+str(available[0]),' MB available')
    # Not enough disk space available
    if available<100.:
        raise Exception('NOT enough disk space available')


def check_killfile(jobs=None,hyperthread=True):
    """ This checks for a kill file and if found kills all the active jobs. """
    if jobs is None: raise ValueError("No jobs structure input")
    killfile = '/tmp/'+os.getlogin()+'/killbricks'
    if os.path.exists(killfile) is True:
        sub, = np.where((jobs['submitted']==1) & (jobs['done']==0))
        nsub = len(sub)
        print('Kill file found.  Killing all '+str(nsub)+' job(s)')
        for i in range(nsub):
            # Killing the job
            print('Killing '+jobs['name'][sub[i]]+'  JobID='+jobs['jobid'][sub[i]])
            if hyperthread is False:
                out = subprocess.run(['qdel',jobs['jobid'][sub[i]]],stderr=subprocess.STDOUT,stdout=subprocess.PIPE,
                                     shell=False,check=False).stdout
            else:
                # Using negative PID will also kill subprocesses
                try:
                    out = subprocess.run(['kill','-9','-'+jobs['jobid'][sub[i]]],stderr=subprocess.STDOUT,stdout=subprocess.PIPE,
                                         shell=False,check=False).stdout
                except:
                    out = subprocess.run(['kill','-9',jobs['jobid'][sub[i]]],stderr=subprocess.STDOUT,stdout=subprocess.PIPE,
                                         shell=False,check=False).stdout
            # Change the status back in the database
            db.setstatus(jobs['brickid'][sub[i]]],'TODO')
        # Also kill any allframes that are running
        print('Killing any dangling allframe and IDL jobs')
        out = subprocess.run(['killall','allframe'],stderr=subprocess.STDOUT,stdout=subprocess.PIPE,
                             shell=False,check=False).stdout        
        out = subprocess.run(['killall','idl'],stderr=subprocess.STDOUT,stdout=subprocess.PIPE,
                             shell=False,check=False).stdout
        # Remove the kill file
        print('Deleting kill file "'+killfile+'"')
        os.remove(killfile)
        return True
    return False

def makescript(inp=None,indir=None,name=None,prefix=None,hyperthread=True,idle=False):
    """This makes job scripts for the dbjob_daemon program.

    Parameters
    ----------
    inp : string list or array
          The command to execute.  Can be an array.  Must input absolute path names.
          idlbatch will be used if it's an IDL command.  If the command is a
          series of unix commands separated by commas then these
          will be put on separate lines.
    indir : string list or array
          The directory to put the job script in.
    name : string list or array, optional
          The name to call the job script (without the '.sh' or '.batch' ending).
          If this is not provided then an autogenerated name will
          be made and returned.
    prefix : string, optional
          The prefix to use for the job script name. "pr" by default.
    hyperthread : bool, optional
          Not on a job server but one with multiple hyperthreaded
                   processors.
    idle : bool, optional
          This is an IDL command, otherwise a SHELL command.

    Returns
    -------
    scriptname : string array
               The absolute names of the scripts.
    Job scripts output to the directories and with the names specified.


    Example
    -------

    .. code-block:: python

        scriptname = makescript(inp,dir,name,hyperthread=True)

    """
  
    # Not enough inputs
    if (inp is None):
        raise ValueError('No input given')

    ninp = dln.size(inp)
    ndir = dln.size(indir)
    nname = dln.size(name)
    if type(name) is not list:
        name = [name]

    # Not enough directories input
    if (ndir>0) & (ndir!=ninp):
        raise ValueError('INPUT and DIRECTORIES are of different size')

    # Current directory
    curdir = os.getcwd()

    # No directories input
    if ndir==0: indir = np.repeat(curdir,ninp)
    if ndir==1: indir = np.repeat(indir,ninp)    # multiple commands in same dir

    # Construct names
    if (nname==0):
        name = np.zeros(ninp,dtype=(str,200))
        if prefix is not None:
            pre = dln.first_el(prefix)
        else:
            pre = 'job'
        for i in range(ninp):
            tid,tfile = tempfile.mkstemp(prefix=pre,dir=indir[i])
            os.close(tid)   # mkstemp opens the file, close it
            name[i] = os.path.basename(tfile)

    # Make scriptnames
    scriptname = np.zeros(nname,dtype=np.dtype('U100'))
    scriptname[:] = np.array(dln.strjoin(dln.pathjoin(indir,name),'.sh'),ndmin=1)
    # Script loop
    for i,input1 in enumerate(np.array(inp,ndmin=1)):
        base = str(name[i])
        sname = str(indir[i])+'/'+base+'.sh'

        #------
        # PBS
        #------
        if hyperthread is False:
            # IDL command
            if idle is True:
                # Making an IDL batch file
                bname = str(indir[i])+'/'+base+'.batch'
                dln.writelines(bname,input1,overwrite=True,raw=True)
                # The execution command
                cmd = 'idl < '+base+'.batch'
            # SHELL command
            else:
                # The execution command
                cmd = input1
                # If there are commas in the line then break it up into multiple lines
                if cmd.find(';') != -1:
                    cmd = cmd.replace(';','\n;')
                    cmd = cmd.split('j')

            # Make the command
            #----------------------
            lines = []
            lines.append('#!/bin/sh\n')
            lines.append('#PBS -l nodes=1:ppn=1\n')
            lines.append('#PBS -l walltime=96:00:00\n')
            lines.append('#PBS -o '+base+'.report.out\n')
            lines.append('#PBS -e '+base+'.error.out\n')
            #lines.append('#PBS -M dln5q@virginia.edu\n')
            lines.append('#PBS -V\n')
            lines.append('\n')
            lines.append('echo Running on host `hostname`\n')
            lines.append('echo Time is `date`\n')
            lines.append('echo "Nodes used for this job:"\n')
            lines.append('echo "------------------------"\n')
            lines.append('cat $PBS_NODEFILE\n')
            lines.append('echo "------------------------"\n')
            lines.append('\n')
            lines.append('cd '+indir[i]+'\n')
            for j in range(len(cmd)): lines.append(cmd[j]+'\n')
            lines.append('\n')
            lines.append('# print end time\n')
            lines.append('echo\n')
            lines.append('echo "Job Ended at `date`"\n')
            lines.append('echo\n')
            # Writing the file
            dln.writelines(scriptname[i],lines,overwrite=True,raw=True)
            # Print info
            print('PBS script written to: '+str(scriptname[i]))

        #----------------
        # Hyperthreaded
        #----------------
        else:
            # Just make batch file
            # treat shell and idl the same
            # The execution command
            cmd = input1
            # If there are commas in the line then break it up into multiple lines
            if cmd.find(';') != -1:
                cmd = cmd.replace(';','\n;')
                cmd = cmd.split(';')
            # IDL files should end in .batch
            if idle is True: scriptname[i] = str(indir[i])+'/'+base+'.batch'
            # Writing the file
            dln.writelines(scriptname[i],cmd,overwrite=True,raw=True)
            # Make SHELL scripts executable
            if idle is False: os.chmod(scriptname[i],0o755)
            # Print info
            print('HYPERTHREAD script written to: '+str(scriptname[i]))

    # Erase the temporary files that mkstemp makes
    dln.remove(dln.pathjoin(indir,name),allow=True)

    if dln.size(scriptname)==1: scriptname=scriptname[0]
    return scriptname


def submitjob(scriptname=None,indir=None,hyperthread=True,idle=False):
    """ This submits new jobs.

    Parameters
    ----------
    scriptname : string
           Name of the script to run.
    indir : string, optional
           Directory to change to before script is executed.
    hyperthread : bool, optional
          Not on a job server but one with multiple hyperthreaded
             processors.  By default, this is True.
    idle : bool, optional
         This is an IDL program.  By default his is False.

    Returns
    -------
    jobid : string
         The job ID for this job.
    logfile : string
         Name of the output log file for this job/script.
    The script `scriptname` will be run.

    Example
    -------

    .. code-block:: python

        jobid,logfile = submitjob(scriptname,indir,hyperthread=True)

    """
    if scriptname is None: raise ValueError("scriptname must be input")

    curdir = os.getcwd()
    # Submitting the job
    if hyperthread is False:
        try:
            out = subprocess.check_output('qsub '+scriptname,stderr=subprocess.STDOUT,shell=True)
        except subprocess.CalledProcessError:
            raise Exception("Problem submitting PBS job")
        jobid = dln.first_el(out)
        logfile = scriptname+'.log'
    else:
        if idle is True:
            batchprog = mkidlbatch()
        else:
            batchprog = mkrunbatch()
        if indir is not None: os.chdir(indir)
        try:
            out = subprocess.check_output(batchprog+' '+scriptname,stderr=subprocess.STDOUT,shell=True)
            #out = subprocess.run(batchprog+' '+scriptname,stderr=subprocess.STDOUT,stdout=subprocess.PIPE,
            #                     shell=True,check=False).stdout
        except:
            raise Exception("Problem submitting shell job")
        if indir is not None: os.chdir(curdir)
        # Get the JOBID
        out = out.decode().split('\n')
        jobid_ind = dln.grep(out,'^JOBID=',index=True)
        njobid_ind = len(jobid_ind)
        jobid = out[jobid_ind[0]].split('=')[1]
        jobid = jobid.strip()
        # Get the logfile
        logfile_ind = dln.grep(out,'^Log file: ',index=True)
        nlogfile_ind = len(logfile_ind)
        logfile = out[logfile_ind[0]].split(':')[1]
        logfile = logfile.strip()

    # Printing info
    print('Submitted '+scriptname+'  JobID='+jobid)

    return jobid, logfile


def getstat(jobid=None,hyperthread=True):
    """Get the status of jobs

    This checks the status of jobs for the job_daemon program.
    If no jobs are found in the queue then an empty
    statstr structure is returned.

    Parameters
    ----------
    jobid : string or int
          Specific JOBID to check.
    hyperthread : bool, optional
          Not on a PBS machine but one with multipe hyperthreaded
          processors running simultaneously.  Default is True.

    Results
    -------
    statstr : numpy structured array
          Status structure.

    Example
    -------

    .. code-block:: python

        stat = getstat(jobid)

    """
    if jobid is None: raise ValueError("Must input jobid")
    njobid = dln.size(jobid)

    # PBS
    #--------
    if hyperthread is False:
        if njobid>0: 
            out = subprocess.check_output(['qstat',jobid[0]],stderr=subprocess.STDOUT,shell=False)
        else:
            out = subprocess.check_output(['qstat'],stderr=subprocess.STDOUT,shell=False)
        nout = np.sum(out.strip() != '')
        gd = dln.grep(out,'^'+jobid,index=True)
        ngd = len(gd)
        if ngd>0:
            statlines = out[gd[0]]
            nstat = len(statlines)
        else:
            statlines = None
        # Some jobs in queue
        if statlines is not None:
            arr = statlines
            arr = arr.split(' ')
            statstr = mkstatstr(nstat)
            statstr['jobid'] = arr[0]
            statstr['name'] = arr[1]
            statstr['user'] = arr[2]
            statstr['timeuse'] = arr[3]
            statstr['status'] = arr[4]
            statstr['queue'] = arr[5]
        # No jobs in queue
        else:
            statstr = mkstatstr(1)

    # Hyperthreaded.  Need a jobid
    #-----------------------------
    else:
        # No JOBID input
        if njobid==0:
            print('Need JOBID with /hyperthread')
            return mkstatstr(1)

        try:
            out = subprocess.run(['ps','-o','pid,user,etime,command','-p',str(jobid)],
                                 stderr=subprocess.STDOUT,stdout=subprocess.PIPE,shell=False,check=False).stdout
        except:
            statstr = mkstatstr(1)
            statstr['jobid'] = jobid
            statstr['queue'] = 'hyperthread'
            return statstr
        # can put in the column that you want
        # ps -o etime -p jobid
        out = dln.strsplit(out.decode(),'\n')
        out = dln.strip(out)
        gd = dln.grep(out,'^'+str(jobid),index=True)
        ngd = len(gd)
        if ngd>0:
            statlines = np.array(out,ndmin=1)[gd[0]]
        else:
            statlines = None

        # Some jobs in queue
        if statlines is not None:
            arr = statlines.split()
            statstr = mkstatstr(1)
            statstr['jobid'] = arr[0]
            statstr['user'] = arr[1]
            # CAN'T get the name.
            statstr['timeuse'] = arr[2]
            statstr['status'] = 'R'
            statstr['queue'] = 'hyperthread'
        # No jobs in queue
        else:
            statstr = mkstatstr(1)
            statstr['jobid'] = jobid
            statstr['queue'] = 'hyperthread'
    return statstr


def checkjobs(jobs=None,hyperthread=True):
    """ Check on the status of the running jobs"""
    if jobs is None: raise ValueError("jobs must be input")
    sub, = np.where((jobs['submitted']==True) & (jobs['done']==False))
    nsub = len(sub)
    nfinished = 0
    for i in range(nsub):
        # Checking status
        jobid = jobs[sub[i]]['jobid']
        statstr = getstat(jobid,hyperthread=hyperthread)
        # Job done
        if statstr['status']=='':
            if nfinished==0: print('')
            print(time.ctime()+'  Input '+str(sub[i]+1)+' '+jobs['name'][sub[i]]+' JobID='+jobs['jobid'][sub[i]]+' FINISHED')
            jobs['done'][sub[i]] = True
            jobs['endtime'][sub[i]] = time.time()/3600/24   # in days
            jobs['duration'][sub[i]] = (jobs['endtime'][sub[i]] - jobs['begtime'][sub[i]])*3600*24   # in sec
            # Check if the job crashed
            #  there are some non-ascii characters, so we need to ignore the errors
            with open(jobs['scriptname'][sub[i]]+'.log','r',errors='ignore') as f:
                lines = f.readlines()
            # The log files should end with these four lines
            # Saving joint catalogs
            # Wed Mar  1 21:03:51 2023
            # Writing meta-data to /net/dl2/dnidever/delve/bricks/0329/0329m662/0329m662_joint_meta.fits
            # dt = 78.0 sec.
            # Success
            if lines[-1].startswith('dt = ') and  lines[-2].startswith('Writing meta-data to'):
                jobs['success'][sub[i]] = True
                print('Job {:} for brick {:} completed successfully'.format(jobs['jobid'][sub[i]],
                                                                                   jobs['brickname'][sub[i]]))
                db.setstatus(jobs['brickid'][sub[i]],'DONE')
            # No update
            elif lines[-1].startswith('Nothing to UPDATE'):
                jobs['success'][sub[i]] = True
                print('Job {:} for brick {:} had no new chips'.format(jobs['jobid'][sub[i]],
                                                                      jobs['brickname'][sub[i]]))
                db.setstatus(jobs['brickid'][sub[i]],'NOUPDATE')            
            # Crashed
            else:
                jobs['success'][sub[i]] = False
                print('Job {:} for brick {:} crashed'.format(jobs['jobid'][sub[i]],jobs['brickname'][sub[i]]))
                db.setstatus(jobs['brickid'][sub[i]],'CRASHED')
            # Check if there's a meta file and figure out the number of chips
            brickname = jobs['brickname'][sub[i]]
            metafile = '/net/dl2/dnidever/delve/bricks/'+brickname[0:4]+'/'+brickname+'/'+brickname+'_meta.fits'
            if os.path.exists(metafile):
                head = fits.getheader(metafile,1)
                nchips = head['NAXIS2']
                cur = db.connection.cursor()
                cur.execute("update delvered_processing.bricks set nchips='"+str(nchips)+"' where brickname='"+brickname+"'")
                cur.close()
            nfinished += 1
            # Check for errors as well!! and put in jobs structure
    return jobs


def status_update(jobs=None):
    """ Print out the status update."""
    if jobs is None: raise ValueError("jobs must be input")
    njobs = dln.size(jobs)                                                   # Number of total jobs
    n_inqueue = np.sum((jobs['submitted']==True) & (jobs['done']==False))    # Number of jobs still in queue  
    n_finished = np.sum(jobs['done']==True)                                  # Number of jobs finished 
    # Print the status
    print('')
    print(time.ctime())
    print('Jobs Summary: %d submitted, %d finished, %d running' % (njobs,n_finished,n_inqueue))


def daemon(scriptsdir=None,nmulti=4,waittime=0.2,statustime=60,redo=False):
    """
    This program is a job manager run off of a database table.

    NOTE:  If you want to "kill" all of the jobs create a file "killjobs"
    in the same directory that JOB_DAEMON is being run in and all of the
    PBS jobs will be killed.

    Parameters
    ----------
    scriptsdir : str, optional
       Directory for the scripts.
    nmulti : int, optional
       How many nodes to run these jobs on.  Default is 4.
    statustime : float or int, optional
       The time between status updates.  However, the status
        will always be updated if something has actually changed.
        Default is 60.
    waittime : float or int, optional
       Time to wait between checking the running jobs.  Default is 0.2 sec.
    redo : boolean, optional
       Redo the bricks.  Default is to use update.

    Results
    -------
    jobs : numpy structured array
       The jobs structure with information on the JOBID, script name, etc.
    Jobs are run in a batch mode.

    Example
    -------

    .. code-block:: python

        jobs = dbjob_daemon(hyperthread=True,nmulti=5)

    """

    idle = True
    hyperthread = True

    # Current directory
    curdir = os.getcwd()

    # Scripts directory
    if scriptsdir is None:
        username = os.getlogin()
        scriptsdir = '/tmp/'+username
        if os.path.exists(scriptsdir)==False:
            os.makedirs(scriptsdir)
        print('Putting scripts in '+scriptsdir)
    
    # Defaults
    if waittime<0.1: waittime=0.1
    if statustime<1: statustime=1

    # Host name
    hostname = socket.gethostname()
    host = hostname.split('.')[0]

    # Which IDL are we using?
    if idle is True:
        try:
            out = subprocess.check_output(['which','idl'],stderr=subprocess.STDOUT,shell=False)
        except subprocess.CalledProcessError:
            raise Exception("IDL program not available")
        idlprog = out.strip()
        if os.path.exists(idlprog) is False:
            raise Exception("IDL program "+idlprog+" not found")

    print('---------------------')
    print(' RUNNING DBJOB_DAEMON ')
    print('---------------------')
    print('Host='+host)
    print('Nmulti='+str(nmulti))
    
    t0 = time.time()
    timesincelaststatus = time.time()  # initializing the update time

    #--------
    # DAEMON
    #--------
    # -Keep submitting jobs until nmulti is reached
    # -Check every minute or so to see how many jobs are still
    #  running.  If it falls below nmulti and more jobs are left then
    #  submit more jobs
    # -Don't return until all jobs are done.

    # Initialize the "jobs" structure
    # id will be the ID from Pleione
    jobs = mkjobstr(1)
    jobs['host'] = host

    # Loop until all jobs are done
    # On each loop check the queue and figure out what to do
    count = 0
    endflag = False
    while (endflag==False):
        # Status update
        dtstatus_sec = time.time()-timesincelaststatus
        if (dtstatus_sec>statustime):
            updatestatus = True
            timesincelaststatus = time.time()
        else:
            updatestatus = False
  
        # Check disk space
        check_diskspace('/data0/')
        # Check for kill file
        if check_killfile(jobs) is True: return jobs

        # Check status of running jobs
        #-----------------------------
        jobs = checkjobs(jobs,hyperthread=hyperthread)

        # Status update
        #---------------
        if updatestatus is True: status_update(jobs)

        # Submit new jobs
        #----------------
        n_inqueue = np.sum((jobs['submitted']==True) & (jobs['done']==False))  # Number of jobs still in queue  
        nnew = nmulti-n_inqueue
        if (nnew>0):
            # Get the indices of new jobs to be submitted
            #nosubmit, = np.where(jobs['submitted']==False)
            #newind = nosubmit[0:nnew]
            # Update immediately if there are new jobs to submit
            print('')
            print(time.ctime())
            print('Updating Queue: '+str(n_inqueue)+' JOB(S) running, out of '+str(nmulti)+' Maximum. Submitting '+str(nnew)+' more job(s)')

            # Loop through the new submits
            for i in range(nnew):
                print('')
                # Get new brick from the database
                brickname,brickid,runid,dbstatus = db.nextbrick()
                name = 'dlvbrcks-'+runid
                if redo or dbstatus=='REDO':
                    cmd = "delvered_forcebrick,'"+brickname+"',/redo"
                else:
                    cmd = "delvered_forcebrick,'"+brickname+"',/update"
                #cmd = "print,'Testing for brick: "+brickname+"'"
                if idle is True:
                    print('Input '+str(len(jobs))+'  Command: >>IDL>'+str(cmd)+'<<')                    
                else:
                    print('Input '+str(len(jobs))+'  Command: >>'+str(cmd)+'<<')
                # Make script
                scriptname = makescript(cmd,indir=scriptsdir,name=name,
                                        hyperthread=hyperthread,idle=idle)
                name = os.path.basename(os.path.splitext(scriptname)[0])
                # Submitting the job
                jobid, logfile = submitjob(scriptname,scriptsdir,hyperthread=hyperthread,idle=idle)
                # Set the logfile in the database
                cur = db.connection.cursor()
                cur.execute("update delvered_processing.bricks set logfile='"+logfile+"' where brickname='"+brickname+"'")
                cur.close()
                newjob = mkjobstr(1)
                # Updating the jobs structure
                newjob['submitted'] = True
                newjob['jobid'] = jobid
                newjob['brickname'] = brickname
                newjob['brickid'] = brickid
                newjob['input'] = cmd
                newjob['name'] = name                
                newjob['dir'] = scriptsdir
                newjob['scriptname'] = scriptname
                newjob['logfile'] = logfile
                newjob['begtime'] = time.time()/3600/24   # in days
                jobs = np.hstack((jobs,newjob))
                # Wait a bit
                #  the first thing delvered_forcebrick.pro does it query the
                #  chips database.  we don't want to hammer it too hard.
                time.sleep(5)

        # Are we done?
        #-------------
        ndone = np.sum(jobs['done']==True)
        #if (ndone==njobs): endflag=True
        # Wait a bit
        #--------------
        if (endflag==0): time.sleep(waittime)
        # Increment the counter
        count += 1

    print('DONE')
    print('dt = {:.1f} sec'.format(time.time-t0))

    # Close the database connection
    cur.close()
    connection.close()
    
    return jobs
