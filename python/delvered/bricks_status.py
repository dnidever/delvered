#!/usr/bin/env python
#
# Script to run doppler.fit() on a spectrum

from __future__ import print_function

import os
import sys
import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from datetime import datetime
from astropy.time import Time
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.table import Table
from dlnpyutils import utils as dln
from argparse import ArgumentParser
import traceback
try:
    import __builtin__ as builtins # Python 2
except ImportError:
    import builtins # Python 3

def statusmap(db=None,plot_nchips=False):
    # Plot distribution of jobs and save to file

    if db is None:
        from delvered import bricks_daemon
        db = bricks_daemon.db

    try:
        from gala.coordinates import MagellanicStreamNidever08
    except:
        print('You need to install the gala library from APW, this is going to crash.')
    print('Retrieving data...')
    all_data = Table(db.query(sql="select ra,dec,priority_points from delvered_processing.bricks"))
    done_data = Table(db.query(sql="select ra,dec from delvered_processing.bricks where status='DONE'"))
    redo_data = Table(db.query(sql="select ra,dec from delvered_processing.bricks where status='REDO'"))
    crashed_data = Table(db.query(sql="select ra,dec from delvered_processing.bricks where status='CRASHED'"))
    noupdate_data = Table(db.query(sql="select ra,dec from delvered_processing.bricks where status='NOUPDATE'"))
    running_data = Table(db.query(sql="select ra,dec from delvered_processing.bricks where status='R'"))
        
    all_coords = SkyCoord(all_data['ra']*u.deg,all_data['dec']*u.deg).transform_to(MagellanicStreamNidever08())

    print('Making the plot...')
        
    mpl.rcParams['font.size'] = 16
    fig, ax = plt.subplots(figsize=(10,8), facecolor='white')

    if plot_nchips:
        bkg_map = ax.scatter(all_coords.L.deg, all_coords.B.deg, s=8, marker='.',c=all_data['nchips'],
                             cmap='viridis',norm=LogNorm())
        plt.colorbar(bkg_map,ax=ax,label='Nchips')
    else:
        bkg_map = ax.scatter(all_coords.L.deg, all_coords.B.deg, s=8, marker='.',c=all_data['priority_points'],cmap='viridis')
        plt.colorbar(bkg_map,ax=ax,label='Priority')

    if plot_nchips==False:
        if len(done_data)>0:
            done_coords = SkyCoord(done_data['ra']*u.deg,done_data['dec']*u.deg).transform_to(MagellanicStreamNidever08())
            ax.scatter(done_coords.L.deg, done_coords.B.deg, s=8, marker='.',c='purple', label='DONE ('+str(len(done_data))+')')
        if len(redo_data)>0:
            redo_coords = SkyCoord(redo_data['ra']*u.deg,redo_data['dec']*u.deg).transform_to(MagellanicStreamNidever08())
            ax.scatter(redo_coords.L.deg, redo_coords.B.deg, s=8, marker='.',c='orange', label='REDO ('+str(len(redo_data))+')')
        if len(noupdate_data)>0:
            noupdate_coords = SkyCoord(noupdate_data['ra']*u.deg,noupdate_data['dec']*u.deg).transform_to(MagellanicStreamNidever08())
            ax.scatter(noupdate_coords.L.deg, noupdate_coords.B.deg, s=8, marker='.',c='#a00498', label='NOUPDATE ('+str(len(noupdate_data))+')')
        if len(crashed_data)>0:
            crashed_coords = SkyCoord(crashed_data['ra']*u.deg,crashed_data['dec']*u.deg).transform_to(MagellanicStreamNidever08())
            ax.scatter(crashed_coords.L.deg, crashed_coords.B.deg, s=8, marker='.',c='red', label='CRASHED ('+str(len(crashed_data))+')')
        if len(running_data)>0:
            running_coords = SkyCoord(running_data['ra']*u.deg,running_data['dec']*u.deg).transform_to(MagellanicStreamNidever08())
            ax.scatter(running_coords.L.deg, running_coords.B.deg, s=8, marker='.',c='#0bf9ea', label='RUNNING ('+str(len(running_data))+')')

    ax.set_xlabel(r'$L_{\rm MS}$ (deg)')
    ax.set_ylabel(r'$B_{\rm MS}$ (deg)')
    ax.set_xlim(26,-42)
    ax.set_ylim(-30,30)
    if plot_nchips==False:
        ax.legend(loc='upper right',markerscale=6)

    rightnow = datetime.now().strftime("%Y-%m-%dT%H:%M:%S")
    plotdir = '/net/dl2/dnidever/delve/bricks/status_plots/'
    if plot_nchips:
        plotfile = plotdir+'brick_status_nchips_'+rightnow+'.png'
    else:
        plotfile = plotdir+'brick_status_'+rightnow+'.png'
    fig.savefig(plotfile,bbox_inches='tight')
    print('Plot saved to: '+plotfile)


def runstatus(db=None):
    # Detailed information on running jobs

    if db is None:
        from delvered import bricks_daemon
        db = bricks_daemon.db

    print('Delvered Brick Processing Status at '+datetime.now().ctime())
    print('{:d} RUNNING'.format(nrunning))
    print('--------------------------------------------------------------------------------------------------------------------')
    print(' IND  BRICK  JOBID    RA    DEC  Rlmc  Rsmc Priority    Start Time      Runtime  Host   User           Logfile')
    print('--------------------------------------------------------------------------------------------------------------------')
    for i in range(nrunning):
        rtab1 = rtab[i]
        runid = rtab1['runid']
        jobid = rtab1['runjobid']
        runstart = rtab1['runstart'][0:21]
        host = rtab1['runhost'].split('.')[0]
        dt = (Time(datetime.now())-Time(runstart)).sec
        if dt < 60:
            dtime = '{:2d}s'.format(int(np.rint(dt)))
        elif dt >= 60 and dt < 3600:
            dtime = '{:4.1f}m'.format(dt/60)
        elif dt >= 3600 and dt < 24*3600:
            dtime = '{:4.1f}h'.format(dt/3600)
        else:
            dtime = '{:4.1f}d'.format(dt/3600/24.)
        data = (i+1,rtab1['brickname'],jobid,rtab1['ra'],rtab1['dec'],rtab1['lmc_radius'],
                rtab1['smc_radius'],rtab1['priority'],runstart,dtime,host,
                rtab1['runuser'],rtab1['logfile'])
        fmt = '{:3d} {:8s} {:6d} {:6.2f} {:6.2f} {:4.1f} {:4.1f} {:5d} {:23s} {:>4s} {:>6s} {:8s} {:60s}'
        print(fmt.format(*data))
    print('----------------------------------------------------------------------------------------------------------')


def status(running_flag=False,plot_flag=False):

    from delvered import bricks_daemon
    db = bricks_daemon.db

    rtab = db.query(sql="select * from delvered_processing.bricks where status='R'")
    nrunning = len(rtab)

    if plot_flag:
        mapfile = statusmap(db)

    if running_flag:
        runstatus(db)

    # Summary information
    ndone = db.query(sql="select count(*) from delvered_processing.bricks where status='DONE'")
    ntodo = db.query(sql="select count(*) from delvered_processing.bricks where status='TODO'")
    nredo = db.query(sql="select count(*) from delvered_processing.bricks where status='REDO'")
    ncrashed = db.query(sql="select count(*) from delvered_processing.bricks where status='CRASHED'")
    # Print out the information
    print('Delvered Brick Processing Status at '+datetime.now().ctime())
    print('{:6d} DONE'.format(ndone[0][0]))
    print('{:6d} RUNNING'.format(nrunning))
    print('{:6d} CRASHED'.format(ncrashed[0][0]))
    print('{:6d} REDO'.format(nredo[0][0]))
    print('{:6d} TODO'.format(ntodo[0][0]))
    
    # Find the groups of processes that are running
    print(' ')
    print('Jobs running ('+str(nrunning)+'):')
    jid = np.char.array(rtab['runhost'])+'-'+np.char.array(rtab['runuser'])
    index = dln.create_index(jid)
    for i in range(len(index['value'])):
        ind = index['index'][index['lo'][i]:index['hi'][i]+1]
        nind = len(ind)
        val = index['value'][i]
        host = val.split('-')[0].split('.')[0]
        user = val.split('-')[1]
        print(' {:3d} {:s} {:s} '.format(nind,host,user))

        
    db.close()

