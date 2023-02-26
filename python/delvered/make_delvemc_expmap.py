import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import SymLogNorm
from gala.coordinates import MagellanicStreamNidever08
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.table import Table, Column, join, unique, setdiff, MaskedColumn
from exposure_map_utils import decam_array_map, tot_2_eff
from noirlab_archive_query import inside_delvemc

def get_processing_data():
    '''
    CHANGE THIS FUNCTION TO USE NEWER PROCESSED DATA.

    Joins data from DELVERED processing together with Noirlab Archive data.
    '''

    processing_data = Table.read('/home/polmassana/Documents/DELVE/DELVEMC_data/mother_summary_20221005.fits')
    #archive_hdus = Table.read('/home/polmassana/Documents/DELVE/DELVEMC_data/delvemc_info_20120901_20220930.fits')
    archive_data = Table.read('/home/polmassana/Documents/DELVE/DELVEMC_data/delvemc_info_20120901_20221009_light.fits')

    tmp=Table()
    tmp.add_column(Column(processing_data['EXPNUM'], name='EXPNUM', dtype=int))
    processing_data.update(tmp)
    processing_data = unique(processing_data, keys='EXPNUM')

    archive_data['G-TRANSP'][np.isnan(archive_data['G-TRANSP'])] = 0.7
    archive_data['G-TRANSP'].fill_value = 0.7

    comb_table = join(processing_data,archive_data.filled(), keys='EXPNUM', join_type='inner')

    return comb_table

def add_teff_col(orig_table):
    '''
    Adds the teff data for processed exposures, so that it can be joined with the data from the DES pipeline.
    '''

    teffs = []
    for row in orig_table:
        teffs.append(tot_2_eff(row['exposure'],row['SKYMODE'],row['G-TRANSP'],row['FWHM'],filter=row['FILTER']))

    orig_table.add_column(Column(np.array(teffs), name='teff'))

    return orig_table


def add_non_processed(processed_table):
    '''
    Adds teff form the DES pipeline to the already processed table for exposures that exist but have not been processed yet.
    '''

    alex_data = Table.read('/home/polmassana/Documents/DELVE/DELVEMC_data/decam-exposures-20220926.fits.gz') #CHANGE THIS LINE FOR YOUR OWN PATH

    alex_coords = SkyCoord(alex_data['ra']*u.deg,alex_data['dec']*u.deg)
    alex_mc = alex_data[inside_delvemc(alex_coords)]
    alex_mc = alex_mc[alex_mc['exptime']>=60]
    alex_mc.rename_column('expnum','EXPNUM')
    #print(alex_mc)

    missing_exps = setdiff(alex_mc, processed_table, keys='EXPNUM')

    processed_table.rename_columns(['ra_center','dec_center', 'FILTER'],['ra','dec', 'filter'])
    missing_exps.add_column(MaskedColumn(missing_exps['teff']*missing_exps['exptime'],name='summed_teff'))
    final_missing = missing_exps['EXPNUM','ra','dec','summed_teff']
    final_missing.rename_column('summed_teff','teff')
    final_missing.fill_value = 0

    missing_filtrs = []
    for fltr in missing_exps['filter']:
        missing_filtrs.append(fltr[0])
    final_missing.add_column(Column(np.array(missing_filtrs), name='filter'))

    return join(processed_table['EXPNUM','ra','dec','teff','filter'], final_missing['EXPNUM','ra','dec','teff','filter'].filled(), keys=['EXPNUM','ra','dec','teff','filter'], join_type='outer')

def run_expmap_vs_SMASH(data):
    '''
    Based on a Table containing coordinates and exposure times, return arrays for all the exposure maps of g, r and i bands compared to the SMASH reference.
    '''

    g_mask = data['filter']=='g'
    r_mask = data['filter']=='r'
    i_mask = data['filter']=='i'
    mask_list = [g_mask, r_mask, i_mask]
    field_coords = SkyCoord(data['ra']*u.degree, data['dec']*u.degree)
    field_coords_mag = field_coords.transform_to(MagellanicStreamNidever08())

    g_expmap, r_expmap, i_expmap = np.array([]), np.array([]), np.array([])
    expmap_list =[g_expmap, r_expmap, i_expmap]

    for i, mask in enumerate(mask_list):
        expmap_list[i] = decam_array_map(field_coords_mag.L.deg[mask],field_coords_mag.B.deg[mask],data['teff'][mask])

    g_teff_thresh = 131
    r_teff_thresh = 298
    i_teff_thresh = 490
    thresh_list = [g_teff_thresh,r_teff_thresh,i_teff_thresh]

    for expmap,thresh in zip(expmap_list,thresh_list):
        expmap[expmap>0] -= thresh
        expmap[expmap==0] = np.nan

    return expmap_list

def make_DELVEvsSMASH_plot(expmap_list, name='DELVE-MCvsSMASH.png'):
    '''
    Given a list of exposure map arrays, make the SMASH vs DELVE-MC plot for g, r and i filters.

    Parameters:
    -----------
    expmap_list: list
        List of 2D arrays from run_expmap_vs_SMASH(), order for the bands needs to be g,r and i.
    name: str
        Name for the file to be saved, or with its full path.

    Output:
    -------
    Saves a file with the plot.

    '''

    range_min = -33
    range_max = 30

    fig = plt.figure(figsize = (14.5,5.5), facecolor='white')
    gs = GridSpec(1, 3, figure=fig)
    gss1 = gs[0,0].subgridspec(1,2, wspace=0, width_ratios=[0.95,0.05])
    gss2 = gs[0,1].subgridspec(1,2, wspace=0, width_ratios=[0.95,0.05])
    gss3 = gs[0,2].subgridspec(1,2, wspace=0, width_ratios=[0.95,0.05])
    ax1 = fig.add_subplot(gss1[0])
    ax1c = fig.add_subplot(gss1[1])
    ax2 = fig.add_subplot(gss2[0])
    ax2c = fig.add_subplot(gss2[1])
    ax3 = fig.add_subplot(gss3[0])
    ax3c = fig.add_subplot(gss3[1])
    axs_list = [ax1, ax2, ax3]
    axsc_list = [ax1c, ax2c, ax3c]
    filter_list = ['g', 'r', 'i']

    fig.suptitle('SMASH vs DELVE \n Effective exposure time difference (s)')

    for ax, axc, filt, expmap in zip(axs_list, axsc_list, filter_list, expmap_list):
        mappable = ax.imshow(expmap.T, extent=(range_min, range_max, range_min, range_max), aspect='auto',
                             cmap = 'RdBu', norm=SymLogNorm(linthresh=10, linscale=0.1, vmin=-250.0, vmax=250.0, base=10), origin='lower', rasterized=True, interpolation='none')
        plt.colorbar(mappable, cax=axc)
        ax.annotate(r'$'+filt+'$', (0.9, 0.9), xycoords = 'axes fraction')
        ax.set_xlim(26,-33)
        ax.set_ylim(-30,30)

    ax1.set_ylabel(r'B$_{\mathrm{MS}}$ (deg)')
    ax1.set_xlabel(r'L$_{\mathrm{MS}}$ (deg)')
    ax2.set_xlabel(r'L$_{\mathrm{MS}}$ (deg)')
    ax3.set_xlabel(r'L$_{\mathrm{MS}}$ (deg)')

    plt.tight_layout()
    fig.savefig(name)

if __name__=="__main__":

    filter_list = ['g', 'r', 'i']
    expmap_list = []
    for fltr in filter_list:
        expmap_list.append(np.genfromtxt(f'../../data/exposuremap_vsSMASH_{fltr}.txt'))
    make_DELVEvsSMASH_plot(expmap_list)
