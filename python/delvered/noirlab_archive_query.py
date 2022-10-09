import numpy as np
import requests
import json
from astropy.coordinates import SkyCoord, Galactic
import astropy.units as u
from astropy.table import Table, Column
from astropy import table
import warnings

def teff_calc(transp, fwhm, skyb, skybdark, exptime):
    # Function that returns the effective exposure time as calculated using Equation 2 from Nielsen et al. (2016).
    tau = ((transp*0.9)/fwhm)**2*(skybdark/skyb)
    return tau*exptime

def inside_delvemc(coords):
    # Function that returns whether a set of coordinates (using SkyCoord from astropy) in (RA,DEC) is inside the DELVE-MC footprint.
    # Returns a boolean array to be used as a mask.
    LMC = SkyCoord.from_name('LMC')
    SMC = SkyCoord.from_name('SMC')
    lmc_sep_deg = coords.separation(LMC).deg
    smc_sep_deg = coords.separation(SMC).deg
    coords_gal_lat = coords.transform_to(Galactic).b.deg
    return (((lmc_sep_deg<=25) | (smc_sep_deg<=15)) & (coords_gal_lat<=-10))

def query_fields(outfields, search_fields, limit=1000000):
    '''
    This function returns the response of the NOIRLab archive (https://astroarchive.noirlab.edu/)
    according to the search parameters specified in the input. The outfields returned belong to
    the core fields of each exposure. Returns an astropy Table if the search was succesful and an
    error message if it was not.

    Parameters
    ----------
    outfields: list
        The fields that the response will have (columns in the returned Table). Available fields
        can be found in https://astroarchive.noirlab.edu/api/fcatalog/.
    search_fields: list
        The fields that the query will search through. The available fields can be found here
        https://astroarchive.noirlab.edu/api/adv_search/fadoc/. Some fields (e.g. will often be
        required or it will raise error).
    limit: int, optional
        Limit of output rows. Default: 1000000

    Returns:
    --------
    response: Table or str
        Response from NOIRLab archive in Table format (if successful).
    '''
    natroot = "https://astroarchive.noirlab.edu"
    adsurl = f'{natroot}/api/adv_search'
    apiurl=f'{adsurl}/find/?limit={limit}'

    jj = {
        "outfields" : outfields,
        "search" : search_fields
    }
    response = requests.post(apiurl,json=jj)

    archive_result = None
    if response.status_code == 200:
        archive_result = Table(response.json()[1:])
        print(f'File records retrieved: {len(archive_result)}')
    else:
        print(response.json()['errorMessage'])
        print('Returning empty table...')
    return archive_result

def query_hdus(outfields, search_fields, limit=5000000):
    '''
    This function returns the response of the NOIRLab archive (https://astroarchive.noirlab.edu/)
    according to the search parameters specified in the input. The outfields returned belong to
    the metadata in the headers of each CCD. Returns an astropy Table if the search was succesful
    and an error message if it was not.

    Parameters
    ----------
    outfields: list
        The fields that the response will have (columns in the returned Table). Available fields
        can be found in https://astroarchive.noirlab.edu/api/fcatalog/.
    search_fields: list
        The fields that the query will search through. Some of the available fields can be found here
        https://astroarchive.noirlab.edu/api/adv_search/hadoc/, but also some of them will have to have
        hdu: appended to them at the beginning. Some fields (e.g. will often be required or it will
        raise error).
    limit: int, optional
        Limit of output rows. Default: 5000000

    Returns:
    --------
    response: Table or str
        Response from NOIRLab archive in Table format (if successful).
    '''
    natroot = "https://astroarchive.noirlab.edu"
    adsurl = f'{natroot}/api/adv_search'
    apiurl=f'{adsurl}/find/?limit={limit}&rectype=hdu'

    jj = {
        "outfields" : outfields,
        "search" : search_fields
    }

    response = requests.post(apiurl,json=jj)

    archive_result = None
    if response.status_code == 200:
        archive_result = Table(response.json()[1:])
        print(f'HDU records retrieved: {len(archive_result)}')
    else:
        print(response.json()['errorMessage'])
        print('Returning empty table...')
    return archive_result


def fix_table_obj_to_float(tablet):
    '''
    Makes a Table which is a subset of the original table and changes the dtype of the columns that are
    object to float types. The Table returned can be used to update the original Table to change the dtypes.

    Parameters:
    -----------
    table: Table
        Astropy Table that has data Columns with Object dtype and can be turned into floats.

    Returns:
    --------
    subset_table: Table
        Astropy Table containing only the changed columns. The user can run tablet.update(subset_table) to
        update the original Table.
    '''

    ign_warn = warnings.filterwarnings('ignore','RuntimeWarning')

    tmp_table = Table()
    for column in tablet.columns:
        if tablet[column].dtype==object:
            try:
                tmp_table.add_column(Column(tablet[column], name=column, dtype=float))
            except:
                tmp_table.add_column(Column(tablet[column], name=column, dtype=str))
    return tmp_table


if __name__=="__main__":
    print('This is the main program.')
