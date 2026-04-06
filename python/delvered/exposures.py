"""
DECam exposure list construction for the DELVERED pipeline.

Python port of:
  - make_exposures_list.pro      → make_exposures_list()
  - make_exposures_list_deep.pro → make_exposures_list_deep()
"""

import gzip
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u


# ---------------------------------------------------------------------------
# LMC / SMC centre coordinates (J2000)
# ---------------------------------------------------------------------------

_LMC_RA = 81.9      # degrees
_LMC_DEC = -69.867
_SMC_RA = 13.183
_SMC_DEC = -72.8283

# DELVE / MagLiteS proposal IDs accepted regardless of release date
_ALLOWED_PROPIDS = ('2019A-0305', '2018A-0242')


def _angular_separation(ra1, dec1, ra2, dec2):
    """Return angular separation in degrees between (ra1,dec1) and (ra2,dec2)."""
    c1 = SkyCoord(ra=ra1 * u.deg, dec=dec1 * u.deg, frame='icrs')
    c2 = SkyCoord(ra=ra2 * u.deg, dec=dec2 * u.deg, frame='icrs')
    return c1.separation(c2).deg


def _load_fits_gz(path):
    """Load the first binary table extension from a (optionally gzipped) FITS file."""
    path = str(path)
    if path.endswith('.gz'):
        with gzip.open(path) as gz:
            return fits.getdata(gz, 1)
    return fits.getdata(path, 1)


def _apply_release_cuts(str_arr):
    """
    Keep rows that are either public (release_date <= today) or belong to
    DELVE (2019A-0305) / MagLiteS (2018A-0242) proposals.
    """
    release_date = np.char.strip(str_arr['release_date'].astype(str))
    prop_id = np.char.strip(str_arr['prop_id'].astype(str))

    today = datetime.now(tz=timezone.utc).strftime('%Y-%m-%d')

    # Compare ISO date strings lexicographically (YYYY-MM-DD format)
    public_mask = release_date <= today
    delve_mask = np.isin(prop_id, list(_ALLOWED_PROPIDS))
    return str_arr[public_mask | delve_mask]


# ---------------------------------------------------------------------------
# Main routines
# ---------------------------------------------------------------------------

def make_exposures_list(input_file, output_dir=None, delvereddir=None):
    """
    Build a DECam Magellanic Cloud exposure list from the NSC master DECam list.

    Python port of make_exposures_list.pro.

    Loads the master DECam instcal FITS list, applies release-date cuts
    (public data up to today, or DELVE / MagLiteS proposals), then keeps
    only exposures with |b| > 10° that fall within 25° of the LMC or 15°
    of the SMC.

    Parameters
    ----------
    input_file : str or Path
        Path to the master DECam exposure list FITS file
        (e.g. ``delvemc_info_*.fits.gz``).
    output_dir : str or Path, optional
        Directory for the output file.  Defaults to
        ``{delvereddir}/data/`` when *delvereddir* is provided, otherwise
        the current directory.
    delvereddir : str or Path, optional
        Root of the delvered project tree.

    Returns
    -------
    numpy structured array
        Filtered exposure table.
    str
        Path to the written output file.
    """
    input_file = Path(input_file)
    print(f"make_exposures_list: loading {input_file}")
    str_arr = _load_fits_gz(input_file)

    # Normalise prop_id whitespace
    str_arr['prop_id'][:] = np.char.strip(str_arr['prop_id'].astype(str))

    # Apply release-date / proposal cuts
    str_arr = _apply_release_cuts(str_arr)
    print(f"make_exposures_list: {len(str_arr)} total public DECam or DELVE/MagLiteS exposures")

    # Galactic coordinates
    coords = SkyCoord(ra=str_arr['ra'] * u.deg,
                      dec=str_arr['dec'] * u.deg, frame='icrs')
    gal = coords.galactic
    glat = gal.b.deg

    # LMC / SMC angular distances
    lmc_rad = _angular_separation(str_arr['ra'], str_arr['dec'], _LMC_RA, _LMC_DEC)
    smc_rad = _angular_separation(str_arr['ra'], str_arr['dec'], _SMC_RA, _SMC_DEC)

    mc_mask = (np.abs(glat) > 10) & ((lmc_rad < 25) | (smc_rad < 15))
    mc = str_arr[mc_mask]
    print(f"make_exposures_list: {len(mc)} MC exposures")

    # Write output
    now = datetime.now(tz=timezone.utc)
    tag = now.strftime('%Y%m%d')
    if output_dir is None:
        if delvereddir is not None:
            output_dir = Path(delvereddir) / 'data'
        else:
            output_dir = Path('.')
    output_dir = Path(output_dir)
    outfile = output_dir / f'decam_mcs_{tag}.fits'
    hdu = fits.BinTableHDU(mc)
    hdu.writeto(str(outfile), overwrite=True)
    print(f"make_exposures_list: writing to {outfile}")
    subprocess.run(['gzip', '-f', str(outfile)], check=False)
    outfile_gz = Path(str(outfile) + '.gz')
    print(f"make_exposures_list: wrote {outfile_gz}")
    return mc, str(outfile_gz)


def make_exposures_list_deep(input_files, output_dir=None, delvereddir=None):
    """
    Build a DECam deep-field exposure list from per-field DECam instcal lists.

    Python port of make_exposures_list_deep.pro.

    Loads one or more per-field DECam FITS lists (e.g. SexB, NGC 55),
    concatenates them, applies release-date / proposal cuts, and writes
    the combined list.

    Parameters
    ----------
    input_files : list of str or Path
        Paths to per-field DECam exposure FITS files
        (e.g. ``decam_sexb_*.fits.gz``).
    output_dir : str or Path, optional
        Directory for the output file.  Defaults to the current directory.
    delvereddir : str or Path, optional
        Root of the delvered project tree.

    Returns
    -------
    numpy structured array
        Filtered exposure table.
    str
        Path to the written output file.
    """
    from numpy.lib.recfunctions import stack_arrays

    parts = []
    for f in input_files:
        f = Path(f)
        print(f"make_exposures_list_deep: loading {f}")
        arr = _load_fits_gz(f)
        parts.append(arr)

    if not parts:
        raise ValueError("make_exposures_list_deep: no input files loaded")

    str_arr = stack_arrays(parts, asrecarray=False, usemask=False, autoconvert=True)
    str_arr['prop_id'][:] = np.char.strip(str_arr['prop_id'].astype(str))

    str_arr = _apply_release_cuts(str_arr)
    print(f"make_exposures_list_deep: {len(str_arr)} total public DECam or DELVE/MagLiteS exposures")

    # Write output
    now = datetime.now(tz=timezone.utc)
    tag = now.strftime('%Y%m%d')
    if output_dir is None:
        if delvereddir is not None:
            output_dir = Path(delvereddir) / 'data'
        else:
            output_dir = Path('.')
    output_dir = Path(output_dir)
    outfile = output_dir / f'decam_deep_{tag}.fits'
    hdu = fits.BinTableHDU(str_arr)
    hdu.writeto(str(outfile), overwrite=True)
    print(f"make_exposures_list_deep: writing to {outfile}")
    subprocess.run(['gzip', '-f', str(outfile)], check=False)
    outfile_gz = Path(str(outfile) + '.gz')
    print(f"make_exposures_list_deep: wrote {outfile_gz}")
    return str_arr, str(outfile_gz)
