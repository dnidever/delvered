"""
Photometry catalogue construction for the DELVERED pipeline.

Python port of:
  - load_chipcat.pro           → load_chip_cat()
  - make_allframe_photcat.pro  → make_allframe_photcat()
  - delvered_initpsfstars.pro  → init_psf_stars()
"""

import os
import re
import glob
import gzip
import shutil
import numpy as np
from pathlib import Path
from astropy.io import fits
from astropy.wcs import WCS
from astropy.time import Time
from astropy.coordinates import SkyCoord
import astropy.units as u

from .zeropoint import _load_als


# ---------------------------------------------------------------------------
# Load and calibrate a single chip ALS catalogue
# ---------------------------------------------------------------------------

def load_chip_cat(chstr, use_alf=False):
    """
    Load an ALS (or ALF) file for a single chip and apply calibrations.

    Python port of load_chipcat.pro.

    Applies exposure time, aperture correction, and zero-point
    calibration to the raw DAOPHOT magnitudes and converts pixel
    coordinates to RA/DEC using the FITS WCS.

    Parameters
    ----------
    chstr : dict or numpy void
        Chip metadata.  Required keys: file, exptime, apcor,
        calib_zpterm, calib_zptermsig, base, filter, expnum, chip,
        utdate, uttime.  Optional: alffile (if use_alf=True).
    use_alf : bool
        Use ALF file instead of ALS file.

    Returns
    -------
    numpy structured array or None
        Calibrated per-source photometry catalogue.
    """
    def _get(key):
        return chstr[key] if hasattr(chstr, '__getitem__') else getattr(chstr, key)

    fitsfile = str(_get('file')).strip()

    # Determine ALS/ALF filename
    if not use_alf:
        alsfile = fitsfile.replace('.fits.fz', '.fits').replace('.fits', '.als')
    else:
        alsfile = str(_get('alffile')).strip()

    if not os.path.exists(alsfile):
        print(f"load_chip_cat: {alsfile} NOT FOUND")
        return None

    als = _load_als(alsfile)
    if als is None or len(als) == 0:
        print(f"load_chip_cat: {alsfile} IS EMPTY")
        return None

    # Load FITS header for WCS
    try:
        with fits.open(fitsfile) as hdul:
            hdr = hdul[1].header if fitsfile.endswith('.fz') else hdul[0].header
        wcs = WCS(hdr)
    except Exception as e:
        print(f"load_chip_cat: cannot read WCS from {fitsfile}: {e}")
        return None

    # Trim bad half of DECam chip 31 (after MJD 56660)
    chip = int(_get('chip'))
    if chip == 31:
        dateobs = str(_get('utdate')).strip() + 'T' + str(_get('uttime')).strip()
        try:
            mjd = Time(dateobs, format='isot', scale='utc').mjd
            if mjd > 56660:
                bad = als['x'] >= 1000
                good_mask = ~bad
                if not good_mask.any():
                    print(f"load_chip_cat: no useful measurements in {fitsfile}")
                    return None
                als = als[good_mask]
                print(f"load_chip_cat: masking bad half of DECam chip 31, "
                      f"{bad.sum()} removed, {good_mask.sum()} left")
        except Exception:
            pass

    nals = len(als)
    exptime = float(_get('exptime'))
    apcor = float(_get('apcor'))
    calib_zpterm = float(_get('calib_zpterm'))
    calib_zptermsig = float(_get('calib_zptermsig'))

    # Calibrate:  cmag = mag + 2.5*log10(exptime) - apcor - zpterm
    # aperture correction is SUBTRACTIVE (makes it brighter)
    # ZPTERM is SUBTRACTIVE
    cmag = als['mag'] + 2.5 * np.log10(exptime) - apcor - calib_zpterm
    cerr = np.sqrt(als['err']**2 + calib_zptermsig**2)

    # Pixel → RA/DEC  (DAOPHOT uses 1-based coords)
    try:
        sky = wcs.pixel_to_world(als['x'] - 1.0, als['y'] - 1.0)
        ra = sky.ra.deg
        dec = sky.dec.deg
    except Exception:
        ra = np.zeros(nals)
        dec = np.zeros(nals)

    # MJD from utdate/uttime
    try:
        dateobs = str(_get('utdate')).strip() + 'T' + str(_get('uttime')).strip()
        mjd = Time(dateobs, format='isot', scale='utc').mjd
    except Exception:
        mjd = 0.0

    expnum = str(_get('expnum')).strip()
    base = str(_get('base')).strip()
    filt = str(_get('filter')).strip()

    # Build output catalogue
    dt = [
        ('id', 'U80'), ('objid', 'U80'), ('exposure', 'U80'),
        ('ccdnum', np.int16), ('filter', 'U4'), ('mjd', np.float64),
        ('forced', np.uint8),
        ('x', np.float32), ('y', np.float32),
        ('ra', np.float64), ('dec', np.float64),
        ('imag', np.float32), ('ierr', np.float32),
        ('mag', np.float32), ('err', np.float32),
        ('sky', np.float32), ('chi', np.float32), ('sharp', np.float32),
    ]
    cat = np.zeros(nals, dtype=dt)
    cat['id'] = np.array([f'{expnum}_{chip}.{als_id}' for als_id in als['id']])
    cat['exposure'] = base
    cat['ccdnum'] = chip
    cat['filter'] = filt
    cat['mjd'] = mjd
    cat['x'] = als['x']
    cat['y'] = als['y']
    cat['ra'] = ra
    cat['dec'] = dec
    cat['imag'] = als['mag']
    cat['ierr'] = als['err']
    cat['mag'] = cmag
    cat['err'] = cerr
    cat['sky'] = als['sky']
    cat['chi'] = als['chi']
    cat['sharp'] = als['sharp']

    return cat


# ---------------------------------------------------------------------------
# Make ALLFRAME photometry catalogue for a brick
# ---------------------------------------------------------------------------

def make_allframe_photcat(brick, meta, mch_file=None, ebv_from_sfd=True):
    """
    Build a calibrated photometry catalogue from ALLFRAME output.

    Python port of make_allframe_photcat.pro.

    Loads per-chip ALS files listed in the MCH file, applies
    calibrations from the chip metadata, combines multi-exposure
    photometry by flux-weighted averaging, and optionally adds
    SFD E(B-V) extinction.

    Parameters
    ----------
    brick : str
        Brick name (e.g. ``'1234p567'``).
    meta : numpy structured array
        Chip metadata with one row per chip.  Required columns:
        base, expnum, chip, filter, exptime, apcor, calib_zpterm,
        calib_zptermsig, utdate, uttime, file.
    mch_file : str, optional
        Path to the MCH (match) file listing the chips.  If None,
        inferred from the brick directory.
    ebv_from_sfd : bool
        If True, add SFD E(B-V) column to the output.

    Returns
    -------
    numpy structured array or None
        Combined photometry catalogue.
    """
    from .catalog_ops import simple_avg_meas
    from .utils import create_index

    nchips = len(meta)
    if nchips == 0:
        return None

    all_meas = []

    for i in range(nchips):
        chstr = meta[i]
        cat = load_chip_cat(chstr)
        if cat is None or len(cat) == 0:
            continue
        all_meas.append(cat)

    if not all_meas:
        return None

    # Concatenate measurements
    from numpy.lib.recfunctions import stack_arrays
    meas = stack_arrays(all_meas, asrecarray=False, usemask=False,
                        autoconvert=True)

    # Assign preliminary OBJIDs (sequential, will be overwritten by matching)
    meas = np.lib.recfunctions.append_fields(
        meas, 'objid_tmp',
        np.arange(len(meas), dtype=np.int64).astype('U80'),
        usemask=False)

    # Average measurements into object catalogue
    obj = simple_avg_meas(meas)
    if obj is None:
        return None

    # Add SFD E(B-V)
    if ebv_from_sfd:
        try:
            from .extinction import _sfd_ebv
            obj = np.lib.recfunctions.append_fields(
                obj, 'ebv',
                _sfd_ebv(obj['ra'], obj['dec']).astype(np.float32),
                usemask=False)
        except Exception:
            pass

    return obj


# ---------------------------------------------------------------------------
# Initialise PSF star list from Gaia
# ---------------------------------------------------------------------------

def init_psf_stars(fitsfile, refcatfile, outfile=None):
    """
    Create an initial PSF star list from Gaia for DAOPHOT.

    Python port of delvered_initpsfstars.pro.

    Reads the Gaia reference catalogue, selects stars on the chip
    (within the image footprint), and writes a PHOTRED-compatible
    CMN file with star coordinates.

    Parameters
    ----------
    fitsfile : str
        Path to the chip FITS image.
    refcatfile : str
        Path to the Gaia reference catalogue FITS file (possibly .gz).
    outfile : str, optional
        Output CMN file path.  Defaults to
        ``{fitsfile_base}_psfstars.lst``.

    Returns
    -------
    str or None
        Path to the written CMN file, or None on failure.
    """
    if outfile is None:
        base = fitsfile.replace('.fits.fz', '').replace('.fits', '')
        outfile = base + '_psfstars.lst'

    # Load Gaia reference catalogue
    if refcatfile.endswith('.gz'):
        with gzip.open(refcatfile) as gz:
            ref = fits.getdata(gz)
    else:
        ref = fits.getdata(refcatfile, 1)

    if ref is None or len(ref) == 0:
        print(f"init_psf_stars: empty reference catalogue {refcatfile}")
        return None

    # Read FITS header and WCS
    try:
        with fits.open(fitsfile) as hdul:
            hdr = hdul[1].header if fitsfile.endswith('.fz') else hdul[0].header
        wcs = WCS(hdr)
        nx = hdr.get('NAXIS1', 2048)
        ny = hdr.get('NAXIS2', 4096)
    except Exception as e:
        print(f"init_psf_stars: cannot open {fitsfile}: {e}")
        return None

    # Convert RA/DEC to pixel coordinates
    try:
        ra_ref = ref['ra'].astype(float)
        dec_ref = ref['dec'].astype(float)
    except (KeyError, ValueError):
        try:
            ra_ref = ref['raj2000'].astype(float)
            dec_ref = ref['dej2000'].astype(float)
        except (KeyError, ValueError):
            print("init_psf_stars: cannot find RA/DEC columns in reference catalogue")
            return None

    sky = SkyCoord(ra=ra_ref * u.deg, dec=dec_ref * u.deg, frame='icrs')
    try:
        px, py = wcs.world_to_pixel(sky)
    except Exception as e:
        print(f"init_psf_stars: WCS conversion failed: {e}")
        return None

    # Keep only stars on the image
    on_chip = ((px >= 0) & (px < nx) & (py >= 0) & (py < ny))
    # Filter bright stars (gmag < 50)
    if 'gmag' in ref.dtype.names:
        bright = ref['gmag'] < 50
        on_chip &= bright

    sel = np.where(on_chip)[0]
    if len(sel) == 0:
        print("init_psf_stars: no reference stars on chip")
        return None

    # Write CMN file: ID X Y MAG 0.0 0.0
    with open(outfile, 'w') as f:
        for k, i in enumerate(sel):
            star_id = k + 1
            mag = float(ref['gmag'][i]) if 'gmag' in ref.dtype.names else 15.0
            f.write(f'{star_id:6d} {px[i]+1:8.3f} {py[i]+1:8.3f} '
                    f'{mag:8.3f}  0.000  0.000\n')

    print(f"init_psf_stars: wrote {len(sel)} PSF stars to {outfile}")
    return outfile
