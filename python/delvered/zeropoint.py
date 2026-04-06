"""
Photometric zero-point calibration for the DELVERED pipeline.

Python port of:
  - delvered_zeropoint.pro           → run_zeropoint()
  - delvered_zeropoint_exposure.pro  → zeropoint_exposure()
"""

import os
import re
import glob
import numpy as np
from pathlib import Path
from astropy.io import fits
from astropy.coordinates import SkyCoord, match_coordinates_sky
import astropy.units as u
from scipy.optimize import least_squares


# ---------------------------------------------------------------------------
# Zero-point for a single exposure
# ---------------------------------------------------------------------------

def zeropoint_exposure(expname, modeleqnfile, refcat_dir='refcat',
                       apcor_file='apcor.lst', redo=False):
    """
    Calculate the photometric zero-point for a single exposure.

    Python port of delvered_zeropoint_exposure.pro.

    Loads all ALS files for the exposure, cross-matches to a reference
    catalogue, computes model magnitudes, and fits the zero-point offset.

    Parameters
    ----------
    expname : str
        Exposure name, e.g. ``'F1-00912345'``.
    modeleqnfile : str
        Path to the model magnitude equation file.
    refcat_dir : str
        Directory containing field reference catalogue FITS files.
    apcor_file : str
        Path to aperture correction list file.
    redo : bool
        If True, redo even if output file already exists.

    Returns
    -------
    dict or None
        Exposure zero-point structure with keys: name, field, filter,
        exptime, ncat, nref, num, zpterm, zptermerr, translines.
        Returns None on failure.
    """
    from .photometry_utils import get_model_mag

    field = expname.split('-')[0]
    outfile = os.path.join(field, f'{expname}_zeropoint.fits')

    if os.path.exists(outfile) and not redo:
        print(f"zeropoint_exposure: {outfile} EXISTS and redo=False")
        return fits.getdata(outfile, 1)

    expstr = {
        'name': expname, 'field': field, 'filter': '', 'exptime': 0.0,
        'ncat': 0, 'nref': 0, 'num': 0, 'zpterm': 99.99, 'zptermerr': 9.99,
        'translines': ['', '', ''],
    }

    # Load aperture corrections
    apcor = {}
    if os.path.exists(apcor_file):
        with open(apcor_file) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    base = os.path.basename(parts[0]).replace('a.del', '')
                    try:
                        apcor[base] = float(parts[1])
                    except ValueError:
                        pass

    # Find ALS files for this exposure
    als_pattern = os.path.join(field, 'chip??', f'{expname}_??.als')
    als_files = sorted(glob.glob(als_pattern))
    fits_pattern = os.path.join(field, 'chip??', f'{expname}_??.fits*')
    fits_files = sorted(glob.glob(fits_pattern))

    if not als_files:
        print(f"zeropoint_exposure: no ALS files found for {expname}")
        return None

    # Get filter and exptime from first FITS header
    fitsfile0 = fits_files[0] if fits_files else None
    filter_name = ''
    exptime = 1.0
    if fitsfile0:
        try:
            hdr = fits.getheader(fitsfile0)
            filter_name = str(hdr.get('FILTER', '')).strip()
            # Normalise filter name (e.g. 'g DECam SDSS c0002...' → 'g')
            filter_name = filter_name.split()[0].lower()
            exptime = float(hdr.get('EXPTIME', 1.0))
        except Exception:
            pass

    expstr['filter'] = filter_name
    expstr['exptime'] = exptime

    # Load and concatenate all chip ALS catalogs
    all_rows = []
    for als_file in als_files:
        base = os.path.splitext(os.path.basename(als_file))[0]
        chip_match = re.search(r'_(\d+)$', base)
        ccdnum = int(chip_match.group(1)) if chip_match else 0

        # Find matching FITS file for WCS
        fits_file = als_file.replace('.als', '.fits')
        if not os.path.exists(fits_file):
            fits_file = als_file.replace('.als', '.fits.fz')
        if not os.path.exists(fits_file):
            continue

        als_cat = _load_als(als_file)
        if als_cat is None or len(als_cat) == 0:
            continue

        # WCS: pixel → RA/DEC
        try:
            from astropy.wcs import WCS
            with fits.open(fits_file) as hdul:
                hdr = hdul[1].header if fits_file.endswith('.fz') else hdul[0].header
            wcs = WCS(hdr)
            sky = wcs.pixel_to_world(als_cat['x'] - 1, als_cat['y'] - 1)
            ra = sky.ra.deg
            dec = sky.dec.deg
        except Exception:
            ra = np.zeros(len(als_cat))
            dec = np.zeros(len(als_cat))

        # Calibrated instrumental magnitude
        apcor_val = apcor.get(base, 0.0)
        cmag = als_cat['mag'] + 2.5 * np.log10(exptime) - apcor_val

        for i in range(len(als_cat)):
            all_rows.append({
                'id': f'{field}_{ccdnum}.{als_cat["id"][i]}',
                'ra': ra[i], 'dec': dec[i],
                'mag': als_cat['mag'][i], 'err': als_cat['err'][i],
                'cmag': cmag[i], 'chi': als_cat['chi'][i],
                'sharp': als_cat['sharp'][i], 'sky': als_cat['sky'][i],
            })

    if not all_rows:
        print(f"zeropoint_exposure: no ALS measurements for {expname}")
        return None

    # Build catalogue array
    cat = np.array(all_rows,
                   dtype=[('id', 'U50'), ('ra', float), ('dec', float),
                          ('mag', float), ('err', float), ('cmag', float),
                          ('chi', float), ('sharp', float), ('sky', float)])
    expstr['ncat'] = len(cat)

    # Determine field centre
    ra_arr = cat['ra']
    if np.ptp(ra_arr) > 100:
        ra_adj = np.where(ra_arr > 180, ra_arr - 360, ra_arr)
        cenra = float(np.mean([ra_adj.min(), ra_adj.max()]))
        if cenra < 0:
            cenra += 360
    else:
        cenra = float(np.mean([ra_arr.min(), ra_arr.max()]))
    cendec = float(np.mean([cat['dec'].min(), cat['dec'].max()]))

    # Load reference catalogue
    reffile = os.path.join(refcat_dir, f'{field}_refcat.fits.gz')
    if not os.path.exists(reffile):
        reffile = os.path.join(refcat_dir, f'{field}_refcat.fits')
    if not os.path.exists(reffile):
        print(f"zeropoint_exposure: {reffile} NOT FOUND")
        return None

    ref = fits.getdata(reffile, 1)
    expstr['nref'] = len(ref)

    # Cross-match catalogue to reference (1 arcsec radius)
    cat_coord = SkyCoord(ra=cat['ra'] * u.deg, dec=cat['dec'] * u.deg)
    ref_coord = SkyCoord(ra=ref['ra'].astype(float) * u.deg,
                         dec=ref['dec'].astype(float) * u.deg)
    idx_ref, sep, _ = match_coordinates_sky(cat_coord, ref_coord)
    match_mask = sep.arcsec < 1.0
    if not match_mask.any():
        print(f"zeropoint_exposure: no cross-matches for {expname}")
        return None

    cat1 = cat[match_mask]
    ref1 = ref[idx_ref[match_mask]]

    # Get model magnitudes
    instfilt = f'c4d-{filter_name}'
    mmags = get_model_mag(ref1, instfilt, cendec, modeleqnfile)

    # Quality cuts
    good = ((cat1['mag'] < 50) & (cat1['err'] < 0.05) &
            (np.abs(cat1['sharp']) < 1) & (cat1['chi'] < 3) &
            (mmags[:, 0] < 50) & (mmags[:, 1] < 5))
    if good.sum() < 10:
        good = ((cat1['mag'] < 50) & (cat1['err'] < 0.08) &
                (np.abs(cat1['sharp']) < 1) & (cat1['chi'] < 3) &
                (mmags[:, 0] < 50) & (mmags[:, 1] < 5))
    if not good.any():
        print(f"zeropoint_exposure: no stars pass quality cuts for {expname}")
        return None

    diff = mmags[good, 0] - cat1['cmag'][good]
    err = np.sqrt(mmags[good, 1]**2 + cat1['err'][good]**2)

    # 3-sigma clip
    med = np.median(diff)
    sig = 1.4826 * np.median(np.abs(diff - med))
    keep = np.abs(diff - med) < 3 * sig
    if keep.sum() < 3:
        keep = np.ones(len(diff), dtype=bool)

    zpterm, zptermerr = _fit_zpterm(diff[keep], err[keep])

    expstr['num'] = int(keep.sum())
    expstr['zpterm'] = float(zpterm)
    expstr['zptermerr'] = float(zptermerr)
    expstr['translines'] = [
        f"{expname}  {filter_name}  {filter_name}-{filter_name}  "
        f"{-zpterm:7.4f}    0.0000    0.0000   0.0000   0.0000",
        f"                     {zptermerr:7.4f}    0.0000    0.0000   0.0000   0.0000",
        '',
    ]

    print(f"zeropoint_exposure: {expname} ZPTERM={zpterm:.4f} +/- {zptermerr:.4f} "
          f"(n={keep.sum()})")

    # Save output
    os.makedirs(field, exist_ok=True)
    _write_expstr(expstr, outfile)

    return expstr


def _fit_zpterm(diff, err):
    """Fit a constant zero-point offset using weighted least squares."""
    wt = 1.0 / (err**2)
    zpterm = np.sum(wt * diff) / np.sum(wt)
    zptermerr = np.sqrt(1.0 / np.sum(wt))
    return float(zpterm), float(zptermerr)


def _write_expstr(expstr, outfile):
    """Write exposure zero-point structure to a FITS file."""
    dt = np.dtype([
        ('name', 'U50'), ('field', 'U20'), ('filter', 'U4'),
        ('exptime', float), ('ncat', int), ('nref', int), ('num', int),
        ('zpterm', float), ('zptermerr', float),
        ('transline0', 'U120'), ('transline1', 'U120'), ('transline2', 'U120'),
    ])
    row = np.zeros(1, dtype=dt)
    row['name'] = expstr['name']
    row['field'] = expstr['field']
    row['filter'] = expstr['filter']
    row['exptime'] = expstr['exptime']
    row['ncat'] = expstr['ncat']
    row['nref'] = expstr['nref']
    row['num'] = expstr['num']
    row['zpterm'] = expstr['zpterm']
    row['zptermerr'] = expstr['zptermerr']
    tlines = expstr.get('translines', ['', '', ''])
    row['transline0'] = tlines[0]
    row['transline1'] = tlines[1]
    row['transline2'] = tlines[2]
    hdu = fits.BinTableHDU(row)
    hdu.writeto(outfile, overwrite=True)


# ---------------------------------------------------------------------------
# Loop over all exposures
# ---------------------------------------------------------------------------

def run_zeropoint(modeleqnfile=None, refcat_dir='refcat',
                  apcor_file='apcor.lst', nmulti=5, redo=False):
    """
    Calculate zero-points for all exposures in the current directory.

    Python port of delvered_zeropoint.pro.

    Scans for ALS files from DAOPHOT processing, groups them by exposure,
    and calls ``zeropoint_exposure`` for each.  Results are written to
    ``{field}/{expname}_zeropoint.fits`` and collated into
    ``delve.trans``.

    Parameters
    ----------
    modeleqnfile : str
        Path to the model magnitude equation file.
    refcat_dir : str
        Directory with field reference catalogues.
    apcor_file : str
        Path to aperture correction list.
    nmulti : int
        Number of parallel jobs (currently sequential; use external
        scheduler for true parallelism).
    redo : bool
        Redo previously processed exposures.
    """
    from concurrent.futures import ProcessPoolExecutor, as_completed

    if modeleqnfile is None:
        print("run_zeropoint: modeleqnfile is required")
        return

    # Discover ALS files
    als_files = sorted(glob.glob('*/chip??/*.als'))
    if not als_files:
        print("run_zeropoint: no ALS files found")
        return

    # Group by exposure name (strip chip suffix)
    exp_to_files = {}
    for f in als_files:
        base = os.path.splitext(os.path.basename(f))[0]
        exp = '_'.join(base.split('_')[:-1])
        exp_to_files.setdefault(exp, []).append(f)

    print(f"run_zeropoint: {len(exp_to_files)} unique exposures")

    # Process exposures
    expstrs = []
    for expname in sorted(exp_to_files):
        es = zeropoint_exposure(
            expname, modeleqnfile,
            refcat_dir=refcat_dir, apcor_file=apcor_file, redo=redo)
        if es is not None:
            expstrs.append(es)

    if not expstrs:
        print("run_zeropoint: no successful zero-points")
        return

    # Fill in failed exposures with band-average
    _fill_failed_zp(expstrs)

    # Write delve.trans
    trans_lines = []
    if os.path.exists('delve.trans') and not redo:
        with open('delve.trans') as f:
            trans_lines = f.readlines()
    for es in expstrs:
        if es['num'] > 0:
            trans_lines.extend([l + '\n' for l in es['translines'] if l])
    with open('delve.trans', 'w') as f:
        f.writelines(trans_lines)
    print(f"run_zeropoint: wrote {len(expstrs)} zero-points to delve.trans")


def _fill_failed_zp(expstrs):
    """Replace failed zero-points with the filter-average."""
    filters = list(set(es['filter'] for es in expstrs))
    for filt in filters:
        good = [es for es in expstrs if es['filter'] == filt and es['num'] > 0]
        bad = [es for es in expstrs if es['filter'] == filt and es['num'] == 0]
        if not good or not bad:
            continue
        mean_zp = np.mean([es['zpterm'] for es in good])
        mean_zperr = np.mean([es['zptermerr'] for es in good])
        for es in bad:
            es['zpterm'] = float(mean_zp)
            es['zptermerr'] = float(mean_zperr)
            es['num'] = 1
            expname = es['name']
            filt = es['filter']
            es['translines'] = [
                f"{expname}  {filt}  {filt}-{filt}  "
                f"{-mean_zp:7.4f}    0.0000    0.0000   0.0000   0.0000",
                f"                     {mean_zperr:7.4f}    0.0000    0.0000   0.0000   0.0000",
                '',
            ]
        print(f"run_zeropoint: filled {len(bad)} failed {filt}-band exposures "
              f"with mean ZP={mean_zp:.4f}")


# ---------------------------------------------------------------------------
# ALS file parser
# ---------------------------------------------------------------------------

def _load_als(als_file):
    """
    Parse a DAOPHOT ALS file into a numpy structured array.

    ALS format (fixed-width):
      ID  X  Y  MAG  ERR  SKY  ITER  CHI  SHARP
    """
    rows = []
    try:
        with open(als_file) as f:
            lines = f.readlines()
    except IOError:
        return None

    # Skip header lines (those that don't start with a number)
    for line in lines:
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) < 9:
            continue
        try:
            rows.append({
                'id': parts[0],
                'x': float(parts[1]), 'y': float(parts[2]),
                'mag': float(parts[3]), 'err': float(parts[4]),
                'sky': float(parts[5]),
                'iter': float(parts[6]),
                'chi': float(parts[7]), 'sharp': float(parts[8]),
            })
        except (ValueError, IndexError):
            continue

    if not rows:
        return None

    dt = [('id', 'U20'), ('x', float), ('y', float), ('mag', float),
          ('err', float), ('sky', float), ('iter', float),
          ('chi', float), ('sharp', float)]
    return np.array(rows, dtype=dt)
