"""
WCS-fitting wrapper for the DELVERED pipeline.

Python port of delvered_wcs.pro → run_wcs()
"""

import re
import subprocess
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u


def _read_setup(setup_file):
    cfg = {}
    p = Path(setup_file)
    if not p.exists():
        return cfg
    with open(p) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split(None, 1)
            if len(parts) == 2:
                cfg[parts[0].upper()] = parts[1].strip()
    return cfg


def _read_list(path):
    p = Path(path)
    if not p.exists():
        return []
    with open(p) as fh:
        return [ln.strip() for ln in fh if ln.strip()]


def _write_list(path, lines):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'a') as fh:
        fh.write('\n'.join(lines) + '\n')


def _has_wcs(fits_file):
    """Return True if the FITS file has a WCS-fitted header (WCSFIT HISTORY)."""
    try:
        if str(fits_file).endswith('.fz'):
            hdr = fits.getheader(str(fits_file), ext=1)
        else:
            hdr = fits.getheader(str(fits_file), ext=0)
        ctype1 = str(hdr.get('CTYPE1', '')).strip()
        if not ctype1 or ctype1 == '0':
            return False
        for card in hdr.get('HISTORY', []):
            if 'WCSFIT: RMS' in str(card):
                return True
        return False
    except Exception:
        return False


def _fit_wcs_astrometry_net(fits_file, refcat_file=None,
                             search_radius=0.5, timeout=300):
    """
    Fit WCS using astrometry.net ``solve-field``.

    Parameters
    ----------
    fits_file : str or Path
        Input chip FITS file.
    refcat_file : str or Path, optional
        Reference catalogue FITS for local index (unused if not provided).
    search_radius : float
        Search radius in degrees for astrometry.net.
    timeout : int
        Subprocess timeout in seconds.

    Returns
    -------
    bool
        True if solved successfully.
    """
    fits_file = Path(fits_file)
    try:
        hdr = fits.getheader(str(fits_file))
        ra = float(hdr.get('RA', hdr.get('CRVAL1', 0)))
        dec = float(hdr.get('DEC', hdr.get('CRVAL2', 0)))
    except Exception:
        ra, dec = 0.0, 0.0

    cmd = [
        'solve-field',
        '--no-plots', '--overwrite',
        '--ra', str(ra), '--dec', str(dec),
        '--radius', str(search_radius),
        str(fits_file),
    ]
    try:
        result = subprocess.run(cmd, capture_output=True,
                                timeout=timeout, check=False)
        # solve-field writes a .new file on success
        solved = fits_file.with_suffix('.new').exists()
        if solved:
            # Copy the solved header back
            new_hdr = fits.getheader(
                str(fits_file.with_suffix('.new')), ext=0)
            with fits.open(str(fits_file), mode='update') as hdul:
                for key in ('CTYPE1', 'CTYPE2', 'CRVAL1', 'CRVAL2',
                            'CRPIX1', 'CRPIX2', 'CD1_1', 'CD1_2',
                            'CD2_1', 'CD2_2', 'CDELT1', 'CDELT2'):
                    if key in new_hdr:
                        hdul[0].header[key] = new_hdr[key]
                hdul[0].header.add_history(
                    'WCSFIT: RMS = solved by astrometry.net')
                hdul.flush()
        return solved
    except Exception as exc:
        print(f"  solve-field failed for {fits_file.name}: {exc}")
        return False


def _fit_wcs_scamp(fits_file, refcat_file, timeout=300):
    """
    Fit WCS using SCAMP cross-matching against a reference catalogue.

    Parameters
    ----------
    fits_file : str or Path
        Input chip FITS file.
    refcat_file : str or Path
        Reference catalogue FITS file (e.g. Gaia DR2).
    timeout : int
        Subprocess timeout in seconds.

    Returns
    -------
    bool
        True if SCAMP produced a WCS header.
    """
    fits_file = Path(fits_file)
    refcat_file = Path(refcat_file)
    if not refcat_file.exists():
        return False

    head_file = fits_file.with_suffix('.head')
    cmd = [
        'scamp', str(fits_file),
        '-ASTREF_CATALOG', 'FILE',
        '-ASTREFCAT_NAME', str(refcat_file),
        '-SAVE_REFCATALOG', 'N',
        '-VERBOSE_TYPE', 'QUIET',
        '-WRITE_XML', 'N',
    ]
    try:
        subprocess.run(cmd, capture_output=True,
                       timeout=timeout, check=False)
        if head_file.exists():
            with fits.open(str(fits_file), mode='update') as hdul:
                hdul[0].header.add_history('WCSFIT: RMS = solved by SCAMP')
                hdul.flush()
            return True
    except Exception as exc:
        print(f"  SCAMP failed for {fits_file.name}: {exc}")
    return False


def _run_one_wcs(args):
    """
    Fit WCS for a single chip FITS file.

    Tries astrometry.net first, then SCAMP.

    Parameters
    ----------
    args : tuple
        (fits_file, refcat_file, search_radius, redo)

    Returns
    -------
    (fits_file_str, success: bool)
    """
    fits_file, refcat_file, search_radius, redo = args
    fits_file = Path(fits_file)

    if not redo and _has_wcs(fits_file):
        return str(fits_file), True

    # Try astrometry.net
    ok = _fit_wcs_astrometry_net(fits_file, search_radius=search_radius)
    if not ok and refcat_file:
        ok = _fit_wcs_scamp(fits_file, refcat_file)

    return str(fits_file), ok


def run_wcs(night_dir, redo=False, nmulti=20):
    """
    Fit WCS solutions for all chip FITS files in a DELVE night.

    Python port of delvered_wcs.pro.

    Reads ``photred.setup`` and the ``logs/WCS.inlist`` file, attempts to
    fit a WCS to each chip using astrometry.net (``solve-field``) or SCAMP,
    and writes ``logs/WCS.success`` / ``logs/WCS.failure`` lists.

    Parameters
    ----------
    night_dir : str or Path
        Night directory containing ``photred.setup`` and ``logs/WCS.inlist``.
    redo : bool
        Re-fit files that already have a WCS solution.
    nmulti : int
        Number of parallel workers.

    Returns
    -------
    tuple of (list, list)
        (success_files, failure_files)
    """
    night_dir = Path(night_dir)
    setup = _read_setup(night_dir / 'photred.setup')
    if not redo:
        redo_str = setup.get('REDO', '0')
        redo = redo_str not in ('0', '', '-1')
    nmulti_cfg = int(setup.get('NMULTI_WCS',
                               setup.get('NMULTI', str(nmulti))))
    nmulti = max(nmulti_cfg, 1)
    search_dist = float(setup.get('SEARCHDIST', '20')) / 60.0  # arcmin→deg

    # Input files
    input_files = _read_list(night_dir / 'logs/WCS.inlist')
    if not input_files:
        print("run_wcs: no files in WCS.inlist")
        return [], []

    abs_files = []
    for f in input_files:
        p = Path(f)
        if not p.is_absolute():
            p = night_dir / p
        abs_files.append(str(p))

    # Build args: look for per-chip refcat
    args_list = []
    for f in abs_files:
        p = Path(f)
        base = p.name
        for ext in ('.fits.fz', '.fits'):
            if base.endswith(ext):
                base = base[:-len(ext)]
                break
        chip_dir = p.parent
        # Per-chip refcat: <base>_refcat.fits.gz or .fits
        refcat = chip_dir / f'{base}_refcat.fits.gz'
        if not refcat.exists():
            refcat = chip_dir / f'{base}_refcat.fits'
        if not refcat.exists():
            refcat = None
        args_list.append((str(p), str(refcat) if refcat else None,
                          search_dist, redo))

    print(f"run_wcs: fitting WCS for {len(args_list)} files")

    if nmulti > 1:
        with ProcessPoolExecutor(max_workers=nmulti) as pool:
            results = list(pool.map(_run_one_wcs, args_list))
    else:
        results = [_run_one_wcs(a) for a in args_list]

    success = [f for f, ok in results if ok]
    failure = [f for f, ok in results if not ok]

    logs_dir = night_dir / 'logs'
    logs_dir.mkdir(exist_ok=True)
    if success:
        _write_list(logs_dir / 'WCS.success', success)
    if failure:
        _write_list(logs_dir / 'WCS.failure', failure)

    print(f"run_wcs: {len(success)} success, {len(failure)} failure")
    return success, failure
