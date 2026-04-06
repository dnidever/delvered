"""
MW stellar model file conversion for the DELVERED pipeline.

Python port of convert_mwmodels.pro → convert_mw_models()
"""

import gzip
import subprocess
from pathlib import Path

import numpy as np
from astropy.io import fits

try:
    import healpy as hp
    HAS_HEALPY = True
except ImportError:
    HAS_HEALPY = False


def convert_mw_models(files, delvereddir=None):
    """
    Convert ASCII MW stellar model files to FITS and add HEALPix ring256 index.

    Python port of convert_mwmodels.pro.

    For each ASCII model file, reads the tabular data, maps it to the
    DELVERED brick schema (renaming ``_fe_h_`` → ``feh``), looks up the
    brick centre coordinates from the brick list, computes the nside=256
    ring-scheme HEALPix pixel index, and writes a gzip-compressed FITS
    file alongside the original.

    Parameters
    ----------
    files : list of str or Path
        ASCII model files to convert.  The brick name is inferred from
        the second underscore-delimited component of the filename
        (e.g. ``model_1234p567_...txt`` → brick ``'1234p567'``).
    delvereddir : str or Path, optional
        Root of the delvered project tree, used to locate
        ``data/delvemc_bricks_0.25deg.fits.gz``.
        Defaults to the package data directory.

    Returns
    -------
    list of str
        Paths to the written ``.fits.gz`` output files.
    """
    if not HAS_HEALPY:
        raise ImportError("healpy is required: pip install healpy")

    # Locate brick list
    if delvereddir is not None:
        bricks_file = Path(delvereddir) / 'data' / 'delvemc_bricks_0.25deg.fits.gz'
    else:
        bricks_file = (Path(__file__).parent.parent.parent /
                       'data' / 'delvemc_bricks_0.25deg.fits.gz')

    brickstr = fits.getdata(str(bricks_file), 1)
    brick_names = np.char.strip(brickstr['brickname'].astype(str))

    output_files = []
    print(f"convert_mw_models: {len(files)} file(s) to process")

    for i, fpath in enumerate(files):
        fpath = Path(fpath)
        print(f"{i+1}  {fpath}")

        if not fpath.exists() or fpath.stat().st_size == 0:
            print(f"  {fpath} is empty or missing")
            continue

        # Infer brick name from filename: model_<brick>_...
        parts = fpath.name.split('_')
        brick = parts[1] if len(parts) > 1 else parts[0]

        # Load ASCII table (first row is a header)
        try:
            arr = np.genfromtxt(str(fpath), names=True, dtype=None,
                                encoding='ascii')
        except Exception as exc:
            print(f"  ERROR reading {fpath}: {exc}")
            continue
        if arr.ndim == 0:
            arr = arr.reshape(1)
        narr = len(arr)
        print(f"  {narr} rows")

        # Find brick in list
        match = np.where(brick_names == brick)[0]
        if len(match) == 0:
            print(f"  No match for brick {brick!r}")
            continue
        brick_ra = float(brickstr['ra'][match[0]])
        brick_dec = float(brickstr['dec'][match[0]])

        # HEALPix ring-256 pixel
        theta = np.radians(90.0 - brick_dec)
        phi = np.radians(brick_ra)
        hpix = int(hp.ang2pix(256, theta, phi, nest=False))

        # Build output structured array
        dt = np.dtype([
            ('brick', 'U20'), ('ra', np.float64), ('dec', np.float64),
            ('ring256', np.int32),
            ('u', np.float32), ('g', np.float32), ('r', np.float32),
            ('i', np.float32), ('z', np.float32),
            ('feh', np.float32), ('age', np.float32), ('popid', np.float32),
        ])
        new = np.zeros(narr, dtype=dt)
        new['brick'] = brick
        new['ra'] = brick_ra
        new['dec'] = brick_dec
        new['ring256'] = hpix

        # Copy fields by name (case-insensitive); rename _fe_h_ → feh
        col_map = {name.lower(): name for name in arr.dtype.names}
        for outfield in ('u', 'g', 'r', 'i', 'z', 'age', 'popid'):
            if outfield in col_map:
                new[outfield] = arr[col_map[outfield]].astype(np.float32)
        for key in col_map:
            if 'fe' in key and 'h' in key:
                new['feh'] = arr[col_map[key]].astype(np.float32)
                break

        # Write FITS
        outfile = Path(str(fpath) + '.fits')
        print(f"  Writing to {outfile}")
        hdu = fits.BinTableHDU(new)
        hdu.writeto(str(outfile), overwrite=True)
        subprocess.run(['gzip', '-f', str(outfile)], check=False)
        outfile_gz = Path(str(outfile) + '.gz')
        output_files.append(str(outfile_gz))

    return output_files
