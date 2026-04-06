"""
Final DELVE-MC catalogue construction.

Python port of delvered_finalcat.pro → build_final_cat()
"""

import time
from pathlib import Path

import numpy as np
from astropy.io import fits

try:
    import healpy as hp
    HAS_HEALPY = True
except ImportError:
    HAS_HEALPY = False


# Final catalogue schema
_FINAL_DTYPE = np.dtype([
    ('brick', 'U20'),
    ('id', 'U80'),
    ('ring256', np.int32),
    ('ra', np.float64),
    ('dec', np.float64),
    ('ndet', np.int32),
    ('gmag', np.float32),
    ('gerr', np.float32),
    ('gscatter', np.float32),
    ('ndetg', np.int32),
    ('rmag', np.float32),
    ('rerr', np.float32),
    ('rscatter', np.float32),
    ('ndetr', np.int32),
    ('imag', np.float32),
    ('ierr', np.float32),
    ('iscatter', np.float32),
    ('ndeti', np.int32),
    ('zmag', np.float32),
    ('zerr', np.float32),
    ('zscatter', np.float32),
    ('ndetz', np.int32),
    ('chi', np.float32),
    ('sharp', np.float32),
    ('prob', np.float32),
    ('sfd_ebv', np.float32),
    ('mag_auto', np.float32),
    ('magerr_auto', np.float32),
    ('asemi', np.float32),
    ('bsemi', np.float32),
    ('theta', np.float32),
    ('ellipticity', np.float32),
    ('fwhm', np.float32),
])

_FINAL_DEFAULTS = {
    'gmag': 99.99, 'gerr': 9.99, 'gscatter': 99.99,
    'rmag': 99.99, 'rerr': 9.99, 'rscatter': 99.99,
    'imag': 99.99, 'ierr': 9.99, 'iscatter': 99.99,
    'zmag': 99.99, 'zerr': 9.99, 'zscatter': 99.99,
    'mag_auto': 999.0, 'magerr_auto': 999.0,
    'asemi': 999.0, 'bsemi': 999.0, 'theta': 999.0,
    'ellipticity': 999.0, 'fwhm': 999.0,
}


def _cat2final(cat, brick, brickstr_row):
    """
    Convert a brick catalogue to the final schema.

    Parameters
    ----------
    cat : numpy structured array
        Brick object catalogue (already filtered to brickuniq=1 objects).
    brick : str
        Brick name.
    brickstr_row : row
        Corresponding row from the brick list (unused but kept for API parity).

    Returns
    -------
    numpy structured array
        Final schema array.
    """
    ncat = len(cat)
    new = np.zeros(ncat, dtype=_FINAL_DTYPE)

    # Apply defaults
    for col, val in _FINAL_DEFAULTS.items():
        new[col] = val

    # Copy matching fields
    cat_cols = {name.lower(): name for name in cat.dtype.names}
    for outname in _FINAL_DTYPE.names:
        if outname in ('brick', 'id', 'ring256', 'sfd_ebv'):
            continue
        if outname in cat_cols:
            try:
                new[outname] = cat[cat_cols[outname]]
            except Exception:
                pass

    new['brick'] = brick
    new['id'] = np.array([f'{brick}.{str(cid).strip()}' for cid in cat['id']])

    # HEALPix ring-256
    if HAS_HEALPY:
        theta = np.radians(90.0 - cat['dec'].astype(float))
        phi = np.radians(cat['ra'].astype(float))
        new['ring256'] = hp.ang2pix(256, theta, phi, nest=False)

    # SFD E(B-V)
    if 'ebv' in cat_cols:
        new['sfd_ebv'] = cat[cat_cols['ebv']]

    return new


def build_final_cat(delve_dir, delvereddir=None, db_file=None,
                    table='object'):
    """
    Build the final DELVE-MC catalogue from per-brick FITS files.

    Python port of delvered_finalcat.pro.

    Iterates over all brick directories, loads the per-brick object
    catalogue, keeps only objects inside the brick's unique area
    (``brickuniq == 1``), converts to the final schema, and writes
    each brick's data into the shared SQLite database via
    :func:`delvered_db.writecat2db`.

    Parameters
    ----------
    delve_dir : str or Path
        Root DELVE data directory (contains ``bricks/``).
    delvereddir : str or Path, optional
        Root of the delvered project tree (for the brick list).
        Defaults to the package data directory.
    db_file : str or Path, optional
        SQLite database file.  Defaults to
        ``{delve_dir}/bricks/db/delvered_final.db``.
    table : str
        Database table name (default ``'object'``).

    Returns
    -------
    int
        Total number of objects written.
    """
    from .delvered_db import writecat2db, createindexdb

    t0 = time.time()

    if not HAS_HEALPY:
        raise ImportError("healpy is required: pip install healpy")

    delve_dir = Path(delve_dir)
    brick_dir = delve_dir / 'bricks'

    if db_file is None:
        db_file = brick_dir / 'db' / 'delvered_final.db'
    db_file = Path(db_file)
    db_file.parent.mkdir(parents=True, exist_ok=True)

    # Load brick list
    if delvereddir is not None:
        bricks_file = Path(delvereddir) / 'data' / 'delvemc_bricks_0.25deg.fits.gz'
    else:
        bricks_file = (Path(__file__).parent.parent.parent /
                       'data' / 'delvemc_bricks_0.25deg.fits.gz')

    brickstr = fits.getdata(str(bricks_file), 1)
    brickstr_names = np.char.strip(brickstr['brickname'].astype(str))

    # Find brick directories
    brick_dirs = sorted(
        d for d in brick_dir.glob('????/????????')
        if d.is_dir()
    )
    nbdirs = len(brick_dirs)
    print(f"build_final_cat: {nbdirs} brick directories")
    if nbdirs == 0:
        return 0

    bricks_found = np.array([d.name for d in brick_dirs])
    in_list = np.isin(brickstr_names, bricks_found)
    matched_names = brickstr_names[in_list]
    matched_bstr = brickstr[in_list]

    total = 0
    for i, bname in enumerate(matched_names):
        idir = brick_dir / bname[:4] / bname
        brow = matched_bstr[i]
        print(f"{i+1}  {bname}  ra={float(brow['ra']):.4f}  dec={float(brow['dec']):.4f}")

        catfile = idir / f'{bname}.fits.gz'
        metafile = idir / f'{bname}_meta.fits'

        missing = []
        if not catfile.exists():
            missing.append(str(catfile))
        if not metafile.exists():
            missing.append(str(metafile))
        if missing:
            for m in missing:
                print(f"  {m} NOT FOUND")
            continue

        try:
            cat = fits.getdata(str(catfile), 1)
        except Exception as exc:
            print(f"  ERROR loading {catfile}: {exc}")
            continue

        ncat = len(cat)
        print(f"  {ncat} objects")

        if 'brickuniq' not in cat.dtype.names:
            print("  WARNING: brickuniq column not found; using all objects")
            inside = np.arange(ncat)
        else:
            inside = np.where(cat['brickuniq'] == 1)[0]

        ninside = len(inside)
        print(f"  Keeping {ninside} objects inside the unique brick area")
        if ninside == 0:
            print("  No objects left to save")
            continue

        cat = cat[inside]
        new = _cat2final(cat, bname, brow)

        print("  Writing to database")
        writecat2db(new, str(db_file), table=table)
        total += ninside

    # Index the database
    print("build_final_cat: indexing the table")
    createindexdb(str(db_file), col='ra', table=table, unique=False)
    createindexdb(str(db_file), col='dec', table=table, unique=False)

    elapsed = time.time() - t0
    print(f"build_final_cat: done, {total} objects in {elapsed:.1f}s")
    return total
