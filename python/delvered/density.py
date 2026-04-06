"""
LMC stellar density analysis for DELVERED bricks.

Python port of:
  - brick_lmc_density.pro  → brick_lmc_density()
  - brick_lmc_density1.pro → brick_lmc_density1()
  - brick_lmc_density2.pro → brick_lmc_density2()
  - brick_lmc_density3.pro → brick_lmc_density3()
"""

import os
import gzip
import shutil
import numpy as np
from pathlib import Path
from astropy.io import fits


# Extinction coefficients (magnitudes per unit E(B-V))
_EXT_G = 3.303
_EXT_R = 2.285
_EXT_I = 1.698


# ---------------------------------------------------------------------------
# Basic gnomonic area helper
# ---------------------------------------------------------------------------

def _sky_area_deg2(ra, dec, cenra, cendec, dx=0.005):
    """Estimate occupied sky area (deg²) using a 2D histogram."""
    from .extinction import _gnomonic_reverse  # noqa: F401
    # Forward gnomonic projection
    ra_r = np.deg2rad(ra)
    dec_r = np.deg2rad(dec)
    ra0 = np.deg2rad(cenra)
    dec0 = np.deg2rad(cendec)
    denom = (np.sin(dec0) * np.sin(dec_r) +
             np.cos(dec0) * np.cos(dec_r) * np.cos(ra_r - ra0))
    xi = np.cos(dec_r) * np.sin(ra_r - ra0) / denom
    eta = (np.cos(dec0) * np.sin(dec_r) -
           np.sin(dec0) * np.cos(dec_r) * np.cos(ra_r - ra0)) / denom
    xi = np.rad2deg(xi)
    eta = np.rad2deg(eta)

    x_edges = np.arange(xi.min() - dx, xi.max() + 2 * dx, dx)
    y_edges = np.arange(eta.min() - dx, eta.max() + 2 * dx, dx)
    H, _, _ = np.histogram2d(xi, eta, bins=[x_edges, y_edges])
    n_pix = np.sum(H > 0)
    return n_pix * dx**2


# ---------------------------------------------------------------------------
# Basic LMC density (brick_lmc_density.pro)
# ---------------------------------------------------------------------------

def brick_lmc_density(brick, bricks_dir, output_dir=None, redo=False):
    """
    Calculate LMC stellar density statistics for a single brick.

    Python port of brick_lmc_density.pro.

    Parameters
    ----------
    brick : str
        Brick name (e.g. ``'1234p567'``).
    bricks_dir : str or Path
        Root directory containing per-brick sub-directories.
    output_dir : str or Path, optional
        Where to write output FITS files.  Defaults to
        ``bricks_dir/summary/lmc/``.
    redo : bool
        Reprocess even if output exists.

    Returns
    -------
    numpy structured array or None
    """
    bricks_dir = Path(bricks_dir)
    if output_dir is None:
        output_dir = bricks_dir / 'summary' / 'lmc'
    output_dir = Path(output_dir)

    outfile = output_dir / f'{brick}.fits'
    if outfile.exists() and not redo:
        print(f"brick_lmc_density: {outfile} EXISTS, skipping")
        return fits.getdata(str(outfile), 1)

    subdir = brick[:4]
    bdir = bricks_dir / subdir / brick
    objfile = bdir / f'{brick}_object.fits.gz'
    if not objfile.exists():
        objfile = bdir / f'{brick}_object.fits'
    if not objfile.exists():
        print(f"brick_lmc_density: {objfile} NOT FOUND")
        return None

    obj = fits.getdata(str(objfile), 1)
    nobj = len(obj)
    print(f"brick_lmc_density: {brick}  {nobj} objects")

    ra = obj['ra'].astype(float)
    dec = obj['dec'].astype(float)

    # Field centre
    if np.ptp(ra) > 100:
        ra_adj = np.where(ra > 180, ra - 360, ra)
        cenra = float(np.mean([ra_adj.min(), ra_adj.max()]))
        if cenra < 0:
            cenra += 360
    else:
        cenra = float(np.mean([ra.min(), ra.max()]))
    cendec = float(np.mean([dec.min(), dec.max()]))

    # Build output structure
    brkstr = np.zeros(1, dtype=[
        ('brick', 'U20'), ('ra', float), ('dec', float), ('area', float),
        ('gdepth95', float), ('rdepth95', float), ('idepth95', float),
        ('nobj', int),
        ('nlmc1', int), ('nlmc2', int),
        ('nlmc220', int), ('nlmc225', int), ('nlmc230', int), ('nlmc235', int),
        ('nback1', int), ('nback2', int),
    ])
    brkstr['brick'] = brick
    brkstr['ra'] = cenra
    brkstr['dec'] = cendec
    brkstr['nobj'] = nobj
    for col in ('gdepth95', 'rdepth95', 'idepth95'):
        brkstr[col] = 99.99
    for col in ('nlmc1', 'nlmc2', 'nlmc220', 'nlmc225', 'nlmc230', 'nlmc235',
                'nback1', 'nback2'):
        brkstr[col] = -1

    # Quality cuts
    gmag = obj['gmag'].astype(float) if 'gmag' in obj.dtype.names else None
    imag = obj['imag'].astype(float) if 'imag' in obj.dtype.names else None
    sharp = obj['sharp'].astype(float) if 'sharp' in obj.dtype.names else None
    chi = obj['chi'].astype(float) if 'chi' in obj.dtype.names else None
    prob = obj['prob'].astype(float) if 'prob' in obj.dtype.names else None

    quality = np.ones(nobj, dtype=bool)
    if gmag is not None:
        quality &= gmag < 50
    if imag is not None:
        quality &= imag < 50
    if sharp is not None:
        quality &= np.abs(sharp) < 1
    if chi is not None:
        quality &= chi < 3
    if prob is not None:
        quality &= prob > 0.5

    if quality.any():
        gm = gmag[quality] if gmag is not None else np.full(quality.sum(), 99.99)
        im = imag[quality] if imag is not None else np.full(quality.sum(), 99.99)
        gi = gm - im

        brkstr['nlmc1'] = int(np.sum((gi >= -0.10) & (gi <= 0.56) &
                                     (gm > 21.8) & (gm < 22.8)))
        brkstr['nlmc2'] = int(np.sum((gi >= -0.10) & (gi <= 0.56) &
                                     (gm > 22.0) & (gm < 24.0)))
        brkstr['nlmc220'] = int(np.sum((gi >= -0.10) & (gi <= 0.50) &
                                       (gm > 22.0) & (gm < 22.5)))
        brkstr['nlmc225'] = int(np.sum((gi >= 0.0) & (gi <= 0.50) &
                                       (gm > 22.5) & (gm < 23.0)))
        brkstr['nlmc230'] = int(np.sum((gi >= 0.10) & (gi <= 0.54) &
                                       (gm > 23.0) & (gm < 23.5)))
        brkstr['nlmc235'] = int(np.sum((gi >= 0.10) & (gi <= 0.65) &
                                       (gm > 23.05) & (gm < 24.0)))
        # Background regions
        brkstr['nback1'] = int(np.sum((gi >= 0.60) & (gi <= 0.80) &
                                      (gm >= 22.5) & (gm <= 23.0)))
        brkstr['nback2'] = int(np.sum((gi >= 0.40) & (gi <= 0.60) &
                                      (gm >= 20.0) & (gm <= 20.5)))

    # Per-band 95th-percentile depth
    for band, col in [('g', 'gdepth95'), ('r', 'rdepth95'), ('i', 'idepth95')]:
        magcol = f'{band}mag'
        if magcol in obj.dtype.names:
            mags = obj[magcol].astype(float)
            good = mags < 50
            if good.sum() > 5:
                brkstr[col] = float(np.percentile(mags[good], 95))

    # Unique sky area
    brkstr['area'] = float(_sky_area_deg2(ra, dec, cenra, cendec))

    # Save
    output_dir.mkdir(parents=True, exist_ok=True)
    hdu = fits.BinTableHDU(brkstr)
    hdu.writeto(str(outfile), overwrite=True)
    print(f"brick_lmc_density: wrote {outfile}")
    return brkstr


# ---------------------------------------------------------------------------
# Enhanced density with extinction correction (brick_lmc_density1.pro)
# ---------------------------------------------------------------------------

def brick_lmc_density1(brick, bricks_dir, output_dir=None, redo=False):
    """
    Enhanced LMC density with extinction corrections and blue star counts.

    Python port of brick_lmc_density1.pro.

    Parameters
    ----------
    brick : str
        Brick name.
    bricks_dir : str or Path
        Root bricks directory.
    output_dir : str or Path, optional
        Output directory.
    redo : bool
        Reprocess if output exists.

    Returns
    -------
    numpy structured array or None
    """
    bricks_dir = Path(bricks_dir)
    if output_dir is None:
        output_dir = bricks_dir / 'summary' / 'lmc1'
    output_dir = Path(output_dir)
    outfile = output_dir / f'{brick}.fits'
    if outfile.exists() and not redo:
        return fits.getdata(str(outfile), 1)

    subdir = brick[:4]
    bdir = bricks_dir / subdir / brick
    objfile = bdir / f'{brick}_object.fits.gz'
    if not objfile.exists():
        objfile = bdir / f'{brick}_object.fits'
    if not objfile.exists():
        print(f"brick_lmc_density1: {objfile} NOT FOUND")
        return None

    obj = fits.getdata(str(objfile), 1)
    nobj = len(obj)
    print(f"brick_lmc_density1: {brick}  {nobj} objects")

    ra = obj['ra'].astype(float)
    dec = obj['dec'].astype(float)
    cenra = float(np.mean([ra.min(), ra.max()]))
    cendec = float(np.mean([dec.min(), dec.max()]))

    gmag = obj['gmag'].astype(float) if 'gmag' in obj.dtype.names else np.full(nobj, 99.99)
    imag = obj['imag'].astype(float) if 'imag' in obj.dtype.names else np.full(nobj, 99.99)
    ebv = obj['ebv'].astype(float) if 'ebv' in obj.dtype.names else np.zeros(nobj)

    # Extinction-corrected magnitudes
    gmag_c = gmag - ebv * _EXT_G
    imag_c = imag - ebv * _EXT_I

    brkstr = np.zeros(1, dtype=[
        ('brick', 'U20'), ('ra', float), ('dec', float), ('area', float),
        ('gdepth95', float), ('idepth95', float),
        ('gdepth5', float), ('idepth5', float),
        ('nobj', int), ('nlmc', int), ('nblue', int),
    ])
    brkstr['brick'] = brick
    brkstr['ra'] = cenra
    brkstr['dec'] = cendec
    brkstr['nobj'] = nobj
    for col in ('gdepth95', 'idepth95', 'gdepth5', 'idepth5'):
        brkstr[col] = 99.99

    gi_c = gmag_c - imag_c
    quality = (gmag < 50) & (imag < 50)
    if quality.any():
        lmc = quality & (gi_c >= -0.10) & (gi_c <= 0.56) & (gmag_c > 21.8) & (gmag_c < 22.8)
        brkstr['nlmc'] = int(lmc.sum())
        blue = quality & (gi_c < 0.0) & (gmag_c > 18) & (gmag_c < 21.8)
        brkstr['nblue'] = int(blue.sum())

    for band, magcol, dcol5, dcol95 in [
        ('g', gmag, 'gdepth5', 'gdepth95'),
        ('i', imag, 'idepth5', 'idepth95'),
    ]:
        good = magcol < 50
        if good.sum() > 5:
            brkstr[dcol95] = float(np.percentile(magcol[good], 95))
            brkstr[dcol5] = float(np.percentile(magcol[good], 5))

    brkstr['area'] = float(_sky_area_deg2(ra, dec, cenra, cendec))

    output_dir.mkdir(parents=True, exist_ok=True)
    fits.writeto(str(outfile), brkstr, overwrite=True)
    return brkstr


# ---------------------------------------------------------------------------
# Density using joint catalogue (brick_lmc_density2.pro)
# ---------------------------------------------------------------------------

def brick_lmc_density2(brick, bricks_dir, output_dir=None, redo=False):
    """
    LMC density using r-band selection from joint-processed catalogue.

    Python port of brick_lmc_density2.pro.

    Parameters
    ----------
    brick : str
    bricks_dir : str or Path
    output_dir : str or Path, optional
    redo : bool

    Returns
    -------
    numpy structured array or None
    """
    bricks_dir = Path(bricks_dir)
    if output_dir is None:
        output_dir = bricks_dir / 'summary' / 'lmc2'
    output_dir = Path(output_dir)
    outfile = output_dir / f'{brick}.fits'
    if outfile.exists() and not redo:
        return fits.getdata(str(outfile), 1)

    subdir = brick[:4]
    bdir = bricks_dir / subdir / brick
    # Prefer joint object catalogue
    for fname in (f'{brick}_joint_object.fits.gz', f'{brick}_object.fits.gz',
                  f'{brick}_object.fits'):
        objfile = bdir / fname
        if objfile.exists():
            break
    else:
        print(f"brick_lmc_density2: no object file for {brick}")
        return None

    obj = fits.getdata(str(objfile), 1)
    nobj = len(obj)

    gmag = obj['gmag'].astype(float) if 'gmag' in obj.dtype.names else np.full(nobj, 99.99)
    rmag = obj['rmag'].astype(float) if 'rmag' in obj.dtype.names else np.full(nobj, 99.99)
    ebv = obj['ebv'].astype(float) if 'ebv' in obj.dtype.names else np.zeros(nobj)

    gmag_c = gmag - ebv * _EXT_G
    rmag_c = rmag - ebv * _EXT_R
    gr_c = gmag_c - rmag_c

    brkstr = np.zeros(1, dtype=[
        ('brick', 'U20'), ('ra', float), ('dec', float),
        ('nobj', int), ('nlmc', int),
        ('gdepth95', float), ('rdepth95', float),
        ('gdepth10sig', float), ('rdepth10sig', float),
    ])
    brkstr['brick'] = brick
    brkstr['nobj'] = nobj
    for col in ('gdepth95', 'rdepth95', 'gdepth10sig', 'rdepth10sig'):
        brkstr[col] = 99.99

    lmc = ((gr_c >= 0.10) & (gr_c <= 0.40) & (rmag_c > 21.8) & (rmag_c < 22.7))
    brkstr['nlmc'] = int(lmc.sum())

    # 95th percentile depth and 10-sigma depth
    for band, magcol, col95, col10 in [
        ('g', gmag, 'gdepth95', 'gdepth10sig'),
        ('r', rmag, 'rdepth95', 'rdepth10sig'),
    ]:
        errcol = f'{band}err'
        good = magcol < 50
        if good.sum() > 5:
            brkstr[col95] = float(np.percentile(magcol[good], 95))
        if errcol in obj.dtype.names:
            err = obj[errcol].astype(float)
            snr = 1.087 / np.where(err > 0, err, np.inf)
            sn10 = (snr >= 9.5) & (snr <= 10.5)
            if sn10.sum() > 3:
                brkstr[col10] = float(np.median(magcol[sn10]))

    output_dir.mkdir(parents=True, exist_ok=True)
    fits.writeto(str(outfile), brkstr, overwrite=True)
    return brkstr


# ---------------------------------------------------------------------------
# Comprehensive multi-band density (brick_lmc_density3.pro)
# ---------------------------------------------------------------------------

def brick_lmc_density3(brick, bricks_dir, output_dir=None, redo=False):
    """
    Comprehensive multi-band density analysis for both standard and joint catalogues.

    Python port of brick_lmc_density3.pro.

    Calculates per-band statistics (95th-percentile depth, 10-sigma depth,
    bright-star photometric RMS) for both the standard and joint object
    catalogues.

    Parameters
    ----------
    brick : str
    bricks_dir : str or Path
    output_dir : str or Path, optional
    redo : bool

    Returns
    -------
    numpy structured array or None
    """
    bands = ['u', 'g', 'r', 'i', 'z', 'y']
    bricks_dir = Path(bricks_dir)
    if output_dir is None:
        output_dir = bricks_dir / 'summary' / 'lmc3'
    output_dir = Path(output_dir)
    outfile = output_dir / f'{brick}.fits'
    if outfile.exists() and not redo:
        return fits.getdata(str(outfile), 1)

    subdir = brick[:4]
    bdir = bricks_dir / subdir / brick

    # Build dtype
    dt = [('brick', 'U20')]
    prefixes = ['', 'j']  # '' for standard, 'j' for joint
    for p in prefixes:
        for b in bands:
            dt += [
                (f'{p}{b}depth95', float),
                (f'{p}{b}depth10sig', float),
                (f'{p}{b}rms', float),
            ]
        dt += [(f'{p}nlmc', int)]
    brkstr = np.zeros(1, dtype=dt)
    brkstr['brick'] = brick
    for col in brkstr.dtype.names:
        if col != 'brick':
            brkstr[col] = 99.99 if 'depth' in col or 'rms' in col else -1

    for prefix, fname_pattern in [
        ('', f'{brick}_object.fits.gz'),
        ('j', f'{brick}_joint_object.fits.gz'),
    ]:
        objfile = bdir / fname_pattern
        if not objfile.exists():
            objfile = bdir / fname_pattern.replace('.gz', '')
        if not objfile.exists():
            continue

        obj = fits.getdata(str(objfile), 1)
        nobj = len(obj)
        gmag = obj['gmag'].astype(float) if 'gmag' in obj.dtype.names else np.full(nobj, 99.99)
        rmag = obj['rmag'].astype(float) if 'rmag' in obj.dtype.names else np.full(nobj, 99.99)
        ebv = obj['ebv'].astype(float) if 'ebv' in obj.dtype.names else np.zeros(nobj)

        # LMC count
        gr = gmag - rmag
        lmc = (gr >= 0.10) & (gr <= 0.40) & (rmag > 21.8) & (rmag < 22.7)
        brkstr[f'{prefix}nlmc'] = int(lmc.sum())

        for b in bands:
            magcol = f'{b}mag'
            errcol = f'{b}err'
            if magcol not in obj.dtype.names:
                continue
            mags = obj[magcol].astype(float)
            good = mags < 50

            # 95th percentile depth
            if good.sum() > 5:
                brkstr[f'{prefix}{b}depth95'] = float(np.percentile(mags[good], 95))

            # 10-sigma depth
            if errcol in obj.dtype.names:
                errs = obj[errcol].astype(float)
                snr = 1.087 / np.where(errs > 0, errs, np.inf)
                sn10 = good & (snr >= 9.5) & (snr <= 10.5)
                if sn10.sum() > 3:
                    brkstr[f'{prefix}{b}depth10sig'] = float(np.median(mags[sn10]))

                # Bright-star RMS (14 < mag < 21)
                bright = good & (mags > 14) & (mags < 21)
                if bright.sum() > 5:
                    # Bin by magnitude and measure scatter
                    mbins = np.arange(14, 21.5, 0.5)
                    rms_vals = []
                    for lo, hi in zip(mbins[:-1], mbins[1:]):
                        mask = bright & (mags > lo) & (mags <= hi)
                        if mask.sum() > 3:
                            rms_vals.append(np.std(errs[mask]))
                    if rms_vals:
                        brkstr[f'{prefix}{b}rms'] = float(np.median(rms_vals))

    output_dir.mkdir(parents=True, exist_ok=True)
    fits.writeto(str(outfile), brkstr, overwrite=True)
    print(f"brick_lmc_density3: wrote {outfile}")
    return brkstr
