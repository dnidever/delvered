"""
HEALPix-based coverage map generation for the DELVERED pipeline.

Python port of delvered_coverage.pro.
"""

import os
import numpy as np
from pathlib import Path
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u

try:
    import healpy as hp
    HAS_HEALPY = True
except ImportError:
    HAS_HEALPY = False


# ---------------------------------------------------------------------------
# Main coverage routine
# ---------------------------------------------------------------------------

def delvered_coverage(pix, data_dir, output_dir=None, nside=128, redo=False):
    """
    Generate a HEALPix coverage map for a single nside=128 pixel.

    Python port of delvered_coverage.pro.

    For each sub-pixel (nside=4096), this function calculates the
    sky coverage fraction, median photometric depth, and number of
    exposures for each filter band.

    Parameters
    ----------
    pix : int
        HEALPix ring-scheme pixel number at nside=128.
    data_dir : str or Path
        Root DELVE data directory (contains ``bricks/`` sub-tree).
    output_dir : str or Path, optional
        Directory for output coverage files.  Defaults to
        ``data_dir/coverage/``.
    nside : int
        HEALPix nside for the outer pixel (default 128).
    redo : bool
        Reprocess even if output already exists.

    Returns
    -------
    numpy structured array or None
        Coverage structure, or None on failure.
    """
    if not HAS_HEALPY:
        raise ImportError("healpy is required for coverage maps: pip install healpy")

    data_dir = Path(data_dir)
    if output_dir is None:
        output_dir = data_dir / 'coverage'
    output_dir = Path(output_dir)

    # Output file
    pix_group = pix // 1000
    covfile = output_dir / str(pix_group) / f'{pix}_coverage.fits'
    if covfile.exists() and not redo:
        print(f"delvered_coverage: {covfile} EXISTS and redo=False")
        return fits.getdata(str(covfile), 1)

    # Get coarse pixel geometry
    nside2 = 4096
    vec = hp.pix2vec(nside, pix, nest=False)
    hcenra, hcendec = _vec2radec(vec)

    # Boundary vectors of the coarse pixel
    bound_vecs = hp.boundaries(nside, pix, nest=False)  # (3, 4)
    hra = np.degrees(np.arctan2(bound_vecs[1], bound_vecs[0])) % 360
    hdec = 90.0 - np.degrees(np.arccos(bound_vecs[2]))

    # Get all nside=4096 sub-pixels belonging to this nside=128 pixel
    listpix = _get_subpixels(nside, pix, nside2)
    nlistpix = len(listpix)
    print(f"delvered_coverage: pixel {pix}, {nlistpix} sub-pixels")

    # Sub-pixel centres
    sub_theta, sub_phi = hp.pix2ang(nside2, listpix, nest=False)
    sub_ra = np.degrees(sub_phi)
    sub_dec = 90.0 - np.degrees(sub_theta)

    # Build coverage structure
    filters = ['u', 'g', 'r', 'i', 'z', 'Y', 'VR']
    dt = [
        ('pix', np.int64), ('pix128', np.int64),
        ('ra', np.float64), ('dec', np.float64),
        ('nobj', np.int32), ('coverage', np.float32), ('nexp', np.int32),
    ]
    for b in filters:
        dt += [(f'{b.lower()}coverage', np.float32),
               (f'{b.lower()}nexp', np.int32),
               (f'{b.lower()}depth', np.float32)]
    covstr = np.zeros(nlistpix, dtype=dt)
    covstr['pix'] = listpix
    covstr['pix128'] = pix
    covstr['ra'] = sub_ra
    covstr['dec'] = sub_dec
    for b in filters:
        covstr[f'{b.lower()}depth'] = -9999.0

    # Load brick database to find overlapping bricks
    bricks_file = Path(__file__).parent.parent.parent / 'data' / 'delvemc_bricks_0.25deg.fits.gz'
    if not bricks_file.exists():
        print(f"delvered_coverage: bricks file not found: {bricks_file}")
        return covstr

    bstr = fits.getdata(str(bricks_file), 1)
    # Find bricks that might overlap (rough proximity cut)
    dist = _angular_distance(bstr['ra'], bstr['dec'], hcenra, hcendec)
    nearby_mask = dist < 3.0  # degrees
    nearby_bricks = bstr[nearby_mask]

    if len(nearby_bricks) == 0:
        _save_coverage(covstr, covfile)
        return covstr

    # Load chip metadata from nearby bricks
    all_chips = _load_chip_metadata(nearby_bricks, data_dir, filters)
    if all_chips is None or len(all_chips) == 0:
        _save_coverage(covstr, covfile)
        return covstr

    # For each sub-pixel, compute coverage and depth per filter
    for i in range(nlistpix):
        cenra1 = float(sub_ra[i])
        cendec1 = float(sub_dec[i])
        hpix_bound = hp.boundaries(nside2, listpix[i], nest=False)
        hra1 = np.degrees(np.arctan2(hpix_bound[1], hpix_bound[0])) % 360
        hdec1 = 90.0 - np.degrees(np.arccos(hpix_bound[2]))

        # Build a small raster grid (~20×20 pixels)
        dx = 1e-3  # degrees
        lon0 = float(np.min(hra1)) - cenra1
        lat0 = float(np.min(hdec1)) - cendec1
        # Healpix region mask (approximated as bounding box for speed)
        nx = max(int(np.ptp(hra1) / dx) + 2, 2)
        ny = max(int(np.ptp(hdec1) / dx) + 2, 2)
        maskpix = max(nx * ny, 1)

        for b in filters:
            filt_mask = all_chips['filter'] == b
            filt_chips = all_chips[filt_mask]
            if len(filt_chips) == 0:
                continue

            num_covered = 0
            depth_sum = 0.0
            num_depth = 0

            for chip in filt_chips:
                # Quick bounding box overlap check
                if not _bbox_overlap(chip['vra'], chip['vdec'],
                                     hra1, hdec1):
                    continue
                num_covered += 1
                d = float(chip.get('depth95', 0.0) or chip.get('depth', 0.0))
                if d > 0:
                    depth_sum += d
                    num_depth += 1

            if maskpix > 0:
                covfrac = min(1.0, num_covered / max(maskpix, 1))
            else:
                covfrac = 0.0

            bl = b.lower()
            covstr[f'{bl}coverage'][i] = covfrac
            covstr[f'{bl}nexp'][i] = num_covered
            if num_depth > 0:
                covstr[f'{bl}depth'][i] = depth_sum / num_depth

            # Update totals
            if covfrac > covstr['coverage'][i]:
                covstr['coverage'][i] = covfrac
            covstr['nexp'][i] += num_covered

    _save_coverage(covstr, covfile)
    return covstr


def _get_subpixels(nside, pix, nside2):
    """Return all nside2 pixels whose centre belongs to nside pixel pix."""
    vec = hp.pix2vec(nside, pix, nest=False)
    radius = hp.nside2resol(nside) * 2
    candidates = hp.query_disc(nside2, vec, radius=radius, nest=False)
    theta, phi = hp.pix2ang(nside2, candidates, nest=False)
    pix1 = hp.ang2pix(nside, theta, phi, nest=False)
    return candidates[pix1 == pix]


def _vec2radec(vec):
    """Convert a unit vector to (RA, DEC) in degrees."""
    dec = 90.0 - np.degrees(np.arccos(vec[2]))
    ra = np.degrees(np.arctan2(vec[1], vec[0])) % 360
    return float(ra), float(dec)


def _angular_distance(ra1, dec1, ra2, dec2):
    """Approximate angular distance between (ra1,dec1) and (ra2,dec2) in degrees."""
    c1 = SkyCoord(ra=ra1 * u.deg, dec=dec1 * u.deg, frame='icrs')
    c2 = SkyCoord(ra=ra2 * u.deg, dec=dec2 * u.deg, frame='icrs')
    return c1.separation(c2).deg


def _bbox_overlap(vra1, vdec1, vra2, vdec2):
    """Simple bounding box overlap check."""
    return not (np.max(vra1) < np.min(vra2) or np.min(vra1) > np.max(vra2) or
                np.max(vdec1) < np.min(vdec2) or np.min(vdec1) > np.max(vdec2))


def _load_chip_metadata(bricks, data_dir, filters):
    """Load chip metadata from brick meta FITS files."""
    chips = []
    for brick in bricks:
        bname = str(brick['brickname']).strip() if 'brickname' in brick.dtype.names else str(brick).strip()
        subdir = bname[:4]
        meta_file = data_dir / 'bricks' / subdir / bname / f'{bname}_meta.fits'
        if not meta_file.exists():
            meta_file = str(meta_file) + '.gz'
            if not Path(meta_file).exists():
                continue
        try:
            chstr = fits.getdata(str(meta_file), 2)
            chips.append(chstr)
        except Exception:
            pass
    if not chips:
        return None
    from numpy.lib.recfunctions import stack_arrays
    return stack_arrays(chips, asrecarray=False, usemask=False, autoconvert=True)


def _save_coverage(covstr, covfile):
    """Write coverage structure to a FITS file."""
    covfile = Path(covfile)
    covfile.parent.mkdir(parents=True, exist_ok=True)
    hdu = fits.BinTableHDU(covstr)
    hdu.writeto(str(covfile), overwrite=True)
    print(f"delvered_coverage: wrote {covfile}")
