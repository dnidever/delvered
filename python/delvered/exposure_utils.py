"""
Exposure selection and teff-map utilities for the DELVERED pipeline.

Python port of:
  - delvered_pickexposures.pro → pick_exposures()
  - delvered_teffmap.pro       → teff_map()  (also in photometry_utils.py)
"""

import numpy as np
from numpy.lib.recfunctions import append_fields
from astropy.time import Time

from .photometry_utils import get_teff, teff_map, deltasnr_sort
from .utils import create_index, remove_tags


# ---------------------------------------------------------------------------
# Exposure selection
# ---------------------------------------------------------------------------

def pick_exposures(chstr, max_exposures=50):
    """
    Select a quality-optimised subset of chip exposures for brick processing.

    Python port of delvered_pickexposures.pro.

    Computes a calibrated depth metric, effective exposure time (teff),
    and chip-brick overlap fraction.  Exposures are ranked by their
    delta-S/N contribution (preferring new sky coverage over redundant
    coverage).  Up to ``max_exposures`` are selected, with at least
    10 of the overall-deepest included first.

    Parameters
    ----------
    chstr : numpy structured array
        Chip structure.  Required fields: expnum, filter, exptime,
        dao_depth, calib_zpterm, apcor, vx, vy (pixel vertices),
        gain, fwhm, skymode, utdate, uttime.
    max_exposures : int
        Maximum total number of exposures to keep (default 50).

    Returns
    -------
    numpy structured array
        Subset of ``chstr`` for the selected exposures, with the
        temporary working columns removed.
    """
    nchips = len(chstr)

    # --- Calibrated depth ---
    if 'depth' not in chstr.dtype.names:
        depth = (chstr['dao_depth'].astype(float) +
                 2.5 * np.log10(chstr['exptime'].astype(float)) -
                 chstr['calib_zpterm'].astype(float) -
                 chstr['apcor'].astype(float))
        chstr = append_fields(chstr, 'depth', depth, usemask=False)

    # --- Effective exposure time ---
    chstr = get_teff(chstr)

    # --- Chip → brick overlap fraction ---
    tilenx = 3600
    tileny = 3600
    fracoverlap = np.zeros(nchips, dtype=float)
    for i in range(nchips):
        vx = np.asarray(chstr['vx'][i], dtype=float)
        vy = np.asarray(chstr['vy'][i], dtype=float)
        if (vx.max() >= 0 and vx.min() <= tilenx - 1 and
                vy.max() >= 0 and vy.min() <= tileny - 1):
            from skimage.draw import polygon as sk_polygon
            rr, cc = sk_polygon(np.clip(vy, 0, tileny - 1),
                                np.clip(vx, 0, tilenx - 1),
                                shape=(tileny, tilenx))
            fracoverlap[i] = len(rr) / (tilenx * float(tileny))
    if 'fracoverlap' not in chstr.dtype.names:
        chstr = append_fields(chstr, 'fracoverlap', fracoverlap, usemask=False)
    else:
        chstr['fracoverlap'] = fracoverlap

    # --- Build exposure-level summary ---
    eindex = create_index(chstr['expnum'])
    uexp = eindex['value']
    nexposures = len(uexp)

    exp_dt = np.dtype([
        ('expnum', 'U30'), ('filter', 'U4'), ('exptime', float),
        ('obsdate', 'U30'), ('mjd', float), ('nchips', int),
        ('airmass', float), ('zpterm', float), ('nmeasavg', float),
        ('nmeas', int), ('eta', float), ('background', float),
        ('fwhm', float), ('depth', float), ('tau', float),
        ('teff', float), ('fracoverlap', float),
    ])
    expstr = np.zeros(nexposures, dtype=exp_dt)

    for i in range(nexposures):
        ind = eindex['index'][eindex['lo'][i]:eindex['hi'][i] + 1]
        ch = chstr[ind]
        expstr['expnum'][i] = str(uexp[i])
        expstr['filter'][i] = str(ch['filter'][0]).strip()
        expstr['exptime'][i] = float(ch['exptime'][0])
        obsdate = str(ch['utdate'][0]).strip() + 'T' + str(ch['uttime'][0]).strip()
        expstr['obsdate'][i] = obsdate
        try:
            expstr['mjd'][i] = Time(obsdate, format='isot', scale='utc').mjd
        except Exception:
            pass
        expstr['nchips'][i] = len(ind)
        if 'airmass' in ch.dtype.names:
            expstr['airmass'][i] = float(np.mean(ch['airmass']))
        expstr['zpterm'][i] = float(np.mean(ch['calib_zpterm']))
        if 'dao_nsources' in ch.dtype.names:
            expstr['nmeasavg'][i] = float(np.mean(ch['dao_nsources']))
            expstr['nmeas'][i] = int(np.sum(ch['dao_nsources']))
        if 'eta' in ch.dtype.names:
            expstr['eta'][i] = float(np.median(ch['eta']))
        if 'background' in ch.dtype.names:
            expstr['background'][i] = float(np.median(ch['background']))
        if 'fwhm' in ch.dtype.names:
            expstr['fwhm'][i] = float(np.median(ch['fwhm']))
        expstr['depth'][i] = float(np.median(ch['depth']))
        if 'tau' in ch.dtype.names:
            expstr['tau'][i] = float(np.median(ch['tau']))
        if 'teff' in ch.dtype.names:
            expstr['teff'][i] = float(np.median(ch['teff']))
        expstr['fracoverlap'][i] = float(np.sum(ch['fracoverlap']))

    # If already within budget, return all chips (clean up temp columns)
    if nexposures <= max_exposures:
        todel = ['depth', 'eta', 'background', 'tau', 'teff',
                 'fracoverlap', 'mnx', 'mny']
        return remove_tags(chstr, todel)

    # --- Delta-S/N maps per exposure ---
    expstr = append_fields(expstr, 'teffmap',
                           np.zeros((nexposures, 36, 36)),
                           usemask=False)
    for i in range(nexposures):
        ind = eindex['index'][eindex['lo'][i]:eindex['hi'][i] + 1]
        expstr['teffmap'][i] = teff_map(chstr[ind])

    expstr = append_fields(expstr, 'deltasnr', np.zeros(nexposures),
                           usemask=False)
    expstr = append_fields(expstr, 'picked', np.zeros(nexposures, dtype=int),
                           usemask=False)

    # --- Sort overall by delta-S/N (excluding u-band) ---
    not_u = expstr['filter'] != 'u'
    if not_u.any():
        deep_expstr = deltasnr_sort(expstr[not_u])
        top10 = deep_expstr['expnum'][:min(10, len(deep_expstr))]
        for en in top10:
            expstr['picked'][expstr['expnum'] == en] = 1

    # --- Per-filter delta-S/N sort ---
    ufilters = list(set(str(f).strip() for f in expstr['filter']))
    for filt in ufilters:
        fmask = expstr['filter'] == filt
        if fmask.sum() == 0:
            continue
        fexpstr = deltasnr_sort(expstr[fmask])
        expstr[fmask] = fexpstr

    # --- Greedy selection up to max_exposures ---
    count = 0
    while True:
        npicked = int(np.sum(expstr['picked'] == 1))
        if npicked >= max_exposures:
            break
        if np.all(expstr['picked'] != 0):
            break
        for filt in ufilters:
            npicked = int(np.sum(expstr['picked'] == 1))
            if npicked >= max_exposures:
                break
            fmask = (expstr['filter'] == filt) & (expstr['picked'] == 0)
            if not fmask.any():
                continue
            # Pick the one with highest deltasnr
            best = np.argmax(expstr['deltasnr'][fmask])
            pick_idx = np.where(fmask)[0][best]
            expstr['picked'][pick_idx] = 1
        count += 1
        if count > 200:
            break

    picked_exps = expstr['expnum'][expstr['picked'] == 1]
    n_picked = len(picked_exps)
    print(f"pick_exposures: {n_picked} exposures selected from {nexposures}")

    # --- Rebuild per-filter counts ---
    for filt in ufilters:
        n = int(np.sum((expstr['filter'] == filt) & (expstr['picked'] == 1)))
        print(f"  {filt}: {n}")

    # --- Extract chips for selected exposures ---
    chip_mask = np.isin(chstr['expnum'].astype(str), picked_exps.astype(str))
    finalchstr = chstr[chip_mask]

    # Clean up temporary columns
    todel = ['depth', 'eta', 'background', 'tau', 'teff',
             'fracoverlap', 'mnx', 'mny']
    finalchstr = remove_tags(finalchstr, todel)

    return finalchstr
