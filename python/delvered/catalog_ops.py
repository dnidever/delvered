"""
Catalog combination and averaging operations for the DELVERED pipeline.

Python port of:
  - add_meas2obj.pro            → add_meas_to_obj()
  - delvered_simpleavgmeas.pro  → simple_avg_meas()
  - delvered_avgmeas.pro        → avg_meas()
  - delvered_avgmeas_loop.pro   → avg_meas_loop()
"""

import numpy as np
from numpy.lib.recfunctions import append_fields
from .utils import create_index


# Default band ordering
_BANDS = ['u', 'g', 'r', 'i', 'z', 'y']


# ---------------------------------------------------------------------------
# Object schema helpers
# ---------------------------------------------------------------------------

def _make_obj_schema_simple():
    """Dtype for the simplified object catalogue."""
    dt = [
        ('objid', 'U50'),
        ('depthflag', np.int16),
        ('ra', np.float64), ('dec', np.float64),
        ('rarms', np.float32), ('decrms', np.float32),
        ('ndet', np.int32),
    ]
    for b in _BANDS:
        bu = b.upper()
        dt += [
            (f'{b}mag', np.float32), (f'{b}err', np.float32),
            (f'{b}rms', np.float32), (f'ndet{b}', np.int32),
        ]
    dt += [('chi', np.float32), ('sharp', np.float32)]
    return np.dtype(dt)


def _make_obj_schema_full():
    """Dtype for the full object catalogue (from delvered_avgmeas)."""
    dt = [
        ('objid', 'U50'), ('brick', 'U20'),
        ('depthflag', np.int16), ('nalfdetiter', np.int16),
        ('neimerged', np.int16),
        ('ra', np.float64), ('dec', np.float64),
        ('rarms', np.float32), ('decrms', np.float32),
        ('ndet', np.int32), ('ndetall', np.int32),
        ('mlon', np.float64), ('mlat', np.float64),
    ]
    for b in _BANDS:
        bu = b.upper()
        dt += [
            (f'{b}mag', np.float32), (f'{b}err', np.float32),
            (f'{b}rms', np.float32), (f'ndet{b}', np.int32),
            (f'{b}magall', np.float32), (f'{b}errall', np.float32),
            (f'{b}rmsall', np.float32), (f'ndetall{b}', np.int32),
        ]
    dt += [
        ('chi', np.float32), ('sharp', np.float32),
        ('prob', np.float32), ('ebv', np.float32),
        ('mag_auto', np.float32), ('magerr_auto', np.float32),
        ('asemi', np.float32), ('bsemi', np.float32),
        ('theta', np.float32), ('ellipticity', np.float32),
        ('fwhm', np.float32),
        ('rmsvar', np.float32), ('madvar', np.float32),
        ('iqrvar', np.float32), ('etavar', np.float32),
        ('jvar', np.float32), ('kvar', np.float32),
        ('chivar', np.float32), ('romsvar', np.float32),
        ('variable10sig', np.int8), ('nsigvar', np.float32),
        ('gaia_match', np.int8), ('gaia_xdist', np.float32),
        ('gaia_sourceid', np.int64), ('gaia_ra', np.float64),
        ('gaia_ra_error', np.float32), ('gaia_dec', np.float64),
        ('gaia_dec_error', np.float32), ('gaia_parallax', np.float32),
        ('gaia_parallax_error', np.float32), ('gaia_pmra', np.float32),
        ('gaia_pmra_error', np.float32), ('gaia_pmdec', np.float32),
        ('gaia_pmdec_error', np.float32), ('gaia_gmag', np.float32),
        ('gaia_gmag_error', np.float32), ('gaia_bpmag', np.float32),
        ('gaia_bpmag_error', np.float32), ('gaia_rpmag', np.float32),
        ('gaia_rpmag_error', np.float32),
        ('brickuniq', np.uint8),
    ]
    return np.dtype(dt)


# ---------------------------------------------------------------------------
# Flux-weighted averaging utilities
# ---------------------------------------------------------------------------

def _flux_weighted_avg(mags, errs):
    """
    Compute flux-weighted average magnitude and propagated error.

    Parameters
    ----------
    mags, errs : 1-D arrays
        Magnitudes and errors.

    Returns
    -------
    avg_mag, avg_err : float
        Averaged magnitude and error.  Returns (99.99, 9.99) if result
        is non-finite.
    """
    wt = 1.0 / errs**2
    flux = 2.5118864**mags
    total_fluxwt = np.sum(flux * wt)
    total_wt = np.sum(wt)
    if total_wt == 0:
        return 99.99, 9.99
    new_flux = total_fluxwt / total_wt
    new_mag = 2.5 * np.log10(new_flux)
    new_err = np.sqrt(1.0 / total_wt)
    if not np.isfinite(new_mag):
        return 99.99, 9.99
    return float(new_mag), float(new_err)


# ---------------------------------------------------------------------------
# Add a single set of measurements to an existing object catalogue
# ---------------------------------------------------------------------------

def add_meas_to_obj(obj, meas, band):
    """
    Merge a single-band measurement array into an object catalogue.

    Python port of add_meas2obj.pro.

    For objects with no prior measurement in this band, values are
    copied directly.  For objects that already have a measurement, a
    flux-weighted combined magnitude is computed.

    Parameters
    ----------
    obj : numpy structured array
        Object catalogue (modified in-place).  Must have columns
        ``{band}mag``, ``{band}err``, ``ndet{band}``, ``chi``,
        ``sharp``.
    meas : numpy structured array
        Measurements for one band.  Must have columns:
        ``objid``, ``mag``, ``err``, ``chi``, ``sharp``.
    band : str
        Band name (e.g. ``'g'``).
    """
    b = band.lower()
    magcol = f'{b}mag'
    errcol = f'{b}err'
    nobscol = f'ndet{b}'

    for m in meas:
        oid = m['objid']
        mask = obj['objid'] == oid
        if not mask.any():
            continue
        idx = np.where(mask)[0][0]

        old_mag = float(obj[magcol][idx])
        new_mag = float(m['mag'])
        new_err = float(m['err'])
        old_ndet = int(obj[nobscol][idx])

        if old_ndet == 0 or old_mag >= 50:
            obj[magcol][idx] = new_mag
            obj[errcol][idx] = new_err
            obj[nobscol][idx] = 1
            # Average chi/sharp
            n = int(max(obj.get('ndet', np.array([1]))[idx], 1))
            if 'chi' in obj.dtype.names:
                obj['chi'][idx] = (obj['chi'][idx] * (n - 1) + float(m['chi'])) / n
            if 'sharp' in obj.dtype.names:
                obj['sharp'][idx] = (obj['sharp'][idx] * (n - 1) + float(m['sharp'])) / n
        else:
            # Combine existing mag with new measurement
            old_err = float(obj[errcol][idx])
            cmag, cerr = _flux_weighted_avg(
                np.array([old_mag, new_mag]),
                np.array([old_err, new_err]))
            obj[magcol][idx] = cmag
            obj[errcol][idx] = cerr
            obj[nobscol][idx] = old_ndet + 1
            # Average chi/sharp
            if 'chi' in obj.dtype.names:
                obj['chi'][idx] = (obj['chi'][idx] * old_ndet + float(m['chi'])) / (old_ndet + 1)
            if 'sharp' in obj.dtype.names:
                obj['sharp'][idx] = (obj['sharp'][idx] * old_ndet + float(m['sharp'])) / (old_ndet + 1)


# ---------------------------------------------------------------------------
# Simple measurement averaging (no exposure metadata needed)
# ---------------------------------------------------------------------------

def simple_avg_meas(meas):
    """
    Average all measurements per object into a simplified catalogue.

    Python port of delvered_simpleavgmeas.pro.

    OBJID must already be set in ``meas``.

    Parameters
    ----------
    meas : numpy structured array
        All measurements.  Required fields: objid, exposure, filter,
        mag, err, chi, sharp, ra, dec.  Optional: forced.

    Returns
    -------
    numpy structured array
        Object catalogue with averaged photometry.
    """
    nmeas = len(meas)
    if nmeas == 0:
        return None

    oindex = create_index(meas['objid'])
    nobj = len(oindex['value'])

    # Build measurement → object index
    measobj = np.empty(nmeas, dtype=np.int64)
    for i in range(nobj):
        idx = oindex['index'][oindex['lo'][i]:oindex['hi'][i] + 1]
        measobj[idx] = i

    # Initialise object array
    obj = np.zeros(nobj, dtype=_make_obj_schema_simple())
    obj['objid'] = oindex['value']
    for b in _BANDS:
        obj[f'{b}mag'] = 99.99
        obj[f'{b}err'] = 9.99
        obj[f'{b}rms'] = 99.99
    obj['chi'] = 99.99
    obj['sharp'] = 99.99

    has_forced = 'forced' in meas.dtype.names

    # Build exposure index (strip ccdnum suffix)
    exp_stripped = np.array([e.rsplit('_', 1)[0] for e in meas['exposure']])
    eindex = create_index(exp_stripped)
    nexposure = len(eindex['value'])

    # Build exposure → filter lookup
    exp_filter = {}
    for j in range(nexposure):
        idx0 = eindex['index'][eindex['lo'][j]]
        exp_filter[eindex['value'][j]] = meas['filter'][idx0]

    ufilter = list(set(exp_filter.values()))

    # --- Flux-weighted averaging per filter ---
    for filt in ufilter:
        filt_exps = [ev for ev, fi in exp_filter.items() if fi == filt]
        b = filt.lower()
        magcol = f'{b}mag'
        errcol = f'{b}err'
        rmscol = f'{b}rms'
        nobscol = f'ndet{b}'

        if len(filt_exps) == 1:
            j = np.where(eindex['value'] == filt_exps[0])[0][0]
            mind = eindex['index'][eindex['lo'][j]:eindex['hi'][j] + 1]
            obj[magcol][measobj[mind]] = meas['mag'][mind]
            obj[errcol][measobj[mind]] = meas['err'][mind]
            obj[nobscol][measobj[mind]] = 1
            if has_forced:
                obj['depthflag'][measobj[mind]] |= (meas['forced'][mind].astype(np.int16) + 1)
        else:
            totalwt = np.zeros(nobj)
            totalfluxwt = np.zeros(nobj)
            numobs = np.zeros(nobj, dtype=int)
            for exp_name in filt_exps:
                j = np.where(eindex['value'] == exp_name)[0][0]
                mind = eindex['index'][eindex['lo'][j]:eindex['hi'][j] + 1]
                wt = 1.0 / meas['err'][mind].astype(float)**2
                totalwt[measobj[mind]] += wt
                totalfluxwt[measobj[mind]] += (2.5118864**meas['mag'][mind].astype(float)) * wt
                numobs[measobj[mind]] += 1
                if has_forced:
                    obj['depthflag'][measobj[mind]] |= (meas['forced'][mind].astype(np.int16) + 1)

            good_wt = totalwt > 0
            new_flux = np.where(good_wt, totalfluxwt / np.where(good_wt, totalwt, 1), 0)
            new_mag = np.where(good_wt, 2.5 * np.log10(np.where(new_flux > 0, new_flux, 1)), 99.99)
            new_err = np.where(good_wt, np.sqrt(1.0 / np.where(good_wt, totalwt, 1)), 9.99)
            bad = ~np.isfinite(new_mag)
            new_mag[bad] = 99.99
            new_err[bad] = 9.99

            # Compute RMS
            totaldiff = np.zeros(nobj)
            for exp_name in filt_exps:
                j = np.where(eindex['value'] == exp_name)[0][0]
                mind = eindex['index'][eindex['lo'][j]:eindex['hi'][j] + 1]
                totaldiff[measobj[mind]] += (new_mag[measobj[mind]] - meas['mag'][mind].astype(float))**2
            new_rms = np.sqrt(totaldiff / np.maximum(numobs, 1))
            new_rms[numobs <= 1] = 99.99
            new_rms[bad] = 99.99

            obj[magcol] = new_mag.astype(np.float32)
            obj[errcol] = new_err.astype(np.float32)
            obj[rmscol] = new_rms.astype(np.float32)
            obj[nobscol] = numobs

    # NDET total
    obj['ndet'] = sum(obj[f'ndet{b}'] for b in _BANDS)

    # --- Morphology: average chi, sharp ---
    totchi = np.zeros(nobj)
    numchi = np.zeros(nobj, dtype=int)
    totsharp = np.zeros(nobj)
    numsharp = np.zeros(nobj, dtype=int)

    for j in range(nexposure):
        mind = eindex['index'][eindex['lo'][j]:eindex['hi'][j] + 1]
        chi1 = np.full(nobj, np.nan)
        chi1[measobj[mind]] = meas['chi'][mind].astype(float)
        valid = np.isfinite(chi1) & (chi1 < 1e5)
        totchi[valid] += chi1[valid]
        numchi[valid] += 1

        sharp1 = np.full(nobj, np.nan)
        sharp1[measobj[mind]] = meas['sharp'][mind].astype(float)
        valid = np.isfinite(sharp1) & (sharp1 < 1e5)
        totsharp[valid] += sharp1[valid]
        numsharp[valid] += 1

    good = numchi > 0
    obj['chi'][good] = (totchi[good] / numchi[good]).astype(np.float32)
    good = numsharp > 0
    obj['sharp'][good] = (totsharp[good] / numsharp[good]).astype(np.float32)

    # --- Weighted RA/DEC ---
    totalwt = np.zeros(nobj)
    totalrawt = np.zeros(nobj)
    totaldecwt = np.zeros(nobj)
    for j in range(nexposure):
        mind = eindex['index'][eindex['lo'][j]:eindex['hi'][j] + 1]
        wt = 1.0 / meas['err'][mind].astype(float)**2
        totalwt[measobj[mind]] += wt
        totalrawt[measobj[mind]] += meas['ra'][mind].astype(float) * wt
        totaldecwt[measobj[mind]] += meas['dec'][mind].astype(float) * wt

    good = totalwt > 0
    obj['ra'][good] = (totalrawt[good] / totalwt[good])
    obj['dec'][good] = (totaldecwt[good] / totalwt[good])

    # RA/DEC RMS
    totalradiff = np.zeros(nobj)
    totaldecdiff = np.zeros(nobj)
    for j in range(nexposure):
        mind = eindex['index'][eindex['lo'][j]:eindex['hi'][j] + 1]
        totalradiff[measobj[mind]] += (obj['ra'][measobj[mind]] - meas['ra'][mind].astype(float))**2
        totaldecdiff[measobj[mind]] += (obj['dec'][measobj[mind]] - meas['dec'][mind].astype(float))**2

    ndet_safe = np.maximum(obj['ndet'], 1)
    newrarms = np.sqrt(totalradiff / ndet_safe) * 3600 * np.cos(np.deg2rad(obj['dec']))
    newdecrms = np.sqrt(totaldecdiff / ndet_safe) * 3600
    newrarms[obj['ndet'] == 0] = 99.99
    newdecrms[obj['ndet'] == 0] = 99.99
    newrarms[obj['ndet'] == 1] = 99.99
    newdecrms[obj['ndet'] == 1] = 99.99
    obj['rarms'] = newrarms.astype(np.float32)
    obj['decrms'] = newdecrms.astype(np.float32)

    return obj


# ---------------------------------------------------------------------------
# Full measurement averaging (uses expstr for quality weights)
# ---------------------------------------------------------------------------

def avg_meas(expstr, meas):
    """
    Average measurements for objects using per-exposure quality weights.

    Python port of delvered_avgmeas.pro.

    Computes both "ALL" photometry (all detections) and "BEST"
    photometry (quality-filtered, downweighting shallow non-forced
    detections).

    Parameters
    ----------
    expstr : numpy structured array
        Exposure catalogue.  Required fields: expnum, filter, exptime.
    meas : numpy structured array
        Measurement catalogue.  Required fields: objid, exposure,
        filter, mag, err, chi, sharp, forced, ra, dec.

    Returns
    -------
    numpy structured array
        Object catalogue with full photometry schema.
    """
    nmeas = len(meas)
    if nmeas == 0:
        return None

    print(f"avg_meas: {nmeas} measurements")

    oindex = create_index(meas['objid'])
    nobj_idx = len(oindex['value'])
    print(f"avg_meas: {nobj_idx} objects")

    # Build measurement → object index
    measobj = np.empty(nmeas, dtype=np.int64)
    for i in range(nobj_idx):
        idx = oindex['index'][oindex['lo'][i]:oindex['hi'][i] + 1]
        measobj[idx] = i

    # Initialise object array
    obj = np.zeros(nobj_idx, dtype=_make_obj_schema_full())
    obj['objid'] = oindex['value']
    for b in _BANDS:
        for suf in ('mag', 'err', 'rms', 'magall', 'errall', 'rmsall'):
            col = f'{b}{suf}'
            if 'mag' in suf or 'rms' in suf:
                obj[col] = 99.99
            else:
                obj[col] = 9.99
    obj['chi'] = 99.99
    obj['sharp'] = 99.99
    obj['prob'] = 99.99
    obj['ebv'] = 99.99
    for col in ('asemi', 'bsemi', 'theta', 'ellipticity', 'fwhm',
                'rmsvar', 'madvar', 'iqrvar', 'etavar', 'jvar',
                'kvar', 'chivar', 'romsvar', 'nsigvar', 'gaia_xdist',
                'gaia_ra', 'gaia_dec', 'gaia_parallax', 'gaia_parallax_error',
                'gaia_pmra', 'gaia_pmra_error', 'gaia_pmdec', 'gaia_pmdec_error',
                'gaia_gmag', 'gaia_gmag_error', 'gaia_bpmag', 'gaia_bpmag_error',
                'gaia_rpmag', 'gaia_rpmag_error', 'mag_auto', 'magerr_auto'):
        if col in obj.dtype.names:
            obj[col] = 999999.0
    obj['variable10sig'] = -1

    # Build exposure index from meas
    dum = np.array([e.rsplit('_', 1)[0] for e in meas['exposure']])
    exp_stripped = np.array([e.split('-', 1)[-1] if '-' in e else e for e in dum])
    eindex = create_index(exp_stripped)
    nexposure = len(eindex['value'])

    # Match expstr to exposure indices
    expstr_num = expstr['expnum'].astype(str)
    exp_order = eindex['value']
    # Reorder expstr to match eindex
    exp_lookup = {str(en): i for i, en in enumerate(expstr_num)}
    mexpstr_idx = [exp_lookup.get(str(ev)) for ev in exp_order]
    mexpstr_valid = [i for i in mexpstr_idx if i is not None]
    mexpstr = expstr[mexpstr_valid] if mexpstr_valid else expstr[:0]

    # Helper: flux-weighted average across exposures
    def _avg_filter(filt, use_best=False):
        b = filt.lower()
        filtmask = np.array([mexpstr[k]['filter'] if k < len(mexpstr) else '' for k in range(nexposure)]) == filt
        filtind = np.where(filtmask)[0]

        if len(filtind) == 0:
            return

        magcol = f'{b}magall'
        errcol = f'{b}errall'
        rmscol = f'{b}rmsall'
        nobscol = f'ndetall{b}'
        if use_best:
            magcol = f'{b}mag'
            errcol = f'{b}err'
            rmscol = f'{b}rms'
            nobscol = f'ndet{b}'

        if len(filtind) == 1:
            k = filtind[0]
            mind = eindex['index'][eindex['lo'][k]:eindex['hi'][k] + 1]
            obj[magcol][measobj[mind]] = meas['mag'][mind].astype(float)
            obj[errcol][measobj[mind]] = meas['err'][mind].astype(float)
            obj[nobscol][measobj[mind]] = 1
            obj['depthflag'][measobj[mind]] |= (meas['forced'][mind].astype(np.int16) + 1)
        else:
            totalwt = np.zeros(nobj_idx)
            totalfluxwt = np.zeros(nobj_idx)
            numobs = np.zeros(nobj_idx, dtype=int)
            numgdobs = np.zeros(nobj_idx, dtype=int)

            for k in filtind:
                mind = eindex['index'][eindex['lo'][k]:eindex['hi'][k] + 1]
                wt = 1.0 / meas['err'][mind].astype(float)**2
                if use_best and k < len(mexpstr) and mexpstr[k]['exptime'] < 90:
                    shallow_bad = (meas['forced'][mind] == 0) & (meas['err'][mind].astype(float) > 0.2)
                    wt[shallow_bad] *= 1e-4
                    good_idx = np.where(~shallow_bad)[0]
                    numgdobs[measobj[mind[good_idx]]] += 1
                else:
                    numgdobs[measobj[mind]] += 1

                totalwt[measobj[mind]] += wt
                totalfluxwt[measobj[mind]] += (2.5118864**meas['mag'][mind].astype(float)) * wt
                numobs[measobj[mind]] += 1
                obj['depthflag'][measobj[mind]] |= (meas['forced'][mind].astype(np.int16) + 1)

            good_wt = totalwt > 0
            new_flux = np.where(good_wt, totalfluxwt / np.where(good_wt, totalwt, 1), 0)
            new_mag = np.where(good_wt & (new_flux > 0), 2.5 * np.log10(np.where(new_flux > 0, new_flux, 1)), 99.99)
            new_err = np.where(good_wt, np.sqrt(1.0 / np.where(good_wt, totalwt, 1)), 9.99)
            bad = ~np.isfinite(new_mag)
            new_mag[bad] = 99.99
            new_err[bad] = 9.99

            if use_best:
                obj[nobscol] = numgdobs
                poor = (numgdobs == 0) & (new_mag < 50)
                new_mag[poor] = 99.99
                new_err[poor] = 9.99

            totaldiff = np.zeros(nobj_idx)
            numobs2 = np.zeros(nobj_idx, dtype=int)
            for k in filtind:
                mind = eindex['index'][eindex['lo'][k]:eindex['hi'][k] + 1]
                if use_best and k < len(mexpstr) and mexpstr[k]['exptime'] < 90:
                    good_m = np.where(~((meas['forced'][mind] == 0) & (meas['err'][mind].astype(float) > 0.2)))[0]
                else:
                    good_m = np.arange(len(mind))
                if len(good_m):
                    mi = mind[good_m]
                    totaldiff[measobj[mi]] += (new_mag[measobj[mi]] - meas['mag'][mi].astype(float))**2
                    numobs2[measobj[mi]] += 1

            new_rms = np.sqrt(totaldiff / np.maximum(numobs2, 1))
            new_rms[numobs2 <= 1] = 99.99
            new_rms[bad] = 99.99

            obj[magcol] = new_mag.astype(np.float32)
            obj[errcol] = new_err.astype(np.float32)
            obj[rmscol] = new_rms.astype(np.float32)
            if not use_best:
                obj[nobscol] = numobs

    # Unique filters present in mexpstr
    if len(mexpstr):
        ufilter = list(set(mexpstr['filter']))
    else:
        ufilter = list(set(meas['filter']))

    print(f"avg_meas: {nexposure} exposures, {len(ufilter)} unique filters")

    # Compute ALL photometry
    for filt in ufilter:
        _avg_filter(filt, use_best=False)

    # Compute BEST photometry
    for filt in ufilter:
        _avg_filter(filt, use_best=True)

    # NDET totals
    obj['ndet'] = sum(obj[f'ndet{b}'] for b in _BANDS)
    obj['ndetall'] = sum(obj[f'ndetall{b}'] for b in _BANDS)

    # --- Morphology ---
    totchi = np.zeros(nobj_idx)
    numchi = np.zeros(nobj_idx, dtype=int)
    totsharp = np.zeros(nobj_idx)
    totwtsharp = np.zeros(nobj_idx)

    for j in range(nexposure):
        mind = eindex['index'][eindex['lo'][j]:eindex['hi'][j] + 1]

        chi1 = np.full(nobj_idx, np.nan)
        chi1[measobj[mind]] = meas['chi'][mind].astype(float)
        valid = np.isfinite(chi1) & (chi1 < 1e5)
        totchi[valid] += chi1[valid]
        numchi[valid] += 1

        sharp1 = np.full(nobj_idx, np.nan)
        sharp1[measobj[mind]] = meas['sharp'][mind].astype(float)
        snr = np.zeros(nobj_idx)
        snr[measobj[mind]] = 1.087 / meas['err'][mind].astype(float)
        forced_arr = np.full(nobj_idx, -1, dtype=int)
        forced_arr[measobj[mind]] = meas['forced'][mind].astype(int)

        valid = np.isfinite(sharp1) & (sharp1 < 1e5)
        wt_sharp = np.where(valid & (snr < 5) & (forced_arr != 1), 1e-4, 1.0) * valid
        totsharp += wt_sharp * np.where(valid, sharp1, 0)
        totwtsharp += wt_sharp

    good = numchi > 0
    obj['chi'][good] = (totchi[good] / numchi[good]).astype(np.float32)
    good = totwtsharp > 0
    obj['sharp'][good] = (totsharp[good] / totwtsharp[good]).astype(np.float32)

    # --- Weighted RA/DEC ---
    totalwt = np.zeros(nobj_idx)
    totalrawt = np.zeros(nobj_idx)
    totaldecwt = np.zeros(nobj_idx)
    for j in range(nexposure):
        mind = eindex['index'][eindex['lo'][j]:eindex['hi'][j] + 1]
        wt = 1.0 / meas['err'][mind].astype(float)**2
        totalwt[measobj[mind]] += wt
        totalrawt[measobj[mind]] += meas['ra'][mind].astype(float) * wt
        totaldecwt[measobj[mind]] += meas['dec'][mind].astype(float) * wt

    good = totalwt > 0
    obj['ra'][good] = (totalrawt[good] / totalwt[good])
    obj['dec'][good] = (totaldecwt[good] / totalwt[good])

    # RA/DEC RMS
    totalradiff = np.zeros(nobj_idx)
    totaldecdiff = np.zeros(nobj_idx)
    for j in range(nexposure):
        mind = eindex['index'][eindex['lo'][j]:eindex['hi'][j] + 1]
        totalradiff[measobj[mind]] += (obj['ra'][measobj[mind]] - meas['ra'][mind].astype(float))**2
        totaldecdiff[measobj[mind]] += (obj['dec'][measobj[mind]] - meas['dec'][mind].astype(float))**2

    ndetall_safe = np.maximum(obj['ndetall'], 1)
    newrarms = np.sqrt(totalradiff / ndetall_safe) * 3600 * np.cos(np.deg2rad(obj['dec']))
    newdecrms = np.sqrt(totaldecdiff / ndetall_safe) * 3600
    newrarms[obj['ndetall'] == 0] = 99.99
    newdecrms[obj['ndetall'] == 0] = 99.99
    newrarms[obj['ndetall'] == 1] = 99.99
    newdecrms[obj['ndetall'] == 1] = 99.99
    obj['rarms'] = newrarms.astype(np.float32)
    obj['decrms'] = newdecrms.astype(np.float32)

    # --- EBV from SFD ---
    try:
        from .extinction import _sfd_ebv
        obj['ebv'] = _sfd_ebv(obj['ra'], obj['dec']).astype(np.float32)
    except Exception:
        pass

    return obj


# ---------------------------------------------------------------------------
# Loop wrapper: avg_meas across multiple bricks/fields
# ---------------------------------------------------------------------------

def avg_meas_loop(brick_list, data_dir, output_dir=None, **kwargs):
    """
    Run avg_meas for a list of bricks.

    Python port of delvered_avgmeas_loop.pro.

    Parameters
    ----------
    brick_list : list of str
        Brick names (e.g. ``['1234p567', ...]``).
    data_dir : str or Path
        Root directory containing per-brick subdirectories.
    output_dir : str or Path, optional
        Where to write per-brick object catalogues.  Defaults to
        ``data_dir``.
    **kwargs
        Additional keyword arguments passed to ``avg_meas``.

    Returns
    -------
    None
        Results are written to FITS files under ``output_dir``.
    """
    import os
    from pathlib import Path
    from astropy.io import fits

    data_dir = Path(data_dir)
    if output_dir is None:
        output_dir = data_dir

    for brick in brick_list:
        brick_dir = data_dir / brick
        meas_file = brick_dir / f'{brick}_meas.fits.gz'
        expstr_file = brick_dir / f'{brick}_expstr.fits.gz'
        if not meas_file.exists() or not expstr_file.exists():
            print(f"avg_meas_loop: skipping {brick} (files not found)")
            continue

        meas = fits.getdata(str(meas_file), 1)
        expstr = fits.getdata(str(expstr_file), 1)

        obj = avg_meas(expstr, meas, **kwargs)
        if obj is None:
            continue

        out_file = Path(output_dir) / brick / f'{brick}_object.fits'
        out_file.parent.mkdir(parents=True, exist_ok=True)
        fits.writeto(str(out_file), obj, overwrite=True)
        print(f"avg_meas_loop: wrote {out_file}")
