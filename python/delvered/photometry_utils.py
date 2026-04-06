"""
Photometric utility functions for the DELVERED pipeline.

Python port of:
  - delvered_getteff.pro     → get_teff()
  - delvered_getmodelmag.pro → get_model_mag()
  - delvered_photvar.pro     → photvar()
  - delvered_deltasnrsort.pro → deltasnr_sort()
"""

import re
import numpy as np
from scipy.ndimage import uniform_filter1d
from scipy.interpolate import interp1d


# ---------------------------------------------------------------------------
# Effective exposure time
# ---------------------------------------------------------------------------

# Fiducial band parameters (Neilsen+2015)
_UFILTERS = ['u', 'g', 'r', 'i', 'z', 'Y']
_MED_ZPTERM = np.array([2.03, -0.07, -0.39, -0.33, -0.03, 1.08])
_MED_BACKGROUND = np.array([0.10, 1.29, 3.77, 13.19, 23.46, 17.02])
_FWHM_FIDUCIAL = 0.9   # arcsec
_BACKGROUND_FIDUCIAL = 1.0


def get_teff(chstr):
    """
    Calculate effective exposure time (teff) for each chip.

    Python port of delvered_getteff.pro (Neilsen+2015 formalism).

    Parameters
    ----------
    chstr : numpy structured array or dict-like
        Chip structure with fields: filter, gain, calib_zpterm, fwhm,
        skymode, exptime.

    Returns
    -------
    numpy structured array
        Input chstr extended with new fields: eta, background, tau, teff.
    """
    from numpy.lib.recfunctions import append_fields

    nchips = len(chstr)
    zpterm = chstr['calib_zpterm'].copy().astype(float)

    # Old DECam cameras correction (gain < 2)
    old_cam = chstr['gain'] < 2
    if old_cam.any():
        zpterm[old_cam] += 1.55

    # Per-band delta zero-point
    delta_zpterm = np.zeros(nchips, dtype=float)
    background_fiducial = np.zeros(nchips, dtype=float)
    for f, fname in enumerate(_UFILTERS):
        mask = chstr['filter'] == fname
        if mask.any():
            delta_zpterm[mask] = zpterm[mask] - _MED_ZPTERM[f]
            background_fiducial[mask] = _MED_BACKGROUND[f]

    eta = np.exp(-delta_zpterm)
    fwhm_arcsec = chstr['fwhm'].astype(float) * 0.235
    background = chstr['skymode'].astype(float) / chstr['exptime'].astype(float)

    tau = (eta**2) * (fwhm_arcsec / _FWHM_FIDUCIAL)**(-2) * (background / _BACKGROUND_FIDUCIAL)**(-1)
    teff = tau * chstr['exptime'].astype(float)

    # Append new fields to the structured array
    for name, data in [('eta', eta), ('background', background),
                       ('tau', tau), ('teff', teff)]:
        if name not in chstr.dtype.names:
            chstr = append_fields(chstr, name, data, usemask=False)
        else:
            chstr[name] = data

    return chstr


# ---------------------------------------------------------------------------
# teff map (low-resolution 36×36 from chip footprints)
# ---------------------------------------------------------------------------

def teff_map(chstr):
    """
    Build a low-resolution (36×36) teff map from a set of chips.

    Python port of delvered_teffmap.pro.

    Parameters
    ----------
    chstr : numpy structured array
        Chip structure.  Must have VX and VY arrays (pixel vertices)
        and TEFF field.

    Returns
    -------
    numpy ndarray, shape (36, 36)
        Effective exposure time map.
    """
    from skimage.draw import polygon as sk_polygon

    tilenx = 3600
    tileny = 3600
    nx_lr = tilenx // 10
    ny_lr = tileny // 10

    teffim = np.zeros((ny_lr, nx_lr), dtype=np.float32)

    for i in range(len(chstr)):
        teffim1 = np.zeros((ny_lr, nx_lr), dtype=np.float32)
        vx = np.asarray(chstr['vx'][i]) / 10.0
        vy = np.asarray(chstr['vy'][i]) / 10.0
        rr, cc = sk_polygon(vy, vx, shape=(ny_lr, nx_lr))
        if len(rr):
            teffim1[rr, cc] = chstr['teff'][i]
        teffim += teffim1

    # Rebin to 36×36
    factor = nx_lr // 36
    teffim = teffim.reshape(36, factor, 36, factor).mean(axis=(1, 3))
    return teffim


# ---------------------------------------------------------------------------
# Delta-S/N sorting of exposures
# ---------------------------------------------------------------------------

def deltasnr_sort(expstr):
    """
    Sort exposures by their incremental S/N contribution.

    Python port of delvered_deltasnrsort.pro (Neilsen+2015).

    The algorithm iteratively sorts exposures so that new coverage is
    preferred over redundant coverage.  Convergence is reached when
    the total delta-S/N no longer changes or after 10 iterations.

    Parameters
    ----------
    expstr : numpy structured array
        Exposure structure.  Must have fields ``teff`` (scalar per exposure)
        and ``teffmap`` (36×36 array per exposure).  Fields ``deltasnr``
        and ``fraccovered`` are added/updated in-place.

    Returns
    -------
    numpy structured array
        Exposure structure sorted by delta-S/N, with updated
        ``deltasnr`` and ``fraccovered`` fields.
    """
    from numpy.lib.recfunctions import append_fields

    # Add deltasnr / fraccovered fields if absent
    for fname, init in [('deltasnr', 0.0), ('fraccovered', 0.0)]:
        if fname not in expstr.dtype.names:
            expstr = append_fields(expstr, fname, np.zeros(len(expstr)), usemask=False)

    nexp = len(expstr)
    nmap = 36 * 36

    totdeltasnr_last = -1
    for count in range(11):
        if count == 0:
            si = np.argsort(expstr['teff'])[::-1]
        else:
            si = np.argsort(expstr['deltasnr'])[::-1]
        expstr = expstr[si]

        sumteffim = np.zeros((36, 36), dtype=float)
        for j in range(nexp):
            snr_old = np.sqrt(sumteffim)
            snr_new = np.sqrt(sumteffim + expstr['teffmap'][j])
            expstr['deltasnr'][j] = np.sum(snr_new - snr_old)
            sumteffim += expstr['teffmap'][j]
            expstr['fraccovered'][j] = np.sum(sumteffim > 0) / nmap

        totdeltasnr = int(np.sum(expstr['deltasnr']))
        if totdeltasnr == totdeltasnr_last or count >= 10:
            break
        totdeltasnr_last = totdeltasnr

    return expstr


# ---------------------------------------------------------------------------
# Model magnitudes from equation file
# ---------------------------------------------------------------------------

def get_model_mag(cat, instfilt, dec, eqnfile):
    """
    Calculate model magnitudes for reference catalogue sources.

    Python port of delvered_getmodelmag.pro.

    Reads a text equation file and applies colour and magnitude
    transformation equations to the input catalogue.

    Parameters
    ----------
    cat : numpy structured array or dict with numpy arrays
        Source catalogue with magnitude columns (GMAG, JMAG, KMAG,
        E_GMAG, etc.).
    instfilt : str
        Instrument-filter string, e.g. ``'c4d-g'``.
    dec : float
        Declination of the exposure in degrees.
    eqnfile : str
        Path to the model magnitude equation file.

    Returns
    -------
    numpy ndarray, shape (N, 3)
        Columns: [model_mag, model_mag_err, color].
        Sources that do not pass quality cuts have model_mag=99.99.
    """
    import os
    import pandas as pd

    ncat = len(cat) if hasattr(cat, '__len__') else 1
    result = np.full((ncat, 3), [99.99, 99.90, 0.0], dtype=float)

    if not os.path.exists(eqnfile):
        print(f"get_model_mag: {eqnfile} NOT FOUND")
        return result

    # Load equation file
    eqndf = pd.read_csv(eqnfile, delim_whitespace=True)
    eqndf.columns = [c.lower() for c in eqndf.columns]

    # Parse colour and DEC ranges stored as "[lo,hi]" strings
    def parse_range(s):
        s = str(s).strip('[]')
        lo, hi = s.split(',')
        return float(lo), float(hi)

    eqndf['colorlim_lo'], eqndf['colorlim_hi'] = zip(*eqndf['colorange'].map(parse_range))
    eqndf['declim_lo'], eqndf['declim_hi'] = zip(*eqndf['decrange'].map(parse_range))

    # Find matching equation
    inst, band = instfilt.rsplit('-', 1)
    mask = ((eqndf['instrument'].str.strip() == inst) &
            (eqndf['band'].str.strip() == band) &
            (eqndf['declim_lo'] <= dec) & (dec <= eqndf['declim_hi']))
    matches = eqndf[mask]
    if matches.empty:
        print(f"get_model_mag: no equation for {instfilt} at dec={dec:.2f}")
        return result
    eqn = matches.iloc[0]

    # Helper: get catalogue column as array (case-insensitive)
    cat_colnames_upper = {k.upper(): k for k in (cat.dtype.names if hasattr(cat, 'dtype') else cat.keys())}

    def get_col(colname):
        key = cat_colnames_upper.get(colname.upper())
        if key is None:
            raise KeyError(f"Column {colname} not found in catalogue")
        return np.asarray(cat[key], dtype=float)

    # Build a local namespace for equation evaluation
    ns = {}
    coloreqn = str(eqn.get('coloreqn', ''))
    modelmageqn = str(eqn.get('modelmageqn', ''))
    colorlim = (eqn['colorlim_lo'], eqn['colorlim_hi'])
    use_color = not (colorlim[0] < -10 and colorlim[1] > 10 and 'COLOR' not in modelmageqn.upper())

    # Identify columns referenced in equations
    def extract_colnames(eqn_str):
        tokens = re.split(r'[-+*^/\s]', eqn_str)
        cols = []
        for t in tokens:
            t = t.strip()
            if t and not _is_number(t) and t.upper() not in ('COLOR', ''):
                cols.append(t.upper())
        return list(set(cols))

    # Load required columns
    required = extract_colnames(modelmageqn)
    if use_color:
        required += extract_colnames(coloreqn)
    required = list(set(required))

    for colname in required:
        try:
            ns[colname] = get_col(colname)
        except KeyError:
            print(f"get_model_mag: missing column {colname}")
            return result

    # Determine good sources (magnitudes < 50, > 0, finite)
    mag_cols = [c for c in required if c.endswith('MAG') and not c.startswith('E_')]
    goodmask = np.ones(ncat, dtype=bool)
    for mc in mag_cols:
        v = ns[mc]
        goodmask &= (v < 50) & (v > 0) & np.isfinite(v)

    # Apply quality cuts from equation file
    qualcuts = str(eqn.get('qualitycuts', '')).strip()
    if qualcuts and qualcuts not in ('', 'None', 'nan'):
        qualcuts_py = (qualcuts
                       .replace('<=', ' <= ').replace('>=', ' >= ')
                       .replace('<', ' < ').replace('>', ' > ')
                       .replace('=', ' == ').replace('&', ' and ')
                       .replace('|', ' or '))
        for colname in required:
            qualcuts_py = qualcuts_py.replace(colname, f"ns['{colname}']")
        try:
            goodmask &= eval(qualcuts_py)  # safe: ns contains only numpy arrays
        except Exception as e:
            print(f"get_model_mag: could not apply quality cuts: {e}")

    # Compute colour
    if use_color:
        color_expr = _replace_cols(coloreqn, ns)
        color = eval(color_expr)  # noqa: S307
        goodmask &= (color >= colorlim[0]) & (color <= colorlim[1])
    else:
        color = np.zeros(ncat, dtype=float)

    ns['COLOR'] = color
    gd = np.where(goodmask)[0]
    if len(gd) == 0:
        print("get_model_mag: no good sources left")
        return result

    # Compute model magnitudes for good sources
    ns_gd = {k: v[gd] if isinstance(v, np.ndarray) else v for k, v in ns.items()}
    mag_expr = _replace_cols(modelmageqn, ns_gd)
    try:
        modelmag_gd = eval(mag_expr)  # noqa: S307
        result[gd, 0] = modelmag_gd
        result[gd, 2] = color[gd]
    except Exception as e:
        print(f"get_model_mag: error evaluating model magnitude equation: {e}")

    # Propagate errors in quadrature (simplified: sum of E_ columns in quadrature)
    err_cols = [c for c in cat_colnames_upper if c.startswith('E_') and c[2:] in required]
    ns_err = {}
    for ec in err_cols:
        try:
            v = get_col(ec)
            v = np.where((v > 10) | ~np.isfinite(v), 9.99, v)
            ns_err[ec] = v
        except KeyError:
            pass

    if ns_err:
        sumsq = np.zeros(ncat, dtype=float)
        for ec, v in ns_err.items():
            sumsq += v**2
        modelmagerr = np.sqrt(sumsq)
        modelmagerr[~goodmask] = 99.90
        result[:, 1] = modelmagerr

    return result


def _is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def _replace_cols(eqn_str, ns):
    """Replace column names in an equation string with ns['COLNAME']."""
    result = eqn_str.upper()
    # Sort by length descending to avoid partial replacements
    for colname in sorted(ns.keys(), key=len, reverse=True):
        if colname != 'COLOR':
            result = re.sub(r'\b' + colname + r'\b', f"ns['{colname}']", result)
    result = result.replace('COLOR', "ns['COLOR']")
    return result


# ---------------------------------------------------------------------------
# Photometric variability indices
# ---------------------------------------------------------------------------

def photvar(meas, obj):
    """
    Calculate photometric variability indices for each object.

    Python port of delvered_photvar.pro.

    Updates ``obj`` in-place with the columns: rmsvar, madvar, iqrvar,
    etavar, jvar, kvar, chivar, romsvar, variable10sig, nsigvar.

    Parameters
    ----------
    meas : numpy structured array
        All measurements with fields: objid, filter, mag, err, mjd.
    obj : numpy structured array
        Object catalogue.  Must have per-filter MAG/ERR columns and
        variability metric columns.
    """
    from .utils import create_index

    nmeas = len(meas)
    nobj = len(obj)
    if nmeas == 0 or nobj == 0:
        return

    print(f"photvar: {nmeas} measurements")

    # Index measurements by objid
    oindex = create_index(meas['objid'])
    nobj_idx = len(oindex['value'])
    print(f"photvar: {nobj_idx} objects")

    # Map obj.objid → position in obj array
    obj_id_map = {oid: i for i, oid in enumerate(obj['objid'])}

    obj_fields = obj.dtype.names
    nobj_all = len(obj)

    # Initialise variability metrics to NaN
    for col in ('rmsvar', 'madvar', 'iqrvar', 'etavar', 'jvar',
                'kvar', 'chivar', 'romsvar'):
        if col in obj_fields:
            obj[col] = np.nan

    fidmag = np.full(nobj_all, np.nan)

    for i in range(nobj_idx):
        ind = oindex['index'][oindex['lo'][i]:oindex['hi'][i] + 1]
        meas1 = meas[ind]
        objid = oindex['value'][i]
        oi = obj_id_map.get(objid)
        if oi is None:
            continue
        obj1 = obj[oi]

        filt_idx = create_index(meas1['filter'])
        nf = len(filt_idx['value'])
        nm1 = len(meas1)
        resid = np.full(nm1, np.nan)
        relresid = np.full(nm1, np.nan)

        for f in range(nf):
            filt = str(filt_idx['value'][f]).upper()
            findx = filt_idx['index'][filt_idx['lo'][f]:filt_idx['hi'][f] + 1]
            magcol = filt + 'MAG'
            errcol = filt + 'ERR'
            if magcol not in obj_fields or errcol not in obj_fields:
                continue
            obj_mag = float(obj1[magcol])
            mags1 = meas1[findx]['mag'].astype(float)
            errs1 = np.maximum(meas1[findx]['err'].astype(float), 0.02)
            gph = np.where(mags1 < 50)[0]
            if len(gph) > 1:
                n_gph = len(gph)
                resid[findx[gph]] = mags1[gph] - obj_mag
                relresid[findx[gph]] = np.sqrt(n_gph / (n_gph - 1.0)) * (mags1[gph] - obj_mag) / errs1[gph]

        # Absolute residual indices
        gdresid = np.where(np.isfinite(resid))[0]
        if len(gdresid) > 0:
            r2 = resid[gdresid]
            n_r = len(r2)
            sumresidsq = np.sum(r2**2)
            tsi = np.argsort(meas1[gdresid]['mjd'].astype(float))
            r2tsi = r2[tsi]
            q25, q50, q75 = np.percentile(r2, [25, 50, 75])
            obj['rmsvar'][oi] = np.sqrt(sumresidsq / n_r)
            obj['madvar'][oi] = 1.4826 * np.median(np.abs(r2 - q50))
            obj['iqrvar'][oi] = 0.741289 * (q75 - q25)
            denom_eta = np.sum((r2tsi[1:] - r2tsi[:-1])**2)
            obj['etavar'][oi] = sumresidsq / denom_eta if denom_eta != 0 else np.nan

        # Relative residual indices
        gdrel = np.where(np.isfinite(relresid))[0]
        if len(gdrel) > 0:
            rr2 = relresid[gdrel]
            n_rr = len(rr2)
            pk = rr2**2 - 1
            obj['jvar'][oi] = np.sum(np.sign(pk) * np.sqrt(np.abs(pk))) / n_rr
            obj['chivar'][oi] = np.sqrt(np.sum(rr2**2)) / n_rr
            kdenom = np.sqrt(np.sum(rr2**2) / n_rr)
            obj['kvar'][oi] = (np.sum(np.abs(rr2)) / n_rr) / kdenom if kdenom != 0 else np.nan
            obj['romsvar'][oi] = np.sum(np.abs(rr2)) / (n_rr - 1) if n_rr > 1 else np.nan

        # Fiducial magnitude (priority: r, g, i, z, Y, u)
        for band in ('r', 'g', 'i', 'z', 'y', 'u'):
            mc = band.upper() + 'MAG'
            if mc in obj_fields:
                v = float(obj1[mc])
                if v < 50:
                    fidmag[oi] = v
                    break

    # Detect variables using MAD variability index
    if 'variable10sig' in obj_fields:
        obj['variable10sig'] = 0
    if 'nsigvar' in obj_fields:
        obj['nsigvar'] = np.nan

    gdvar = np.where(np.isfinite(obj['madvar']) & np.isfinite(fidmag))[0]
    if len(gdvar) < 2:
        print("photvar: not enough good MADVAR values to detect variables")
        return

    # Bin by magnitude and smooth to get median MAD vs magnitude
    binsize = 0.25
    fmag_gd = fidmag[gdvar]
    madvar_gd = obj['madvar'][gdvar]

    bins = np.arange(fmag_gd.min(), fmag_gd.max() + binsize, binsize)
    bin_idx = np.digitize(fmag_gd, bins)
    nbins = len(bins)
    varmed = np.full(nbins, np.nan)
    fidmagmed = np.full(nbins, np.nan)
    numhist = np.zeros(nbins, dtype=int)
    for b in range(nbins):
        mask = bin_idx == b
        if mask.sum() > 0:
            varmed[b] = np.nanmedian(madvar_gd[mask])
            fidmagmed[b] = np.nanmedian(fmag_gd[mask])
            numhist[b] = mask.sum()

    # Smooth and interpolate
    gv = np.where(np.isfinite(varmed))[0]
    if len(gv) <= 1:
        print("photvar: not enough good MADVAR bins")
        return

    smlen = 5
    smvarmed = _smooth1d(varmed[gv], smlen)
    interp_fn = interp1d(fidmagmed[gv], smvarmed, bounds_error=False,
                         fill_value=(smvarmed[0], smvarmed[-1]))
    objvarmed = interp_fn(fidmag)
    objvarmed = np.maximum(objvarmed, np.min(smvarmed))

    # MAD of residuals from the median line
    varsig = np.full(nbins, np.nan)
    for b in range(nbins):
        mask = bin_idx == b
        if mask.sum() >= 3:
            varsig[b] = np.nanmedian(np.abs(madvar_gd[mask] - interp_fn(fmag_gd[mask]))) * 1.4826

    gv2 = np.where(np.isfinite(varsig))[0]
    if len(gv2) <= 1:
        return
    smvarsig = _smooth1d(varsig[gv2], smlen)
    interp_sig = interp1d(fidmagmed[gv2], smvarsig, bounds_error=False,
                          fill_value=(smvarsig[0], smvarsig[-1]))
    objvarsig = interp_sig(fidmag)
    objvarsig = np.maximum(objvarsig, np.min(smvarsig))

    # Flag variables at 10-sigma
    nsigvar = (obj['madvar'] - objvarmed) / objvarsig
    if 'nsigvar' in obj_fields:
        obj['nsigvar'][gdvar] = nsigvar[gdvar]
    isvar = gdvar[nsigvar[gdvar] > 10.0]
    n_isvar = len(isvar)
    print(f"photvar: {n_isvar} variables detected")
    if n_isvar > 0 and 'variable10sig' in obj_fields:
        obj['variable10sig'][isvar] = 1


def _smooth1d(arr, width):
    """Simple boxcar smoothing, ignoring NaNs."""
    out = np.full_like(arr, np.nan, dtype=float)
    arr = np.asarray(arr, dtype=float)
    for i in range(len(arr)):
        lo = max(0, i - width // 2)
        hi = min(len(arr), i + width // 2 + 1)
        vals = arr[lo:hi]
        good = vals[np.isfinite(vals)]
        if len(good):
            out[i] = np.mean(good)
    return out
