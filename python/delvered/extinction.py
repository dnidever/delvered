"""
Galactic extinction utilities for the DELVERED pipeline.

Python port of:
  - delvered_getexttype.pro  → get_ext_type()
  - delvered_getreddening.pro → get_reddening()
"""

import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _sfd_ebv(ra, dec):
    """
    Query the SFD98 E(B-V) dust map.

    Requires the ``dustmaps`` package with SFD map fetched.

    Parameters
    ----------
    ra, dec : array-like
        Equatorial coordinates (degrees, J2000).

    Returns
    -------
    numpy array
        E(B-V) values.
    """
    try:
        from dustmaps.sfd import SFDQuery
        coords = SkyCoord(ra=np.atleast_1d(ra) * u.deg,
                          dec=np.atleast_1d(dec) * u.deg, frame='icrs')
        sfd = SFDQuery()
        return sfd(coords)
    except ImportError:
        raise ImportError(
            "The 'dustmaps' package is required for SFD queries.  "
            "Install with: pip install dustmaps && python -c \"from dustmaps.sfd import fetch; fetch()\""
        )


def _gal_coords(ra, dec):
    """Convert equatorial to Galactic coordinates."""
    c = SkyCoord(ra=np.atleast_1d(ra) * u.deg,
                 dec=np.atleast_1d(dec) * u.deg, frame='icrs')
    return c.galactic.l.deg, c.galactic.b.deg


def _sphdist(ra1, dec1, ra2, dec2):
    """Angular separation in degrees between two sky positions."""
    c1 = SkyCoord(ra=ra1 * u.deg, dec=dec1 * u.deg, frame='icrs')
    c2 = SkyCoord(ra=ra2 * u.deg, dec=dec2 * u.deg, frame='icrs')
    return c1.separation(c2).deg


# ---------------------------------------------------------------------------
# Extinction type selection
# ---------------------------------------------------------------------------

def get_ext_type(cenra, cendec, radius):
    """
    Determine the appropriate dust extinction method for a field.

    Python port of delvered_getexttype.pro.

    Extinction types:
      1 – SFD (default, high Galactic latitude, low dust)
      2 – RJCE ALLWISE (moderate dust, no GLIMPSE/SAGE)
      3 – RJCE GLIMPSE (Galactic plane with Spitzer GLIMPSE data)
      4 – RJCE SAGE    (LMC/SMC with Spitzer SAGE data)

    Parameters
    ----------
    cenra : float
        Central RA of the field (degrees).
    cendec : float
        Central DEC of the field (degrees).
    radius : float
        Field radius (degrees).

    Returns
    -------
    int
        Extinction type (1–4).
    """
    try:
        from astroquery.vizier import Vizier
    except ImportError:
        Vizier = None
        print("get_ext_type: astroquery not available; GLIMPSE/SAGE checks skipped")

    cengl, cengb = _gal_coords(cenra, cendec)
    cengl = float(cengl)
    cengb = float(cengb)

    # Build a 10×10 grid across the field to sample SFD E(B-V)
    x = np.linspace(-radius, radius, 100)
    xx, yy = np.meshgrid(x, x)
    # Gnomonic reverse projection → RA/DEC
    rr, dd = _gnomonic_reverse(xx.ravel(), yy.ravel(), cenra, cendec)
    gl_grid, gb_grid = _gal_coords(rr, dd)
    ebv_grid = _sfd_ebv(rr, dd)
    maxebv = float(np.max(ebv_grid))

    ext_type = 0

    # Check for GLIMPSE data (Galactic plane, |b| < 5, 65 < l < 290)
    if abs(cengb) < 5 and (cengl < 65 or cengl > 290):
        if Vizier is not None:
            x3 = np.linspace(-radius, radius, 3)
            xx3, yy3 = np.meshgrid(x3, x3)
            rr3, dd3 = _gnomonic_reverse(xx3.ravel(), yy3.ravel(), cenra, cendec)
            for ra_q, dec_q in zip(rr3, dd3):
                try:
                    Vizier.ROW_LIMIT = 1
                    result = Vizier.query_region(
                        SkyCoord(ra=ra_q * u.deg, dec=dec_q * u.deg, frame='icrs'),
                        radius=1.0 * u.arcmin,
                        catalog='II/293/glimpse')
                    if result:
                        ext_type = 3
                        break
                except Exception:
                    pass

    # Check for SAGE data (near LMC/SMC)
    lmcrad = _sphdist(81.9, -69.867, cenra, cendec)
    smcrad = _sphdist(13.183, -72.8283, cenra, cendec)
    if lmcrad < 5.0 or smcrad < 4.0:
        if Vizier is not None:
            x3 = np.linspace(-radius, radius, 3)
            xx3, yy3 = np.meshgrid(x3, x3)
            rr3, dd3 = _gnomonic_reverse(xx3.ravel(), yy3.ravel(), cenra, cendec)
            for ra_q, dec_q in zip(rr3, dd3):
                try:
                    Vizier.ROW_LIMIT = 1
                    result = Vizier.query_region(
                        SkyCoord(ra=ra_q * u.deg, dec=dec_q * u.deg, frame='icrs'),
                        radius=1.0 * u.arcmin,
                        catalog='II/305/archive')
                    if result:
                        ext_type = 4
                        break
                except Exception:
                    pass

    # RJCE ALLWISE: |b| < 16 or high dust, no GLIMPSE/SAGE
    if ext_type == 0 and (abs(cengb) < 16 or maxebv > 0.2):
        ext_type = 2

    # SFD: default
    if ext_type == 0:
        ext_type = 1

    return ext_type


# ---------------------------------------------------------------------------
# Reddening calculation
# ---------------------------------------------------------------------------

def get_reddening(ref, ext_type):
    """
    Calculate E(J-Ks) reddening for reference catalogue sources.

    Python port of delvered_getreddening.pro.

    Updates the following columns in ``ref`` in-place:
      ebv_sfd, ejk, e_ejk, ext_type.

    Parameters
    ----------
    ref : numpy structured array
        Reference catalogue.  Must have fields ra, dec.
        For RJCE methods also needs: jmag, hmag, kmag, and one of
        w2mag (type 2), gl_45mag (type 3), sage_45mag (type 4).
    ext_type : int
        Extinction method: 1=SFD, 2=RJCE ALLWISE, 3=RJCE GLIMPSE,
        4=RJCE SAGE.
    """
    # --- SFD baseline ---
    ebv = _sfd_ebv(ref['ra'], ref['dec'])
    ref['ebv_sfd'] = ebv
    ejk_sfd = 1.5 * 0.302 * ebv
    ref['ejk'] = ejk_sfd.copy()
    ref['e_ejk'] = 0.1
    high_dust = ebv > 0.3
    ref['e_ejk'][high_dust] = 1.0
    ref['ext_type'] = 1

    if ext_type < 2:
        _flag_bad_reddening(ref, ejk_sfd)
        return

    # --- RJCE GLIMPSE (type 3) ---
    if ext_type == 3 and 'gl_45mag' in ref.dtype.names:
        good = ((ref['jmag'] < 50) & (ref['hmag'] < 50) &
                (ref['kmag'] < 50) & (ref['gl_45mag'] < 50))
        if good.any():
            ejk = 1.5 * 0.918 * (ref['hmag'][good] - ref['gl_45mag'][good] - 0.08)
            e_ejk = 1.5 * 0.918 * np.sqrt(ref['e_hmag'][good]**2 + ref['e_gl_45mag'][good]**2)
            better = ejk < ejk_sfd[good]
            idx = np.where(good)[0][better]
            ref['ejk'][idx] = np.maximum(ejk[better], 0.0)
            ref['e_ejk'][idx] = e_ejk[better]
            ref['ext_type'][idx] = 3

    # --- RJCE SAGE (type 4) ---
    if ext_type == 4 and 'sage_45mag' in ref.dtype.names:
        good = ((ref['jmag'] < 50) & (ref['hmag'] < 50) &
                (ref['kmag'] < 50) & (ref['sage_45mag'] < 50))
        if good.any():
            ejk = 1.5 * 0.918 * (ref['hmag'][good] - ref['sage_45mag'][good] - 0.08)
            e_ejk = 1.5 * 0.918 * np.sqrt(ref['e_hmag'][good]**2 + ref['e_sage_45mag'][good]**2)
            better = ejk < ejk_sfd[good]
            idx = np.where(good)[0][better]
            ref['ejk'][idx] = np.maximum(ejk[better], 0.0)
            ref['e_ejk'][idx] = e_ejk[better]
            ref['ext_type'][idx] = 4

    # --- RJCE ALLWISE (type 2) — fills in where GLIMPSE/SAGE unavailable ---
    if 'w2mag' in ref.dtype.names:
        good = ((ref['jmag'] < 50) & (ref['hmag'] < 50) &
                (ref['kmag'] < 50) & (ref['w2mag'] < 50) &
                (ref['ext_type'] <= 1))
        if good.any():
            ejk = 1.5 * 0.918 * (ref['hmag'][good] - ref['w2mag'][good] - 0.05)
            e_ejk = 1.5 * 0.918 * np.sqrt(ref['e_hmag'][good]**2 + ref['e_w2mag'][good]**2)
            better = ejk < ejk_sfd[good]
            idx = np.where(good)[0][better]
            ref['ejk'][idx] = np.maximum(ejk[better], 0.0)
            ref['e_ejk'][idx] = e_ejk[better]
            ref['ext_type'][idx] = 2

    _flag_bad_reddening(ref, ejk_sfd)

    # Fix NaNs in e_ejk
    bad_err = ~np.isfinite(ref['e_ejk'])
    ref['e_ejk'][bad_err] = 9.99


def _flag_bad_reddening(ref, ejk_sfd):
    """Mark sources where only high-EBV SFD is available as unreliable."""
    bad = (ref['ext_type'] == 1) & (ref['ebv_sfd'] > 0.3)
    if bad.any():
        ref['ejk'][bad] = 999999.0
        ref['e_ejk'][bad] = 999999.0
        ref['ext_type'][bad] = 0


# ---------------------------------------------------------------------------
# Gnomonic projection helper
# ---------------------------------------------------------------------------

def _gnomonic_reverse(xi, eta, ra0, dec0):
    """
    Convert gnomonic (tangent-plane) offsets (degrees) to RA/DEC.

    Parameters
    ----------
    xi, eta : array-like
        Tangent-plane offsets in degrees.
    ra0, dec0 : float
        Tangent point (degrees).

    Returns
    -------
    ra, dec : numpy arrays (degrees)
    """
    xi = np.deg2rad(np.asarray(xi, dtype=float))
    eta = np.deg2rad(np.asarray(eta, dtype=float))
    ra0 = np.deg2rad(ra0)
    dec0 = np.deg2rad(dec0)

    denom = np.cos(dec0) - eta * np.sin(dec0)
    ra = ra0 + np.arctan2(xi, denom)
    dec = np.arctan2(np.sin(dec0) + eta * np.cos(dec0),
                     np.sqrt(xi**2 + denom**2))
    return np.rad2deg(ra) % 360, np.rad2deg(dec)
