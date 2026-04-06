"""
Reference catalogue querying and preparation for the DELVERED pipeline.

Python port of:
  - delvered_getrefcat.pro    → get_ref_cat()
  - delvered_refcat_prep.pro  → refcat_prep()
  - redownload_refcats.pro    → redownload_refcats()
"""

import os
import gzip
import shutil
import subprocess
import tempfile
import time
import numpy as np
from pathlib import Path
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u


# ---------------------------------------------------------------------------
# Catalogue name normalisation
# ---------------------------------------------------------------------------

_REFNAME_MAP = {
    'II/312/AIS': 'GALEX',
    '2MASS-PSC': 'TMASS',
    '2MASS': 'TMASS',
    'GAIA/GAIA': 'GAIA',
    'Skymapper': 'SKYMAPPER',
    'skymapperdr4': 'SKYMAPPERDR4',
    'GLIMPSE': 'II/293/glimpse',
    'SAGE': 'II/305/archive',
    'ATLASREFCAT2': 'ATLAS',
}

# DataLab / psql catalogue configurations
_DATALAB_CATS = {
    'TMASS': {
        'table': 'twomass.psc',
        'cols': ('designation,ra as raj2000,dec as dej2000,'
                 'j_m as jmag,j_cmsig as e_jmag,'
                 'h_m as hmag,h_cmsig as e_hmag,'
                 'k_m as kmag,k_cmsig as e_kmag,ph_qual as qflg'),
        'server': 'db02.datalab.noirlab.edu',
        'user': 'dlquery',
        'racol': 'ra', 'deccol': 'dec',
    },
    'GAIA': {
        'table': 'gaia_dr1.gaia_source',
        'cols': ('source_id as source,ra as ra_icrs,ra_error as e_ra_icrs,'
                 'dec as de_icrs,dec_error as e_de_icrs,'
                 'phot_g_mean_flux as fg,phot_g_mean_flux_error as e_fg,'
                 'phot_g_mean_mag as gmag'),
        'server': 'db02.datalab.noirlab.edu',
        'user': 'dlquery',
        'racol': 'ra', 'deccol': 'dec',
    },
    'GAIADR2': {
        'table': 'gaia_dr2.gaia_source',
        'cols': ('source_id as source,ra,ra_error,dec,dec_error,'
                 'parallax,parallax_error,pmra,pmra_error,pmdec,pmdec_error,'
                 'phot_g_mean_flux as fg,phot_g_mean_flux_error as e_fg,'
                 'phot_g_mean_mag as gmag,'
                 'phot_bp_mean_mag as bp,phot_bp_mean_flux as fbp,phot_bp_mean_flux_error as e_fbp,'
                 'phot_rp_mean_mag as rp,phot_rp_mean_flux as frp,phot_rp_mean_flux_error as e_frp'),
        'server': 'db02.datalab.noirlab.edu',
        'user': 'dlquery',
        'racol': 'ra', 'deccol': 'dec',
    },
    'GAIADR3': {
        'table': 'gaia_dr3.gaia_source',
        'cols': ('source_id as source,ra,ra_error,dec,dec_error,'
                 'parallax,parallax_error,pmra,pmra_error,pmdec,pmdec_error,'
                 'phot_g_mean_flux as fg,phot_g_mean_flux_error as e_fg,'
                 'phot_g_mean_mag as gmag,'
                 'phot_bp_mean_mag as bp,phot_bp_mean_flux as fbp,phot_bp_mean_flux_error as e_fbp,'
                 'phot_rp_mean_mag as rp,phot_rp_mean_flux as frp,phot_rp_mean_flux_error as e_frp'),
        'server': 'db02.datalab.noirlab.edu',
        'user': 'dlquery',
        'racol': 'ra', 'deccol': 'dec',
    },
    'PS': {
        'table': 'public.ps1',
        'cols': 'ra, dec, g as gmag, r as rmag, i as imag, z as zmag, y as ymag',
        'server': 'gp02.datalab.noirlab.edu',
        'user': 'datalab',
        'racol': 'ra', 'deccol': 'dec',
    },
    'SKYMAPPER': {
        'table': 'skymapper_dr1.master',
        'cols': ('raj2000, dej2000, u_psf as sm_umag, e_u_psf as e_sm_umag,'
                 'g_psf as sm_gmag, e_g_psf as e_sm_gmag,'
                 'r_psf as sm_rmag, e_r_psf as e_sm_rmag,'
                 'i_psf as sm_imag, e_i_psf as e_sm_imag,'
                 'z_psf as sm_zmag, e_z_psf as e_sm_zmag'),
        'server': 'db02.datalab.noirlab.edu',
        'user': 'dlquery',
        'racol': 'raj2000', 'deccol': 'dej2000',
    },
    'SKYMAPPERDR4': {
        'table': 'skymapper_dr4.master',
        'cols': ('raj2000, dej2000, u_psf as sm_umag, e_u_psf as e_sm_umag,'
                 'g_psf as sm_gmag, e_g_psf as e_sm_gmag,'
                 'r_psf as sm_rmag, e_r_psf as e_sm_rmag,'
                 'i_psf as sm_imag, e_i_psf as e_sm_imag,'
                 'z_psf as sm_zmag, e_z_psf as e_sm_zmag'),
        'server': 'db02.datalab.noirlab.edu',
        'user': 'dlquery',
        'racol': 'raj2000', 'deccol': 'dej2000',
    },
    'ALLWISE': {
        'table': 'allwise.source',
        'cols': 'ra, dec, w1mpro as w1mag, w1sigmpro as e_w1mag, w2mpro as w2mag, w2sigmpro as e_w2mag',
        'server': 'db02.datalab.noirlab.edu',
        'user': 'dlquery',
        'racol': 'ra', 'deccol': 'dec',
    },
    'ATLAS': {
        'table': 'atlasrefcat2',
        'cols': ('objid,ra,dec,plx as parallax,dplx as parallax_error,'
                 'pmra,dpmra as pmra_error,pmdec,dpmdec as pmdec_error,'
                 'gaia,dgaia as gaiaerr,bp,dbp as bperr,rp,drp as rperr,'
                 'teff,agaia,dupvar,ag,rp1,r1,r10,'
                 'g as gmag,dg as gerr,gchi,gcontrib,'
                 'r as rmag,dr as rerr,rchi,rcontrib,'
                 'i as imag,di as ierr,ichi,icontrib,'
                 'z as zmag,dz as zerr,zchi,zcontrib,nstat,'
                 'j as jmag,dj as jerr,h as hmag,dh as herr,k as kmag,dk as kerr'),
        'server': 'gp10.datalab.noirlab.edu',
        'user': 'datalab',
        'racol': 'ra', 'deccol': 'dec',
    },
}


# ---------------------------------------------------------------------------
# Main reference catalogue query
# ---------------------------------------------------------------------------

def get_ref_cat(cenra, cendec, radius, refcat, file=None,
                saveref=False, silent=False, logfile=None,
                tmpdir='/tmp'):
    """
    Query a reference catalogue for a circular sky region.

    Python port of delvered_getrefcat.pro.

    For major catalogues (Gaia, 2MASS, PS1, SkyMapper, ALLWISE, ATLAS)
    a PostgreSQL Q3C cone search on the NOIRLab DataLab servers is used.
    Other catalogues fall back to astroquery.vizier.

    If ``file`` already exists on disk the cached version is returned.

    Parameters
    ----------
    cenra : float
        Central RA (degrees, J2000).
    cendec : float
        Central DEC (degrees, J2000).
    radius : float
        Search radius (degrees).
    refcat : str
        Catalogue name.  Supported: 'GAIADR3', 'GAIADR2', 'GAIA', 'TMASS'
        (2MASS), 'PS', 'SKYMAPPER', 'SKYMAPPERDR4', 'ALLWISE', 'ATLAS',
        'GALEX', and any VizieR catalogue ID.
    file : str, optional
        Path to save/load a cached FITS version of the result.
    saveref : bool
        If True, save the result to ``file``.
    silent : bool
        Suppress printed output.
    logfile : str, optional
        File path for log messages (currently unused; messages go to stdout).
    tmpdir : str
        Directory for temporary files.

    Returns
    -------
    numpy structured array or None
        Reference catalogue, or None if the query fails.
    """
    t0 = time.time()

    # Normalise catalogue name
    refname = _REFNAME_MAP.get(refcat, refcat.upper())

    # Default cache file name
    if file is None:
        file = os.path.join(tmpdir,
                            f'ref_{cenra:.5f}_{cendec:.5f}_{radius:.3f}_{refname}.fits')

    if not silent:
        print(f"get_ref_cat: querying {refname}: RA={cenra:.5f} DEC={cendec:.5f} "
              f"Radius={radius:.3f}")

    # Return cached result if available
    if os.path.exists(file):
        if not silent:
            print(f"get_ref_cat: loading cached file {file}")
        return fits.getdata(file, 1)

    # --- DataLab PostgreSQL query ---
    if refname in _DATALAB_CATS:
        ref = _query_datalab(cenra, cendec, radius, refname,
                             tmpdir=tmpdir, silent=silent)
        if ref is not None and refname == 'ATLAS':
            ref = _fix_atlas(ref)
    else:
        # --- VizieR fallback ---
        ref = _query_vizier(cenra, cendec, radius, refname,
                            tmpdir=tmpdir, silent=silent)

    if ref is None:
        return None

    if not silent:
        dt = time.time() - t0
        print(f"get_ref_cat: {len(ref)} sources found   dt={dt:.1f} sec.")

    if saveref and file:
        if not silent:
            print(f"get_ref_cat: saving catalogue to {file}")
        hdu = fits.BinTableHDU(ref)
        hdu.writeto(file, overwrite=True)

    return ref


# ---------------------------------------------------------------------------
# DataLab psql query
# ---------------------------------------------------------------------------

def _query_datalab(cenra, cendec, radius, refname, tmpdir='/tmp', silent=False):
    """Run a Q3C cone search on a NOIRLab DataLab Postgres server."""
    cfg = _DATALAB_CATS[refname]
    table = cfg['table']
    cols = cfg['cols']
    server = cfg['server']
    user = cfg['user']
    racol = cfg['racol']
    deccol = cfg['deccol']

    # Check psql is available
    psql = shutil.which('psql')
    if psql is None:
        print("get_ref_cat: psql not found; cannot query DataLab")
        return None

    with tempfile.NamedTemporaryFile(suffix='.txt', dir=tmpdir,
                                     delete=False, mode='w') as tf:
        txtfile = tf.name

    sql = (f"SELECT {cols} FROM {table} "
           f"WHERE q3c_radial_query({racol},{deccol},"
           f"{cenra:.4f},{cendec:.4f},{radius:.3f})")
    cmd = [psql, '-h', server, '-U', user, '-d', 'tapdb',
           '-w', '--pset', 'footer', '-c', sql]
    try:
        with open(txtfile, 'w') as fout:
            subprocess.run(cmd, stdout=fout, stderr=subprocess.DEVNULL,
                           check=True, timeout=600)
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError) as e:
        if not silent:
            print(f"get_ref_cat: psql query failed: {e}")
        os.unlink(txtfile)
        return None

    try:
        ref = _parse_psql_output(txtfile)
    finally:
        if os.path.exists(txtfile):
            os.unlink(txtfile)

    return ref


def _parse_psql_output(txtfile):
    """Parse the pipe-delimited psql output into a numpy structured array."""
    import pandas as pd

    try:
        df = pd.read_csv(txtfile, sep='|', skiprows=1, skipinitialspace=True)
        # Remove leading/trailing whitespace from column names
        df.columns = [c.strip().lower() for c in df.columns]
        # Remove the separator line (all dashes) that psql inserts
        df = df[~df.iloc[:, 0].astype(str).str.fullmatch(r'-+')]
        # Drop last empty row if present
        df = df.dropna(how='all')
        df = df.reset_index(drop=True)
        if len(df) == 0:
            return None
        # Convert to numpy structured array via astropy Table
        from astropy.table import Table
        tbl = Table.from_pandas(df.infer_objects())
        return np.array(tbl)
    except Exception as e:
        print(f"_parse_psql_output: failed to parse {txtfile}: {e}")
        return None


def _fix_atlas(ref):
    """Replace 0.0 magnitudes/errors in ATLAS with 99.99/9.99."""
    mag_cols = ['gaia', 'bp', 'rp', 'gmag', 'rmag', 'imag', 'zmag',
                'jmag', 'hmag', 'kmag']
    err_cols = ['gaiaerr', 'bperr', 'rperr', 'gerr', 'rerr', 'ierr',
                'zerr', 'jerr', 'herr', 'kerr']
    names_lower = [n.lower() for n in ref.dtype.names]
    ref = ref.copy()
    for mc in mag_cols:
        if mc in names_lower:
            idx = names_lower.index(mc)
            col = ref.dtype.names[idx]
            bad = ref[col] <= 0.0
            ref[col][bad] = 99.99
    for ec in err_cols:
        if ec in names_lower:
            idx = names_lower.index(ec)
            col = ref.dtype.names[idx]
            bad = ref[col] <= 0.0
            ref[col][bad] = 9.99
    return ref


# ---------------------------------------------------------------------------
# VizieR fallback
# ---------------------------------------------------------------------------

def _query_vizier(cenra, cendec, radius, refname, tmpdir='/tmp', silent=False):
    """Query VizieR using astroquery."""
    try:
        from astroquery.vizier import Vizier
    except ImportError:
        print("get_ref_cat: astroquery not available; cannot query VizieR")
        return None

    coord = SkyCoord(ra=cenra * u.deg, dec=cendec * u.deg, frame='icrs')
    Vizier.ROW_LIMIT = -1
    Vizier.TIMEOUT = 600

    try:
        result = Vizier.query_region(coord,
                                     width=radius * 60 * u.arcmin,
                                     catalog=refname)
    except Exception as e:
        if not silent:
            print(f"get_ref_cat: VizieR query failed: {e}")
        return None

    if not result or len(result) == 0:
        if not silent:
            print("get_ref_cat: no results from VizieR")
        return None

    tbl = result[0]
    return np.array(tbl)


# ---------------------------------------------------------------------------
# Reference catalogue preparation per chip
# ---------------------------------------------------------------------------

def refcat_prep(exposure, refcatfile, offset=0.2):
    """
    Prepare per-chip reference catalogues from a master field catalogue.

    Python port of delvered_refcat_prep.pro.

    Reads the master reference catalogue, then for each of the 62 DECam
    chips extracts the sources within the chip footprint (plus ``offset``
    margin) in gnomonic coordinates and writes a gzip-compressed
    chip-specific FITS file.

    Parameters
    ----------
    exposure : str
        Exposure name, e.g. ``'F1-00912345'``.
    refcatfile : str
        Path to the master field reference catalogue FITS file.
    offset : float
        Spatial margin in degrees to add around each chip's footprint.
    """
    from astropy.wcs import WCS
    from .extinction import _gnomonic_reverse

    print(f"refcat_prep: exposure = {exposure}")
    print(f"refcat_prep: refcatfile = {refcatfile}")
    print(f"refcat_prep: offset = {offset} deg")

    refcat = fits.getdata(refcatfile, 1)
    ra_ref = refcat['ra'].astype(float)
    dec_ref = refcat['dec'].astype(float)

    # Field centre from catalogue extent
    if np.ptp(ra_ref) > 100:
        ra_adj = np.where(ra_ref > 180, ra_ref - 360, ra_ref)
        cenra = float(np.mean([ra_adj.min(), ra_adj.max()]))
        if cenra < 0:
            cenra += 360
    else:
        cenra = float(np.mean([ra_ref.min(), ra_ref.max()]))
    cendec = float(np.mean([dec_ref.min(), dec_ref.max()]))

    # Gnomonic projection of reference catalogue
    from astropy.coordinates import SkyCoord
    from astropy.wcs import utils as wcsutils

    def gnomic_fwd(ra, dec, cenra, cendec):
        """RA/DEC → gnomonic offsets (degrees)."""
        ra_r = np.deg2rad(ra)
        dec_r = np.deg2rad(dec)
        ra0 = np.deg2rad(cenra)
        dec0 = np.deg2rad(cendec)
        denom = (np.sin(dec0) * np.sin(dec_r) +
                 np.cos(dec0) * np.cos(dec_r) * np.cos(ra_r - ra0))
        xi = np.cos(dec_r) * np.sin(ra_r - ra0) / denom
        eta = (np.cos(dec0) * np.sin(dec_r) -
               np.sin(dec0) * np.cos(dec_r) * np.cos(ra_r - ra0)) / denom
        return np.rad2deg(xi), np.rad2deg(eta)

    reflon, reflat = gnomic_fwd(ra_ref, dec_ref, cenra, cendec)

    field = exposure.split('-')[0]

    for c in range(1, 63):
        schip = f'{c:02d}'
        chipfile = f'{field}/chip{schip}/{exposure}_{schip}.fits'
        if not os.path.exists(chipfile):
            print(f"  {chipfile} NOT FOUND")
            continue

        with fits.open(chipfile) as hdul:
            hdr = hdul[0].header
            nx = hdr.get('NAXIS1', 2048)
            ny = hdr.get('NAXIS2', 4096)
            wcs = WCS(hdr)

        # Chip corner coordinates
        px = np.array([0, nx - 1, nx - 1, 0], dtype=float)
        py = np.array([0, 0, ny - 1, ny - 1], dtype=float)
        world = wcs.pixel_to_world(px, py)
        vra = np.array([w.ra.deg for w in world])
        vdec = np.array([w.dec.deg for w in world])
        vlon, vlat = gnomic_fwd(vra, vdec, cenra, cendec)

        # Select reference stars within chip footprint + margin
        inside = ((reflon >= vlon.min() - offset) &
                  (reflon <= vlon.max() + offset) &
                  (reflat >= vlat.min() - offset) &
                  (reflat <= vlat.max() + offset))
        if not inside.any():
            continue

        refcat1 = refcat[inside]
        chipreffile = f'{field}/chip{schip}/{exposure}_{schip}_refcat.fits'
        hdu = fits.BinTableHDU(refcat1)
        hdu.writeto(chipreffile, overwrite=True)

        # Gzip compress
        gz_file = chipreffile + '.gz'
        with open(chipreffile, 'rb') as f_in, gzip.open(gz_file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        os.unlink(chipreffile)
        print(f"  wrote {gz_file}  ({inside.sum()} sources)")


# ---------------------------------------------------------------------------
# Redownload reference catalogues
# ---------------------------------------------------------------------------

def redownload_refcats(nights_dir, refcat_name='GAIADR3', **kwargs):
    """
    Rebuild per-chip reference catalogues from freshly-downloaded masters.

    Python port of redownload_refcats.pro.

    Parameters
    ----------
    nights_dir : str or Path
        Root directory containing night sub-directories.
    refcat_name : str
        Reference catalogue name passed to ``get_ref_cat``.
    **kwargs
        Additional keyword arguments passed to ``get_ref_cat``.
    """
    import glob

    nights_dir = Path(nights_dir)
    night_dirs = sorted([p for p in nights_dir.iterdir()
                         if p.is_dir() and p.name.isdigit()])

    for night_dir in night_dirs:
        print(f"redownload_refcats: processing night {night_dir.name}")
        refcat_dir = night_dir / 'refcat'
        refcat_dir.mkdir(exist_ok=True)

        # Find field directories
        field_dirs = sorted([p for p in night_dir.iterdir() if p.is_dir()])
        for field_dir in field_dirs:
            fitsfiles = list(field_dir.glob('*.fits')) + list(field_dir.glob('*.fits.fz'))
            if not fitsfiles:
                continue
            # Get field name from directory
            field = field_dir.name
            refcat_file = refcat_dir / f'{field}_refcat.fits.gz'
            if refcat_file.exists():
                print(f"  {field}: cached refcat found")
                continue
            # Need RA/DEC centre: read first FITS header
            try:
                hdr = fits.getheader(str(fitsfiles[0]))
                cenra = float(hdr.get('RA', hdr.get('CRVAL1', 0)))
                cendec = float(hdr.get('DEC', hdr.get('CRVAL2', 0)))
            except Exception:
                continue

            ref = get_ref_cat(cenra, cendec, radius=1.1,
                              refcat=refcat_name, silent=False, **kwargs)
            if ref is None:
                continue
            # Write gzip-compressed FITS
            tmp_fits = str(refcat_file).replace('.gz', '')
            fits.writeto(tmp_fits, ref, overwrite=True)
            with open(tmp_fits, 'rb') as f_in, gzip.open(str(refcat_file), 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
            os.unlink(tmp_fits)
            print(f"  {field}: wrote {refcat_file}")
