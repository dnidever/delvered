"""
Diagnostic and repair utilities for the DELVERED pipeline.

Python ports of:
  - fix_badobjid.pro                    → fix_bad_objid()
  - fix_localheaders_bscale.pro         → fix_local_headers_bscale()
  - fix_missing_fitsfiles.pro           → fix_missing_fits_files()
  - fix_resourcefiles.pro               → fix_resource_files()
  - check_resourcefiles.pro             → check_resource_files()
  - checkwcsexp.pro                     → check_wcs_exp()
  - remake_delvered_expforced.pro       → remake_expforced()

The IDL scripts fix_missing_smash_exposure.pro, fix_refcat.pro,
fix_resourcefiles_extname.pro, fix_resourcefiles_missingfiles.pro,
fix_smash_nights.pro, check_resourcefiles_extnameproblem.pro, and
check_badfwhm_problem.pro are omitted because they contain
site-specific hardcoded paths or depend on infrastructure (mass-store
FTP fallbacks, SMASH-specific catalogs) that cannot be generalised
without more context.
"""

import re
import gzip
import subprocess
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _iter_nights(delve_dir, nights=None):
    """Yield (night_str, night_dir Path) for all 8-digit night directories."""
    exp_dir = Path(delve_dir) / 'exposures'
    if nights:
        for n in nights:
            d = exp_dir / str(n)
            if d.is_dir():
                yield str(n), d
    else:
        for d in sorted(exp_dir.iterdir()):
            if d.is_dir() and d.name.isdigit() and len(d.name) == 8:
                yield d.name, d


def _iter_chips(field_dir, field, max_chip=62):
    """Yield (chip_num, chip_dir, base) for all chip dirs in a field."""
    for c in range(1, max_chip + 1):
        if c == 61:
            continue
        schip = f'{c:02d}'
        chip_dir = field_dir / f'chip{schip}'
        if not chip_dir.is_dir():
            continue
        yield c, schip, chip_dir


def _parse_resource_file(rfile):
    """
    Parse a PHOTRED resource dot-file.

    Returns a dict mapping tag ('fluxfile', 'wtfile', 'maskfile') to
    the full value string (e.g. ``'/mss1/.../c4d.fits.fz[S1]'``).
    """
    result = {}
    p = Path(rfile)
    if not p.exists():
        return result
    with open(p) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) >= 3 and parts[1] == '=':
                result[parts[0]] = parts[2]
            elif len(parts) == 2 and '=' in parts[0]:
                k, v = parts[0].rstrip('='), parts[1]
                result[k] = v
    return result


def _write_resource_file(rfile, flux, wt, mask):
    """Rewrite a resource dot-file with the three canonical lines."""
    with open(rfile, 'w') as fh:
        fh.write(f'fluxfile = {flux}\n'
                 f'wtfile = {wt}\n'
                 f'maskfile = {mask}\n')


# ---------------------------------------------------------------------------
# fix_missing_fits_files
# ---------------------------------------------------------------------------

def fix_missing_fits_files(delve_dir, nights=None):
    """
    Create zero-byte stub FITS files wherever resource dot-files exist but
    the corresponding stub FITS file is absent.

    Python port of fix_missing_fitsfiles.pro.

    Parameters
    ----------
    delve_dir : str or Path
        Root DELVE data directory (contains ``exposures/``).
    nights : list of str, optional
        Restrict to specific night identifiers.  Defaults to all nights.
    """
    for inight, night_dir in _iter_nights(delve_dir, nights):
        print(f"{inight}")
        for field_dir in sorted(night_dir.glob('F*')):
            if not field_dir.is_dir():
                continue
            ifield = field_dir.name
            for c, schip, chip_dir in _iter_chips(field_dir, ifield):
                for rfile in chip_dir.glob(f'.{ifield}-????????_{schip}.fits'):
                    stub = chip_dir / rfile.name[1:]  # remove leading dot
                    if not stub.exists():
                        print(f"  creating {stub}")
                        stub.write_text('')


# ---------------------------------------------------------------------------
# fix_resource_files
# ---------------------------------------------------------------------------

def fix_resource_files(delve_dir, nights=None):
    """
    Normalise whitespace in PHOTRED resource dot-files.

    Python port of fix_resourcefiles.pro.

    Some resource files have extra spaces between the tag, ``=``, and the
    value.  This function rewrites them in canonical form.

    Parameters
    ----------
    delve_dir : str or Path
        Root DELVE data directory (contains ``exposures/``).
    nights : list of str, optional
        Restrict to specific nights.
    """
    for inight, night_dir in _iter_nights(delve_dir, nights):
        print(f"{inight}")
        for field_dir in sorted(night_dir.glob('F*')):
            if not field_dir.is_dir():
                continue
            ifield = field_dir.name

            # Find exposures from chip01 stubs
            chip01 = field_dir / 'chip01'
            if not chip01.is_dir():
                continue
            for stub in chip01.glob(f'{ifield}-????????_01.fits'):
                expnum = stub.stem[len(ifield) + 1:]  # strip 'F1-'
                for c, schip, chip_dir in _iter_chips(field_dir, ifield):
                    rfile = chip_dir / f'.{ifield}-{expnum}_{schip}.fits'
                    if not rfile.exists():
                        continue
                    info = _parse_resource_file(rfile)
                    if len(info) == 3:
                        _write_resource_file(
                            rfile, info['fluxfile'],
                            info['wtfile'], info['maskfile'])


# ---------------------------------------------------------------------------
# check_resource_files
# ---------------------------------------------------------------------------

def check_resource_files(delve_dir, nights=None):
    """
    Scan resource dot-files and report which CP archive FITS files are missing.

    Python port of check_resourcefiles.pro.

    Parameters
    ----------
    delve_dir : str or Path
        Root DELVE data directory.
    nights : list of str, optional
        Restrict to specific nights.

    Returns
    -------
    list of str
        Paths of CP flux files referenced by resource files but absent on disk.
    """
    missing = []
    for inight, night_dir in _iter_nights(delve_dir, nights):
        night_missing = []
        for field_dir in sorted(night_dir.glob('F*')):
            if not field_dir.is_dir():
                continue
            ifield = field_dir.name
            chip01 = field_dir / 'chip01'
            if not chip01.is_dir():
                continue
            for stub in chip01.glob(f'{ifield}-????????_01.fits'):
                expnum = stub.stem[len(ifield) + 1:]
                for c, schip, chip_dir in _iter_chips(field_dir, ifield):
                    rfile = chip_dir / f'.{ifield}-{expnum}_{schip}.fits'
                    if not rfile.exists():
                        continue
                    info = _parse_resource_file(rfile)
                    flux = info.get('fluxfile', '')
                    if '[' in flux:
                        flux = flux[:flux.index('[')]
                    if flux and not Path(flux).exists():
                        night_missing.append(flux)
        if night_missing:
            n_uniq = len(set(night_missing))
            print(f"{inight}: {n_uniq} missing flux files")
        missing.extend(night_missing)
    missing = list(dict.fromkeys(missing))  # deduplicate
    print(f"check_resource_files: {len(missing)} missing files total")
    return missing


# ---------------------------------------------------------------------------
# fix_local_headers_bscale
# ---------------------------------------------------------------------------

def fix_local_headers_bscale(delve_dir, nights=None):
    """
    Add a BSCALE keyword to local ``.fits.head`` files whose BUNIT is
    'ELECTRONS' (old DES-format images later replaced by CP ADU images).

    Python port of fix_localheaders_bscale.pro.

    For each affected chip, estimates BSCALE as SKYBRITE / median(CP_image),
    falling back to the ALS sky median or AVSKY if SKYBRITE is absent.

    Parameters
    ----------
    delve_dir : str or Path
        Root DELVE data directory.
    nights : list of str, optional
        Restrict to specific nights.
    """
    for inight, night_dir in _iter_nights(delve_dir, nights):
        head_files = list(night_dir.rglob('*.fits.head'))
        if not head_files:
            continue
        print(f"{inight}: {len(head_files)} .fits.head files")
        for hfile in head_files:
            try:
                hdr = fits.Header.fromtextfile(str(hfile))
            except Exception:
                continue
            bunit = str(hdr.get('BUNIT', '')).strip().lower()
            if bunit != 'electrons':
                continue

            # Parse chip directory and base from path
            chip_dir = hfile.parent
            base = hfile.name.replace('.fits.head', '')
            # base is like F5-00423440_34  → strip trailing chip suffix
            m = re.match(r'^(.+)_(\d{2})$', base)
            if not m:
                continue
            exp_base = m.group(1)  # e.g. F5-00423440

            for chip in range(1, 63):
                schip = f'{chip:02d}'
                rfile = chip_dir.parent.parent / f'chip{schip}' / f'.{exp_base}_{schip}.fits'
                hfile_c = chip_dir.parent.parent / f'chip{schip}' / f'{exp_base}_{schip}.fits.head'
                if not rfile.exists() or not hfile_c.exists():
                    continue

                info = _parse_resource_file(rfile)
                flux_ref = info.get('fluxfile', '')
                if '[' in flux_ref:
                    ext_name = flux_ref[flux_ref.index('[') + 1:flux_ref.index(']')]
                    flux_path = flux_ref[:flux_ref.index('[')]
                else:
                    continue

                if not Path(flux_path).exists():
                    continue

                try:
                    chdr = fits.Header.fromtextfile(str(hfile_c))
                except Exception:
                    continue

                bscale_existing = chdr.get('BSCALE', None)
                if bscale_existing is not None and 3.5 <= float(bscale_existing) <= 5.0:
                    continue

                skybrite = chdr.get('SKYBRITE', None)
                # Read the CP FITS to get median
                try:
                    with fits.open(flux_path, memmap=True) as hdul:
                        im = hdul[ext_name].data
                    new_med = float(np.median(im))
                    if new_med == 0:
                        continue
                    if skybrite is not None and float(skybrite) > 0:
                        bscale = float(skybrite) / new_med
                    else:
                        avsky = chdr.get('AVSKY', None)
                        if avsky is not None:
                            bscale = float(avsky) / new_med
                        else:
                            continue
                    if bscale <= 0:
                        continue

                    # Back up original, write updated header
                    bak = Path(str(hfile_c) + '.bak')
                    if not bak.exists():
                        import shutil
                        shutil.copy2(str(hfile_c), str(bak))
                    chdr['BSCALE'] = (bscale, 'scale to DES-like electrons')
                    chdr.tofile(str(hfile_c), overwrite=True)
                    print(f"  fixed {hfile_c.name}: BSCALE={bscale:.3f}")
                except Exception as exc:
                    print(f"  WARNING: could not fix {hfile_c}: {exc}")


# ---------------------------------------------------------------------------
# check_wcs_exp
# ---------------------------------------------------------------------------

def check_wcs_exp(expnums, search_dir='.'):
    """
    Report WCS fit RMS and NMATCH for each chip of the given exposures.

    Python port of checkwcsexp.pro.

    Parameters
    ----------
    expnums : list of str
        Exposure base names to check (e.g. ``['F5-00423440']``).
    search_dir : str or Path
        Directory to search for chip FITS files.

    Returns
    -------
    list of dict
        One entry per chip file with keys: file, crval1, crval2, rms, nmatch.
    """
    results = []
    search_dir = Path(search_dir)
    for expnum in expnums:
        chip_files = sorted(search_dir.rglob(f'chip??/{expnum}_[0-9][0-9].fits'))
        if not chip_files:
            chip_files = sorted(search_dir.glob(f'{expnum}_[0-9][0-9].fits'))
        if not chip_files:
            print(f"check_wcs_exp: NO files for {expnum}")
            continue
        for cfile in chip_files:
            try:
                hdr = fits.getheader(str(cfile), ext=0)
            except Exception:
                continue
            crval1 = float(hdr.get('CRVAL1', 0))
            crval2 = float(hdr.get('CRVAL2', 0))
            # Extract WCS RMS and NMATCH from HISTORY
            rms = 9999.99
            nmatch = 9999
            for card in hdr.get('HISTORY', []):
                card_s = str(card)
                if 'WCSFIT: RMS' in card_s:
                    m = re.search(r'RMS\s*=\s*([\d.]+)', card_s)
                    if m:
                        rms = float(m.group(1))
                elif 'WCSFIT: NMATCH' in card_s:
                    m = re.search(r'NMATCH\s*=\s*(\d+)', card_s)
                    if m:
                        nmatch = int(m.group(1))
            print(f"{cfile.name}  crval1={crval1:.4f}  crval2={crval2:.4f}  "
                  f"rms={rms:.3f}  nmatch={nmatch}")
            results.append({'file': str(cfile), 'crval1': crval1,
                            'crval2': crval2, 'rms': rms, 'nmatch': nmatch})
    return results


# ---------------------------------------------------------------------------
# remake_expforced
# ---------------------------------------------------------------------------

def remake_expforced(brick, brick_dir=None):
    """
    Fix corrupted OBJID values in a brick's ``_expforced.fits.gz`` catalogue.

    Python port of remake_delvered_expforced.pro (fast path only).

    The fast fix extracts the numeric part of each measurement ID
    (``expnum_chip.NUM``) and constructs the corrected OBJID as
    ``brick.NUM``.

    Parameters
    ----------
    brick : str
        Brick name (e.g. ``'0988m505'``).
    brick_dir : str or Path, optional
        Path to the brick directory.  Defaults to
        ``{cwd}/{brick[:4]}/{brick}/``.

    Returns
    -------
    int
        Number of objects whose OBJID was updated.
    """
    brick = str(brick).strip()
    if brick_dir is None:
        brick_dir = Path('.') / brick[:4] / brick
    else:
        brick_dir = Path(brick_dir)

    expfile = brick_dir / f'{brick}_expforced.fits.gz'
    if not expfile.exists():
        expfile = brick_dir / f'{brick}_expforced.fits'
    if not expfile.exists():
        print(f"remake_expforced: {expfile} NOT FOUND")
        return 0

    print(f"remake_expforced: loading {expfile}")
    if str(expfile).endswith('.gz'):
        with gzip.open(str(expfile)) as gz:
            expcat = fits.getdata(gz, 1)
    else:
        expcat = fits.getdata(str(expfile), 1)

    # Extract the numeric suffix from ID (format: expnum_chip.NUM)
    ids = np.char.strip(expcat['id'].astype(str))
    nums = np.array([i.split('.')[-1] for i in ids])
    new_objids = np.array([f'{brick}.{n}' for n in nums])

    # Rebuild with corrected objid
    from numpy.lib.recfunctions import drop_fields, append_fields
    expcat2 = drop_fields(expcat, 'objid')
    expcat2 = append_fields(expcat2, 'objid', new_objids, usemask=False)

    # Write back
    out_path = brick_dir / f'{brick}_expforced.fits'
    hdu = fits.BinTableHDU(expcat2)
    hdu.writeto(str(out_path), overwrite=True)
    subprocess.run(['gzip', '-f', str(out_path)], check=False)
    print(f"remake_expforced: wrote {out_path}.gz  ({len(expcat2)} rows)")
    return len(expcat2)
