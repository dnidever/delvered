"""
DELVERED pre-processing: set up nightly PHOTRED directories.

Python port of:
  - delvered_prep.pro       → prep()
  - delvered_prep_night.pro → prep_night()
"""

import os
import re
import subprocess
import gzip
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u


# CTIO timezone offset (hours behind UTC)
_CTIO_TZ_HOURS = 4.0

# Accepted filter prefixes
_GOOD_FILTERS = ('u', 'g', 'r', 'i', 'z', 'Y')

# PHOTRED setup template
_PHOTRED_SETUP = """\
##### REQUIRED #####
scriptsdir  {scriptsdir}
irafdir     {irafdir}
telescope   Blanco
instrument  DECAM
observatory CTIO
nmulti      10
nmulti_wcs       20
nmulti_daophot   20
nmulti_allframe  10
filtref     g,i,r,z,u,Y
modeleqnfile {modeleqnfile}
trans       delve.trans
##### OPTIONAL #####
sepfielddir  1
sepchipdir   1
keepmef      0
catformat    FITS
workdir      {workdir}
clean        1
skipcheck    1
redo         0
wcsrefname   GAIADR2
searchdist   20
hyperthread   1
daopsfva      1
daofitradfwhm 1.0
psfstars      1
psfcomsrc     0
psfcomglobal  0
psfcomgauss   0
finditer      1
alfdetprog  sextractor
alftrimcomb   0
cmbforce      1
keepinstr     1
avgmag        1
avgonlymag    0
todered       u,g,r,i,z,g-i
sumquick      1
##### STAGES #####
#rename
#split
 wcs
 daophot
 match
#allframe
 apcor
 astrom
 zeropoint
 calib
 combine
 deredden
 save
#html
"""


def _date_to_night(date_obs, tz_hours=_CTIO_TZ_HOURS):
    """
    Convert an ISO date-obs string to an 8-digit YYYYMMDD night identifier.

    CTIO is UTC-4; exposures taken before local midnight belong to the
    previous calendar date.
    """
    from astropy.time import Time
    import datetime

    date_obs = str(date_obs).strip().replace(' ', 'T')
    try:
        t = Time(date_obs, format='isot', scale='utc')
        local_dt = (t - tz_hours * u.hour).to_datetime()
        return local_dt.strftime('%Y%m%d')
    except Exception:
        return ''


def _load_fields_file(path):
    """Return list of (shname, longname) from an ASCII fields file."""
    entries = []
    p = Path(path)
    if not p.exists():
        return entries
    with open(p) as fh:
        for line in fh:
            parts = line.strip().split()
            if len(parts) >= 2 and not parts[0].startswith('#'):
                entries.append((parts[0], parts[1]))
    return entries


def _read_log(path):
    """Return non-empty lines from a log file, or [] if absent."""
    p = Path(path)
    if not p.exists():
        return []
    with open(p) as fh:
        return [ln.strip() for ln in fh if ln.strip()]


def _gnomonic_proj(ra, dec, cen_ra, cen_dec):
    """
    Gnomonic (tangent-plane) projection of (ra, dec) around (cen_ra, cen_dec).
    Returns (xi, eta) in degrees.
    """
    ra = np.radians(ra)
    dec = np.radians(dec)
    cra = np.radians(cen_ra)
    cdec = np.radians(cen_dec)
    denom = (np.sin(cdec) * np.sin(dec) +
             np.cos(cdec) * np.cos(dec) * np.cos(ra - cra))
    xi = np.degrees(np.cos(dec) * np.sin(ra - cra) / denom)
    eta = np.degrees((np.cos(cdec) * np.sin(dec) -
                      np.sin(cdec) * np.cos(dec) * np.cos(ra - cra)) / denom)
    return xi, eta


def _group_exposures_by_field(expstr, radius_deg=0.1):
    """
    Group exposures into fields by positional proximity.

    Uses a simple friends-of-friends-like approach: each un-assigned
    exposure seeds a new field; all exposures within *radius_deg* of
    any member are added.

    Returns
    -------
    numpy array of str
        Field label (e.g. 'F1', 'F2', ...) for each exposure.
    """
    n = len(expstr)
    labels = np.full(n, '', dtype='U20')
    field_count = 1

    coords = SkyCoord(ra=expstr['ra'].astype(float) * u.deg,
                      dec=expstr['dec'].astype(float) * u.deg, frame='icrs')

    for i in range(n):
        if labels[i] != '':
            continue
        fname = f'F{field_count}'
        field_count += 1
        labels[i] = fname
        # Grow field
        changed = True
        while changed:
            changed = False
            in_field = np.where(labels == fname)[0]
            for j in range(n):
                if labels[j] != '':
                    continue
                sep = coords[j].separation(coords[in_field]).deg
                if np.any(sep <= radius_deg):
                    labels[j] = fname
                    changed = True
    return labels


def prep_night(expfile, delve_dir, delvereddir, decam_txt,
               scriptsdir='', irafdir='', workdir='/tmp',
               redo=False):
    """
    Prepare one night's DECam exposures for PHOTRED processing.

    Python port of delvered_prep_night.pro.

    For each exposure listed in *expfile*:

    1. Groups exposures spatially into fields.
    2. Downloads (or reuses) a reference catalogue for each field.
    3. Creates per-chip stub FITS resource files and dot-files with
       CP archive paths for each chip/extension.
    4. Writes per-chip reference catalogue subsets.
    5. Writes ``logs/WCS.inlist``, ``photred.setup``, and ``fields``.

    Parameters
    ----------
    expfile : str or Path
        Per-night FITS file listing exposures (produced by :func:`prep`).
    delve_dir : str or Path
        Root DELVE data directory (contains ``exposures/``).
    delvereddir : str or Path
        Root of the delvered project tree.
    decam_txt : str or Path
        ``data/decam.txt`` mapping extension names to ccdnums.
    scriptsdir, irafdir, workdir : str
        Written into ``photred.setup``.
    redo : bool
        Re-create files even if they already exist.
    """
    from .refcat import get_ref_cat

    expfile = Path(expfile)
    delve_dir = Path(delve_dir)
    delvereddir = Path(delvereddir)
    modeleqnfile = str(delvereddir / 'params' / 'modelmag_equations.txt')

    if not expfile.exists():
        print(f"prep_night: {expfile} NOT FOUND")
        return

    expstr = fits.getdata(str(expfile), 1)
    nexp = len(expstr)
    inight = expfile.parent.name
    print(f"prep_night: {nexp} exposures for night {inight}")

    # Load DECam extension→ccdnum mapping
    decam = np.genfromtxt(str(decam_txt), names=True, dtype=None, encoding='ascii')
    ext_to_ccd = {str(r['name']).strip(): int(r['ccdnum']) for r in decam}

    night_dir = delve_dir / 'exposures' / inight
    setup_file = night_dir / 'photred.setup'
    night_dir.mkdir(parents=True, exist_ok=True)
    (night_dir / 'logs').mkdir(exist_ok=True)

    setup_text = _PHOTRED_SETUP.format(
        scriptsdir=scriptsdir, irafdir=irafdir,
        modeleqnfile=modeleqnfile, workdir=workdir)
    if not setup_file.exists() or redo:
        setup_file.write_text(setup_text)

    # Check which exposures have already been set up
    existing_expnums = set()
    for logname in ('logs/WCS.inlist', 'logs/WCS.success'):
        for line in _read_log(night_dir / logname):
            # Extract exposure number from path like F1/chip34/F1-00806749_34.fits
            m = re.search(r'-(\d{8})_', line)
            if m:
                existing_expnums.add(m.group(1))

    to_add_mask = np.ones(nexp, dtype=bool)
    if existing_expnums and not redo:
        expnums = np.char.strip(expstr['expnum'].astype(str))
        to_add_mask = ~np.isin(expnums, list(existing_expnums))

    exptoadd = expstr[to_add_mask]
    ntoadd = len(exptoadd)
    print(f"prep_night: {ntoadd} exposures to add")
    if ntoadd == 0:
        return

    # Load existing fields file
    old_fields = _load_fields_file(night_dir / 'fields')
    old_fieldnums = set()
    for sn, _ in old_fields:
        m = re.match(r'F(\d+)$', sn)
        if m:
            old_fieldnums.add(int(m.group(1)))
    field_counter_start = (max(old_fieldnums) + 1) if old_fieldnums else 1

    # Group new exposures into fields
    if ntoadd > 1:
        field_labels = _group_exposures_by_field(exptoadd)
        # Re-number starting from field_counter_start
        seen = {}
        counter = field_counter_start
        new_labels = np.empty(ntoadd, dtype='U20')
        for i, lbl in enumerate(field_labels):
            if lbl not in seen:
                seen[lbl] = f'F{counter}'
                counter += 1
            new_labels[i] = seen[lbl]
        field_labels = new_labels
    else:
        field_labels = np.array([f'F{field_counter_start}'])

    # Build field index
    unique_fields, inv = np.unique(field_labels, return_inverse=True)
    new_field_entries = []  # (shname, longname)

    refcat_dir = night_dir / 'refcat'
    refcat_dir.mkdir(exist_ok=True)

    outfiles = []

    for fi, fname in enumerate(unique_fields):
        find = np.where(field_labels == fname)[0]
        nfind = len(find)
        fexp = exptoadd[find]
        field_dir = night_dir / fname
        field_dir.mkdir(parents=True, exist_ok=True)

        cen_ra = float(np.mean(fexp['ra'].astype(float)))
        cen_dec = float(np.mean(fexp['dec'].astype(float)))

        # Get or reuse reference catalogue
        refcat_file = refcat_dir / f'{fname}_refcat.fits.gz'
        refcat = None
        if refcat_file.exists() and not redo:
            try:
                with gzip.open(str(refcat_file)) as gz:
                    refcat = fits.getdata(gz, 1)
            except Exception:
                refcat = None

        if refcat is None:
            print(f"  {fname}: fetching reference catalogue")
            raw_file = str(refcat_dir / f'{fname}_refcat.fits')
            try:
                refcat = get_ref_cat(cen_ra, cen_dec, 1.5,
                                     refcat='GAIADR2', outfile=raw_file)
                subprocess.run(['gzip', '-f', raw_file], check=False)
            except Exception as exc:
                print(f"  WARNING: could not get refcat for {fname}: {exc}")
                refcat = None

        # Gnomonic projection for reference catalogue
        if refcat is not None and len(refcat) > 0:
            ref_ra = refcat['ra'].astype(float)
            ref_dec = refcat['dec'].astype(float)
            ref_xi, ref_eta = _gnomonic_proj(ref_ra, ref_dec, cen_ra, cen_dec)
        else:
            ref_xi = ref_eta = None

        # Get field long name from first exposure FITS header OBJECT keyword
        long_name = fname
        first_flux = str(fexp['fluxfile'][0]).strip().replace('/net/mss1/', '/mss1/')
        if Path(first_flux).exists():
            try:
                hdr0 = fits.getheader(first_flux, ext=0)
                obj = str(hdr0.get('OBJECT', '')).strip()
                if obj:
                    long_name = re.sub(r'\s+', '', obj)
            except Exception:
                pass
        new_field_entries.append((fname, long_name))

        print(f"  {fname}/{long_name}: {nfind} exposure(s)")

        # Exposure loop
        for ei in range(nfind):
            exp = fexp[ei]
            flux_file = str(exp['fluxfile']).strip().replace('/net/mss1/', '/mss1/')
            mask_file = str(exp['maskfile']).strip().replace('/net/mss1/', '/mss1/')
            wt_file = str(exp['wtfile']).strip().replace('/net/mss1/', '/mss1/')
            expnum = str(exp['expnum']).strip()

            if not Path(flux_file).exists():
                print(f"    {expnum}: {flux_file} NOT FOUND, skipping")
                continue

            # Open FITS to read extension names
            try:
                with fits.open(flux_file, memmap=True) as hdul:
                    ext_names = [h.name for h in hdul
                                 if h.name and h.name != 'PRIMARY']
            except Exception as exc:
                print(f"    {expnum}: cannot open {flux_file}: {exc}")
                continue

            # Chip loop
            for ext_name in ext_names:
                ccd = ext_to_ccd.get(ext_name)
                if ccd is None:
                    continue
                schip = f'{ccd:02d}'
                chip_dir = field_dir / f'chip{schip}'
                chip_dir.mkdir(exist_ok=True)

                base = f'{fname}-{expnum}_{schip}'
                stub = chip_dir / f'{base}.fits'
                dot = chip_dir / f'.{base}.fits'
                stub.write_text('')
                dot.write_text(
                    f'fluxfile = {flux_file}[{ext_name}]\n'
                    f'wtfile = {wt_file}[{ext_name}]\n'
                    f'maskfile = {mask_file}[{ext_name}]\n'
                )

                # Per-chip reference catalogue
                chip_refcat_file = chip_dir / f'{base}_refcat.fits.gz'
                if (not chip_refcat_file.exists() or redo) and refcat is not None:
                    try:
                        with fits.open(flux_file, memmap=True) as hdul:
                            chdr = hdul[ext_name].header.copy()
                        if 'ZNAXIS1' in chdr:
                            chdr['NAXIS1'] = chdr['ZNAXIS1']
                            chdr['NAXIS2'] = chdr['ZNAXIS2']
                        wcs = WCS(chdr)
                        nx = int(chdr.get('NAXIS1', 2048))
                        ny = int(chdr.get('NAXIS2', 4096))
                        corners = np.array([[0, 0], [nx-1, 0],
                                            [nx-1, ny-1], [0, ny-1]], float)
                        sky = wcs.pixel_to_world(corners[:, 0], corners[:, 1])
                        vra = sky.ra.deg
                        vdec = sky.dec.deg
                        vxi, veta = _gnomonic_proj(vra, vdec, cen_ra, cen_dec)
                        offset = 0.02
                        if abs(cen_dec) > 70:
                            offset = 0.2
                        chip_mask = (
                            (ref_xi >= vxi.min() - offset) &
                            (ref_xi <= vxi.max() + offset) &
                            (ref_eta >= veta.min() - offset) &
                            (ref_eta <= veta.max() + offset)
                        )
                        chip_ref = refcat[chip_mask]
                        raw = str(chip_refcat_file).replace('.gz', '')
                        fits.BinTableHDU(chip_ref).writeto(raw, overwrite=True)
                        subprocess.run(['gzip', '-f', raw], check=False)
                    except Exception as exc:
                        print(f"    WARNING: refcat trim failed {base}: {exc}")

                outfiles.append(
                    f'{fname}/chip{schip}/{base}.fits')

    # Write WCS.inlist
    wcs_inlist = night_dir / 'logs' / 'WCS.inlist'
    with open(wcs_inlist, 'w') as fh:
        fh.write('\n'.join(outfiles) + '\n')
    print(f"prep_night: wrote {len(outfiles)} entries to {wcs_inlist}")

    # Merge with old fields and ensure unique long names
    all_fields = old_fields + new_field_entries
    seen_names = {}
    for i, (sn, ln) in enumerate(all_fields):
        if ln in seen_names:
            seen_names[ln] += 1
            all_fields[i] = (sn, f'{ln}_{seen_names[ln]}')
        else:
            seen_names[ln] = 1

    fields_file = night_dir / 'fields'
    if fields_file.exists():
        import shutil
        shutil.copy2(str(fields_file), str(fields_file.with_suffix('.bak')))
    with open(fields_file, 'w') as fh:
        for sn, ln in all_fields:
            fh.write(f'{sn}   {ln}\n')


def prep(expfile, delve_dir, delvereddir, decam_txt,
         scriptsdir='', irafdir='', workdir='/tmp',
         nightmin=None, newonly=False, redo=False, nmulti=5,
         smash_dir=None):
    """
    Prepare all DELVE MC nightly directories for PHOTRED processing.

    Python port of delvered_prep.pro.

    Loads the master DECam MC exposure list, applies filter / exposure-time
    cuts, groups exposures by observing night, optionally skips nights that
    have already been processed, and calls :func:`prep_night` for each new
    night (in parallel when *nmulti* > 1).

    Parameters
    ----------
    expfile : str or Path
        The master DECam MC exposure FITS file
        (output of :func:`make_exposures_list`).
    delve_dir : str or Path
        Root DELVE data directory.
    delvereddir : str or Path
        Root of the delvered project tree.
    decam_txt : str or Path
        ``data/decam.txt`` mapping extension names to ccdnums.
    scriptsdir, irafdir, workdir : str
        Written into each night's ``photred.setup``.
    nightmin : str or int, optional
        Minimum night to process (inclusive, YYYYMMDD).
    newonly : bool
        Only process nights not already present in the output tree.
    redo : bool
        Re-process nights even if already set up.
    nmulti : int
        Number of parallel workers.
    smash_dir : str or Path, optional
        Root of the SMASH processed data.  Nights also present in
        *smash_dir* are skipped (they are handled by
        :func:`make_smash_symlinks`).
    """
    from .utils import create_index

    delve_dir = Path(delve_dir)
    exp_dir = delve_dir / 'exposures'
    exp_dir.mkdir(parents=True, exist_ok=True)

    print('-' * 68)
    print(' PREPARING DELVE MC EXPOSURES FOR PHOTRED EXPOSURE-LEVEL PROCESSING')
    print('-' * 68)

    allexpstr = fits.getdata(str(expfile), 1)
    nallexp = len(allexpstr)
    print(f"prep: {nallexp} DECam exposures loaded from {expfile}")

    # Strip whitespace from string columns
    for col in ('fluxfile', 'maskfile', 'wtfile', 'expnum'):
        if col in allexpstr.dtype.names:
            allexpstr[col][:] = np.char.strip(allexpstr[col].astype(str))

    # Fix /net/mss1 paths
    for col in ('fluxfile', 'maskfile', 'wtfile'):
        if col in allexpstr.dtype.names:
            allexpstr[col][:] = np.char.lstrip(
                allexpstr[col].astype(str), '/net')

    # Filter cuts: good filter and exptime >= 90 s
    filt = np.char.upper(
        np.char.strip(allexpstr['filter'].astype(str))[:, :1])
    exp_col = 'exposure' if 'exposure' in allexpstr.dtype.names else 'exptime'
    good = np.isin(np.char.lower(filt.astype(str)),
                   [f[0].lower() for f in _GOOD_FILTERS])
    good &= allexpstr[exp_col].astype(float) >= 90
    expstr = allexpstr[good]
    print(f"prep: {len(expstr)} exposures pass the cuts")

    # Assign nights
    night_arr = np.array([
        _date_to_night(str(r['date_obs']).strip()) for r in expstr
    ])

    # Nightmin cut
    if nightmin is not None:
        mask = np.array([int(n) >= int(nightmin) for n in night_arr],
                        dtype=bool)
        expstr = expstr[mask]
        night_arr = night_arr[mask]
        print(f"prep: {len(expstr)} exposures with night >= {nightmin}")

    # Get unique nights
    night_idx = create_index(night_arr)
    nnights = len(night_idx['value'])
    print(f"prep: {nnights} unique nights of data")

    # SMASH nights to skip
    smash_nights = set()
    if smash_dir is not None:
        smash_dir = Path(smash_dir)
        smash_nights = {
            d.name for d in smash_dir.iterdir()
            if d.is_dir() and d.name.isdigit() and len(d.name) == 8
        }

    # New-only filter
    if newonly:
        done_nights = {
            d.name for d in exp_dir.iterdir()
            if d.is_dir() and d.name.isdigit() and len(d.name) == 8
        }
        keep = np.array([
            night_idx['value'][k] not in done_nights
            for k in range(nnights)
        ])
        kept_values = night_idx['value'][keep]
        if len(kept_values) == 0:
            print("prep: no new nights to process")
            return
        print(f"prep: {len(kept_values)} new nights to process")
        keep_exps = np.zeros(len(expstr), dtype=bool)
        for v in kept_values:
            keep_exps |= (night_arr == v)
        expstr = expstr[keep_exps]
        night_arr = night_arr[keep_exps]
        night_idx = create_index(night_arr)
        nnights = len(night_idx['value'])

    # Build per-night tasks
    night_expfiles = []
    for n in range(nnights):
        inight = str(night_idx['value'][n])
        if inight in smash_nights:
            print(f"  {inight}: SMASH night, skipping")
            continue
        ind = night_idx['index'][night_idx['lo'][n]:night_idx['hi'][n] + 1]
        nexp_night = expstr[ind]

        night_dir = exp_dir / inight
        night_dir.mkdir(parents=True, exist_ok=True)
        night_expfile = night_dir / f'{inight}_exposures.fits'

        if night_expfile.exists() and not redo:
            # Possibly add new exposures
            existing = fits.getdata(str(night_expfile), 1)
            ex_expnums = set(np.char.strip(existing['expnum'].astype(str)))
            new_mask = ~np.isin(
                np.char.strip(nexp_night['expnum'].astype(str)),
                list(ex_expnums))
            if new_mask.any():
                from numpy.lib.recfunctions import stack_arrays
                combined = stack_arrays([existing, nexp_night[new_mask]],
                                        asrecarray=False, usemask=False,
                                        autoconvert=True)
                fits.BinTableHDU(combined).writeto(
                    str(night_expfile), overwrite=True)
                print(f"  {inight}: added {new_mask.sum()} exposures")
            else:
                print(f"  {inight}: already processed, skipping")
                continue
        else:
            fits.BinTableHDU(nexp_night).writeto(
                str(night_expfile), overwrite=True)

        night_expfiles.append(str(night_expfile))

    if not night_expfiles:
        print("prep: nothing to do")
        return

    # Run prep_night for each night
    kwargs = dict(
        delve_dir=delve_dir, delvereddir=delvereddir, decam_txt=decam_txt,
        scriptsdir=scriptsdir, irafdir=irafdir, workdir=workdir, redo=redo,
    )
    if nmulti > 1:
        with ProcessPoolExecutor(max_workers=nmulti) as pool:
            futures = [pool.submit(prep_night, f, **kwargs)
                       for f in night_expfiles]
            for fut in futures:
                fut.result()
    else:
        for f in night_expfiles:
            prep_night(f, **kwargs)
