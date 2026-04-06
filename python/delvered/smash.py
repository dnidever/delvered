"""
SMASH → DELVE symlink creation utilities.

Python port of:
  - make_smash_symlinks.pro       → make_smash_symlinks()
  - make_smash_symlinks_night.pro → make_smash_symlinks_night()
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


# Default PHOTRED setup template written into each DELVE night directory
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
filtref     g,i,r,z,u
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

# DECam chip-number parser: extracts the trailing NN from names like
# F5-00423440_34 → 34
_CHIP_RE = re.compile(r'_(\d{2})$')


def _chip_num_from_base(base):
    """Extract the DECam chip number from a file basename (no extension)."""
    m = _CHIP_RE.search(base)
    return int(m.group(1)) if m else None


def _load_fields_file(path):
    """Return list of (shname, fullname) pairs from an ASCII ``fields`` file."""
    entries = []
    p = Path(path)
    if not p.exists():
        return entries
    with open(p) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 2:
                entries.append((parts[0], parts[1]))
    return entries


def _read_log(path):
    """Return non-empty lines from a log file, or [] if absent."""
    p = Path(path)
    if not p.exists():
        return []
    with open(p) as fh:
        return [ln.strip() for ln in fh if ln.strip()]


def _corner_radec_from_header(hdr):
    """Return corner RA/DEC arrays from a FITS header using WCS."""
    wcs = WCS(hdr)
    nx = int(hdr.get('NAXIS1', hdr.get('ZNAXIS1', 2048)))
    ny = int(hdr.get('NAXIS2', hdr.get('ZNAXIS2', 4096)))
    corners = np.array([[0, 0], [nx - 1, 0], [nx - 1, ny - 1], [0, ny - 1]],
                       dtype=float)
    sky = wcs.pixel_to_world(corners[:, 0], corners[:, 1])
    return sky.ra.deg, sky.dec.deg


def make_smash_symlinks_night(night, smash_dir, delve_dir, delvereddir,
                              instcal_list, decam_txt, redo=False,
                              scriptsdir='', irafdir='', workdir='/tmp'):
    """
    Create DELVE directory structure and symlinks for one processed SMASH night.

    Python port of make_smash_symlinks_night.pro.

    Mirrors the processed SMASH chip files (ALS, PSF, OPT, etc.) into the
    DELVE exposure tree using symbolic links, creates resource files that
    point to the original CP compressed FITS images, trims the reference
    catalogue to each chip footprint, and writes PHOTRED list files so
    that the DELVE pipeline can pick up from the WCS stage.

    Parameters
    ----------
    night : str
        Night identifier (8-digit, e.g. ``'20160101'``).
    smash_dir : str or Path
        Root of the SMASH processed data (contains per-night subdirectories).
    delve_dir : str or Path
        Root of the DELVE exposure tree (output).
    delvereddir : str or Path
        Root of the delvered project tree (for ``params/`` and ``data/``).
    instcal_list : str or Path
        Path to the DECam instcal master list FITS file.
    decam_txt : str or Path
        Path to ``data/decam.txt`` mapping extension names to ccdnums.
    redo : bool
        Re-create symlinks even if they already exist.
    scriptsdir : str
        Path written into ``photred.setup``.
    irafdir : str
        Path written into ``photred.setup``.
    workdir : str
        Working directory written into ``photred.setup``.

    Returns
    -------
    None
    """
    smash_dir = Path(smash_dir)
    delve_dir = Path(delve_dir)
    delvereddir = Path(delvereddir)
    modeleqnfile = str(delvereddir / 'params' / 'modelmag_equations.txt')

    night_src = smash_dir / night
    night_dst = delve_dir / night

    if not night_src.is_dir():
        print(f"make_smash_symlinks_night: {night_src} NOT FOUND")
        return

    print(f"\nMaking SMASH symlinks for night {night}")
    print('-' * 41)

    night_dst.mkdir(parents=True, exist_ok=True)

    # Load master DECam instcal list
    expstr = fits.getdata(str(instcal_list), 1)
    expstr_key = np.array(
        [f"{str(e).strip()}-{str(p).strip()}"
         for e, p in zip(expstr['expnum'], expstr['plver'])]
    )

    # Load DECam extension-name → ccdnum mapping
    decam = np.genfromtxt(str(decam_txt), names=True, dtype=None, encoding='ascii')
    decam_name_map = {str(r['name']).strip(): int(r['ccdnum']) for r in decam}

    # Load fields file
    fields = _load_fields_file(night_src / 'fields')

    # Collect FITS files from WCS.success and DAOPHOT.success logs
    fits_files = []
    for logname in ('WCS.success', 'DAOPHOT.success'):
        fits_files += _read_log(night_src / 'logs' / logname)

    # Normalise paths to be relative to night_src
    def _rel(fpath):
        p = Path(fpath)
        night_str = str(night)
        s = str(p)
        pos = s.find(night_str)
        if pos >= 0:
            return s[pos + len(night_str) + 1:]
        return fpath

    fits_files = list(dict.fromkeys(_rel(f) for f in fits_files))  # unique, ordered

    # Build field index
    field_map = {}  # short_name → [rel_paths]
    for rel in fits_files:
        base = Path(rel).name
        short = base.split('-')[0]
        field_map.setdefault(short, []).append(rel)

    # Load / build reference catalogue per field
    refcat_dir = night_dst / 'refcat'
    refcat_dir.mkdir(parents=True, exist_ok=True)

    wcs_lines = []
    daophot_lines = []
    match_lines = []
    apcor_lines = []

    for short_name, rel_files in field_map.items():
        # Find full field name from fields file
        long_name = None
        for sn, fn in fields:
            if sn == short_name:
                long_name = fn
                break
        if long_name is None:
            print(f"  {short_name} NOT FOUND in fields file")
            continue

        # Field summary for centre coordinates
        sumfile = night_src / f'{long_name}_summary.fits'
        if not sumfile.exists():
            print(f"  {sumfile} NOT FOUND")
            continue
        try:
            sumstr = fits.getdata(str(sumfile), 1)
            cen_ra = float(np.median(sumstr['ra']))
            cen_dec = float(np.median(sumstr['dec']))
        except Exception:
            cen_ra, cen_dec = 0.0, 0.0

        # Reference catalogue
        refcat_file = refcat_dir / f'{short_name}_refcat.fits.gz'
        if refcat_file.exists() and not redo:
            try:
                with gzip.open(str(refcat_file)) as gz:
                    refcat = fits.getdata(gz, 1)
            except Exception:
                refcat = None
        else:
            refcat = None
            print(f"  WARNING: no reference catalogue for {short_name}; "
                  f"run get_ref_cat({cen_ra:.4f}, {cen_dec:.4f}, 1.5) manually")

        # Get unique exposure numbers for this field
        expnums_seen = {}
        for rel in rel_files:
            base_noext = Path(rel).name
            for ext in ('.fits.fz', '.fits'):
                if base_noext.endswith(ext):
                    base_noext = base_noext[:-len(ext)]
            m = re.match(r'^([^-]+)-(\d{8})_(\d{2})$', base_noext)
            if m:
                expnum = m.group(2)
                if expnum not in expnums_seen:
                    expnums_seen[expnum] = rel

        field_dir = night_dst / short_name
        field_dir.mkdir(parents=True, exist_ok=True)

        print(f"  {short_name}/{long_name}: {len(expnums_seen)} exposures")

        for expnum, sample_rel in expnums_seen.items():
            # Get pipeline version from FITS header
            src_fits = night_src / sample_rel
            if not src_fits.exists() and not (night_src / (sample_rel + '.fz')).exists():
                print(f"    {expnum}: source FITS not found")
                continue
            if not src_fits.exists():
                src_fits = night_src / (sample_rel + '.fz')

            try:
                hdr0 = fits.getheader(str(src_fits), ext=0)
                plver = str(hdr0.get('PLVER', '')).strip()
            except Exception:
                plver = ''

            key = f'{expnum}-{plver}'
            idx = np.where(expstr_key == key)[0]
            if len(idx) == 0:
                print(f"    {expnum} not in master DECam list")
                continue
            iexp = idx[0]
            flux_file = str(expstr['fluxfile'][iexp]).strip().replace('/net/mss1/', '/mss1/')
            mask_file = str(expstr['maskfile'][iexp]).strip().replace('/net/mss1/', '/mss1/')
            wt_file = str(expstr['wtfile'][iexp]).strip().replace('/net/mss1/', '/mss1/')

            # Process all chips for this exposure
            chip_rels = [r for r in rel_files
                         if re.search(rf'-{expnum}_', Path(r).name)]
            for rel in chip_rels:
                base_stem = Path(rel).name
                for ext in ('.fits.fz', '.fits'):
                    if base_stem.endswith(ext):
                        base_stem = base_stem[:-len(ext)]
                chip_num = _chip_num_from_base(base_stem)
                if chip_num is None:
                    continue

                chip_dir = field_dir / f'chip{chip_num:02d}'
                chip_dir.mkdir(parents=True, exist_ok=True)

                src_opt = night_src / Path(rel).parent / f'{base_stem}.opt'
                if not src_opt.exists():
                    print(f"    {base_stem} NOT FOUND")
                    continue

                # Symlink DAOPHOT output files
                suffixes = ['_cat.dat', '.opt', '.als.opt', '.als', '.ap',
                            '.coo', '.plst', '.psf', '.psf.log', '.log',
                            'a.als', 'a.ap']
                src_base_dir = night_src / Path(rel).parent
                for suf in suffixes:
                    src_f = src_base_dir / f'{base_stem}{suf}'
                    dst_f = chip_dir / f'{base_stem}{suf}'
                    if dst_f.exists() or dst_f.is_symlink():
                        dst_f.unlink()
                    if src_f.exists():
                        dst_f.symlink_to(src_f)

                # Resource FITS stub
                res_fits = chip_dir / f'{base_stem}.fits'
                if res_fits.exists() or res_fits.is_symlink():
                    res_fits.unlink()
                res_fits.write_text('')

                # Resource dot-file with CP archive paths
                # Find extension name from decam mapping
                ext_name = next(
                    (n for n, c in decam_name_map.items() if c == chip_num),
                    str(chip_num)
                )
                dot_res = chip_dir / f'.{base_stem}.fits'
                dot_res.write_text(
                    f'fluxfile = {flux_file}[{ext_name}]\n'
                    f'wtfile = {wt_file}[{ext_name}]\n'
                    f'maskfile = {mask_file}[{ext_name}]\n'
                )

                # Trim reference catalogue to chip footprint
                if refcat is not None:
                    try:
                        with fits.open(str(src_fits)) as hdul:
                            chdr = hdul[ext_name].header if ext_name in [
                                h.name for h in hdul] else hdul[1].header
                        # Patch compressed sizes
                        if 'ZNAXIS1' in chdr:
                            chdr['NAXIS1'] = chdr['ZNAXIS1']
                            chdr['NAXIS2'] = chdr['ZNAXIS2']
                        vra, vdec = _corner_radec_from_header(chdr)
                        chip_ref_mask = (
                            (refcat['ra'] >= np.min(vra) - 0.02) &
                            (refcat['ra'] <= np.max(vra) + 0.02) &
                            (refcat['dec'] >= np.min(vdec) - 0.02) &
                            (refcat['dec'] <= np.max(vdec) + 0.02)
                        )
                        chip_ref = refcat[chip_ref_mask]
                        chip_ref_file = chip_dir / f'{base_stem}_refcat.fits'
                        fits.BinTableHDU(chip_ref).writeto(
                            str(chip_ref_file), overwrite=True)
                        subprocess.run(['gzip', '-f', str(chip_ref_file)],
                                       check=False)
                    except Exception as exc:
                        print(f"    WARNING: refcat trim failed for {base_stem}: {exc}")

                rel_chip = f'{short_name}/chip{chip_num:02d}/{base_stem}.fits'
                wcs_lines.append(rel_chip)
                daophot_lines.append(rel_chip)

        # MCH / RAW / TFR files
        mch_files = list((night_src / long_name).glob(f'{short_name}-????????_??.mch'))
        for mch in mch_files:
            stem = mch.stem
            chip_m = _chip_num_from_base(stem)
            if chip_m is None:
                continue
            chip_dst = field_dir / f'chip{chip_m:02d}'
            chip_dst.mkdir(parents=True, exist_ok=True)
            for suf in ('.mch', '.raw', '.tfr'):
                src = mch.with_suffix(suf)
                dst = chip_dst / src.name
                if dst.exists() or dst.is_symlink():
                    dst.unlink()
                if src.exists():
                    dst.symlink_to(src)
            match_lines.append(f'{short_name}/chip{chip_m:02d}/{mch.name}')

    # Copy apcor.lst
    src_apcor = night_src / 'apcor.lst'
    if src_apcor.exists():
        import shutil
        shutil.copy2(str(src_apcor), str(night_dst / 'apcor.lst'))
        apcor_raw = _read_log(night_src / 'logs' / 'APCOR.success')
        for line in apcor_raw:
            rel = _rel(line)
            base_stem = Path(rel).name
            for ext in ('.fits.fz', '.fits'):
                if base_stem.endswith(ext):
                    base_stem = base_stem[:-len(ext)]
            chip_num = _chip_num_from_base(base_stem)
            short = base_stem.split('-')[0]
            if chip_num is not None:
                apcor_lines.append(
                    f'{short}/chip{chip_num:02d}/{base_stem}.fits')

    # Write log files
    logs_dir = night_dst / 'logs'
    logs_dir.mkdir(parents=True, exist_ok=True)

    def _write_lines(path, lines):
        with open(path, 'w') as fh:
            fh.write('\n'.join(lines) + '\n')

    _write_lines(logs_dir / 'WCS.inlist', wcs_lines)
    _write_lines(logs_dir / 'DAOPHOT.success', daophot_lines)
    _write_lines(logs_dir / 'MATCH.outlist', match_lines)
    _write_lines(logs_dir / 'MATCH.success',
                 [ln.replace('.fits', '.als') for ln in daophot_lines])
    _write_lines(logs_dir / 'APCOR.success', apcor_lines)

    # Copy fields file
    src_fields = night_src / 'fields'
    if src_fields.exists():
        import shutil
        shutil.copy2(str(src_fields), str(night_dst / 'fields'))

    # Write photred.setup
    setup_text = _PHOTRED_SETUP.format(
        scriptsdir=scriptsdir, irafdir=irafdir,
        modeleqnfile=modeleqnfile, workdir=workdir)
    (night_dst / 'photred.setup').write_text(setup_text)

    print(f"  Wrote {len(wcs_lines)} chip entries")


def make_smash_symlinks(smash_dir, delve_dir, delvereddir, instcal_list,
                        decam_txt, redo=False, nmulti=3,
                        scriptsdir='', irafdir='', workdir='/tmp'):
    """
    Create DELVE symlinks for all processed SMASH nights.

    Python port of make_smash_symlinks.pro.

    Discovers all 8-digit night directories under *smash_dir* and calls
    :func:`make_smash_symlinks_night` for each, optionally in parallel.

    Parameters
    ----------
    smash_dir : str or Path
        Root of the SMASH processed data.
    delve_dir : str or Path
        Root of the DELVE exposure tree (output).
    delvereddir : str or Path
        Root of the delvered project tree.
    instcal_list : str or Path
        Master DECam instcal FITS list.
    decam_txt : str or Path
        ``data/decam.txt`` mapping extension names to ccdnums.
    redo : bool
        Re-create symlinks even if they exist.
    nmulti : int
        Number of parallel workers (default 3).
    scriptsdir, irafdir, workdir : str
        Paths written into ``photred.setup``.
    """
    smash_dir = Path(smash_dir)
    nights = sorted(
        d.name for d in smash_dir.iterdir()
        if d.is_dir() and d.name.isdigit() and len(d.name) == 8
    )
    print(f"{'#' * 38}")
    print(f"{len(nights)} SMASH nights to create symlinks for")
    print(f"{'#' * 38}")

    kwargs = dict(
        smash_dir=smash_dir, delve_dir=delve_dir,
        delvereddir=delvereddir, instcal_list=instcal_list,
        decam_txt=decam_txt, redo=redo,
        scriptsdir=scriptsdir, irafdir=irafdir, workdir=workdir,
    )

    if nmulti > 1:
        with ProcessPoolExecutor(max_workers=nmulti) as pool:
            futures = [pool.submit(make_smash_symlinks_night, n, **kwargs)
                       for n in nights]
            for f in futures:
                f.result()
    else:
        for night in nights:
            make_smash_symlinks_night(night, **kwargs)
