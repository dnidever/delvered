"""
DAOMATCH / DAOMASTER wrapper for the DELVERED pipeline.

Python port of delvered_match.pro → run_match()
"""

import re
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np


def _read_setup(setup_file):
    cfg = {}
    p = Path(setup_file)
    if not p.exists():
        return cfg
    with open(p) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split(None, 1)
            if len(parts) == 2:
                cfg[parts[0].upper()] = parts[1].strip()
    return cfg


def _read_list(path):
    p = Path(path)
    if not p.exists():
        return []
    with open(p) as fh:
        return [ln.strip() for ln in fh if ln.strip()]


def _write_list(path, lines):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'a') as fh:
        fh.write('\n'.join(lines) + '\n')


def _mch_ok(mch_path, expected_lines):
    p = Path(mch_path)
    if not p.exists():
        return False
    with open(p) as fh:
        n = sum(1 for _ in fh)
    return n == expected_lines


def _raw_ok(raw_path):
    p = Path(raw_path)
    if not p.exists():
        return False
    with open(p) as fh:
        return sum(1 for _ in fh) > 3


def _get_filter(fits_file):
    """Return the filter from a FITS header."""
    try:
        from astropy.io import fits as _fits
        hdr = _fits.getheader(str(fits_file), ext=0)
        return str(hdr.get('FILTER', '')).strip()[:1]
    except Exception:
        return ''


def _get_exptime(fits_file):
    """Return exposure time from a FITS header."""
    try:
        from astropy.io import fits as _fits
        hdr = _fits.getheader(str(fits_file), ext=0)
        return float(hdr.get('EXPTIME', hdr.get('EXPOSURE', 0)))
    except Exception:
        return 0.0


def _pick_reference(als_bases, filtref_list, chip_dir):
    """
    Pick the reference image (longest exptime in preferred filter).

    Returns the index into *als_bases*.
    """
    best_idx = None
    best_exptime = -1.0
    for filt_pref in filtref_list:
        for i, base in enumerate(als_bases):
            fits_file = chip_dir / f'{base}.fits'
            if not fits_file.exists():
                fits_file = chip_dir / f'{base}.fits.fz'
            if not fits_file.exists():
                continue
            f = _get_filter(fits_file)
            if f.lower() == filt_pref.lower():
                et = _get_exptime(fits_file)
                if et > best_exptime:
                    best_exptime = et
                    best_idx = i
        if best_idx is not None:
            return best_idx
    # Fallback: just pick first
    return 0


def _run_daomatch_chip(args):
    """
    Run daomatch (+ daomaster) for one chip in a field.

    Parameters
    ----------
    args : tuple
        (chip_dir, als_files, ref_base, mch_base, filtref_list, maxshift)

    Returns
    -------
    (mch_file, success: bool, als_list)
    """
    chip_dir, als_files, ref_base, mch_base, filtref_list, maxshift = args
    chip_dir = Path(chip_dir)
    nals = len(als_files)

    if nals < 2:
        return None, False, als_files

    # Build DAOMATCH input: first file is reference, rest follow
    inlist_content = '\n'.join(als_files) + '\n'
    inlist_file = chip_dir / f'.daomatch_in_{mch_base}.tmp'
    inlist_file.write_text(inlist_content)

    mch_file = chip_dir / f'{mch_base}.mch'

    # Build daomatch command
    cmd = ['daomatch']
    if maxshift:
        cmd += ['-maxshift', str(maxshift)]
    cmd += [als_files[0]]  # reference first

    try:
        with open(inlist_file) as stdin_fh:
            result = subprocess.run(
                cmd, cwd=str(chip_dir), stdin=stdin_fh,
                capture_output=True, timeout=300)
    except Exception as exc:
        print(f"  daomatch FAILED {mch_base}: {exc}")
        inlist_file.unlink(missing_ok=True)
        return str(mch_file), False, als_files

    inlist_file.unlink(missing_ok=True)
    ok = _mch_ok(mch_file, nals) and _raw_ok(chip_dir / f'{mch_base}.raw')
    return str(mch_file), ok, als_files


def run_match(night_dir, redo=False, nmulti=10):
    """
    Run DAOMATCH for all chip groups in a DELVE night.

    Python port of delvered_match.pro.

    Groups ALS files by chip (from ``logs/DAOPHOT.success``), picks a
    reference exposure per chip using the preferred filter list from
    ``photred.setup``, calls ``daomatch`` in parallel, and writes
    ``logs/MATCH.success``, ``logs/MATCH.failure``, and
    ``logs/MATCH.outlist``.

    Parameters
    ----------
    night_dir : str or Path
        Night directory.
    redo : bool
        Re-run groups that already have MCH output.
    nmulti : int
        Number of parallel workers.

    Returns
    -------
    tuple of (list, list)
        (success_als_files, failure_als_files)
    """
    night_dir = Path(night_dir)
    setup = _read_setup(night_dir / 'photred.setup')
    filtref_raw = setup.get('FILTREF', 'g,i,r,z,u')
    filtref_list = [f.strip() for f in filtref_raw.split(',')]
    maxshift = setup.get('MCHMAXSHIFT', None)
    if maxshift in (None, '0', '', '-1'):
        maxshift = None
    nmulti_cfg = int(setup.get('NMULTI_MATCH',
                               setup.get('NMULTI', str(nmulti))))
    nmulti = max(nmulti_cfg, 1)
    if not redo:
        redo_str = setup.get('REDO', '0')
        redo = redo_str not in ('0', '', '-1')

    # Input: ALS files from DAOPHOT.success
    input_als = _read_list(night_dir / 'logs/DAOPHOT.success')
    if not input_als:
        # Try outlist
        input_als = _read_list(night_dir / 'logs/DAOPHOT.outlist')
    if not input_als:
        print("run_match: no DAOPHOT ALS input files found")
        return [], []

    # Resolve to absolute paths
    abs_als = []
    for f in input_als:
        p = Path(f)
        if not p.is_absolute():
            p = night_dir / p
        abs_als.append(str(p))

    # Group by chip directory
    chip_groups = {}  # chip_dir_str → [als_base_list]
    for als_path in abs_als:
        p = Path(als_path)
        chip_dir = str(p.parent)
        base = p.stem  # e.g. F5-00423440_34
        chip_groups.setdefault(chip_dir, []).append(base)

    args_list = []
    for chip_dir_str, bases in chip_groups.items():
        chip_dir = Path(chip_dir_str)
        if len(bases) < 2:
            continue
        # Pick reference image
        ref_idx = _pick_reference(bases, filtref_list, chip_dir)
        ref_base = bases[ref_idx]
        # Re-order: reference first
        ordered = [bases[ref_idx]] + [b for i, b in enumerate(bases) if i != ref_idx]
        als_files = [b + '.als' for b in ordered]
        mch_base = ref_base  # MCH named after reference

        if not redo and _mch_ok(chip_dir / f'{mch_base}.mch', len(als_files)):
            continue  # already done

        args_list.append((chip_dir_str, als_files, ref_base, mch_base,
                          filtref_list, maxshift))

    print(f"run_match: {len(args_list)} chip groups to match")

    if not args_list:
        return abs_als, []

    if nmulti > 1:
        with ProcessPoolExecutor(max_workers=nmulti) as pool:
            results = list(pool.map(_run_daomatch_chip, args_list))
    else:
        results = [_run_daomatch_chip(a) for a in args_list]

    success_mch = []
    success_als = []
    failure_als = []
    for mch_file, ok, als_in in results:
        chip_dir = str(Path(mch_file).parent) if mch_file else ''
        if ok:
            success_mch.append(mch_file)
            success_als += [str(Path(chip_dir) / f) if chip_dir else f
                            for f in als_in]
        else:
            failure_als += [str(Path(chip_dir) / f) if chip_dir else f
                            for f in als_in]

    # Write log files
    logs_dir = night_dir / 'logs'
    logs_dir.mkdir(exist_ok=True)
    if success_als:
        _write_list(logs_dir / 'MATCH.success', success_als)
    if failure_als:
        _write_list(logs_dir / 'MATCH.failure', failure_als)
    if success_mch:
        _write_list(logs_dir / 'MATCH.outlist', success_mch)

    print(f"run_match: {len(success_mch)} groups matched, "
          f"{len([r for r in results if not r[1]])} failed")
    return success_als, failure_als
