"""
DAOPHOT/ALLSTAR pipeline wrapper for the DELVERED pipeline.

Python port of delvered_daophot.pro → run_daophot()
"""

import subprocess
import shutil
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np


def _read_setup(setup_file):
    """
    Parse a PHOTRED ``photred.setup`` file into a dict.

    Lines beginning with ``#`` are ignored.  Each key-value pair is
    split on the first run of whitespace.
    """
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
    """Return non-empty lines from a log/list file."""
    p = Path(path)
    if not p.exists():
        return []
    with open(p) as fh:
        return [ln.strip() for ln in fh if ln.strip()]


def _write_list(path, lines):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'a') as fh:
        fh.write('\n'.join(lines) + '\n')


def _als_ok(als_path, min_lines=3):
    """Return True if an ALS file exists and has enough data lines."""
    p = Path(als_path)
    if not p.exists():
        return False
    with open(p) as fh:
        return sum(1 for _ in fh) > min_lines


def _run_one_chip(args):
    """
    Run daophot + allstar on a single chip FITS file.

    Parameters
    ----------
    args : tuple
        (fits_file, scripts_dir, workdir, psfstars, daopsfva, daofitradfwhm)

    Returns
    -------
    tuple of (fits_file_str, success: bool)
    """
    fits_file, scripts_dir, workdir, psfstars, daopsfva, daofitradfwhm = args
    chip_dir = Path(fits_file).parent
    base = Path(fits_file).name
    for ext in ('.fits.fz', '.fits'):
        if base.endswith(ext):
            base = base[:-len(ext)]
            break

    # Copy scripts to chip directory
    scripts_dir = Path(scripts_dir)
    for script in ('daophot.sh', 'photo.opt', 'apcor.opt',
                   'lstfilter.py', 'goodpsf.pro', 'srcfilter.pro'):
        src = scripts_dir / script
        dst = chip_dir / script
        if src.exists() and not dst.exists():
            shutil.copy2(str(src), str(dst))

    # Run daophot.sh
    cmd = [str(chip_dir / 'daophot.sh'), base]
    if workdir:
        cmd.append(str(workdir))
    try:
        subprocess.run(cmd, cwd=str(chip_dir), check=True,
                       capture_output=True, timeout=600)
    except Exception as exc:
        print(f"  daophot.sh FAILED for {fits_file}: {exc}")
        return str(fits_file), False

    als_ok = _als_ok(chip_dir / f'{base}.als')
    aals_ok = _als_ok(chip_dir / f'{base}a.als')
    return str(fits_file), als_ok and aals_ok


def run_daophot(night_dir, redo=False, nmulti=10, workdir=None):
    """
    Run DAOPHOT and ALLSTAR on all chip FITS files for a DELVE night.

    Python port of delvered_daophot.pro.

    Reads the night's PHOTRED setup and ``logs/WCS.success`` list,
    calls ``daophot.sh`` for each chip file in parallel, and writes
    ``logs/DAOPHOT.success`` / ``logs/DAOPHOT.failure`` lists.

    Parameters
    ----------
    night_dir : str or Path
        Night directory containing ``photred.setup`` and ``logs/``.
    redo : bool
        Re-run files that already have ALS output.
    nmulti : int
        Number of parallel workers.
    workdir : str or Path, optional
        Temporary working directory for daophot.

    Returns
    -------
    tuple of (list, list)
        (success_files, failure_files)
    """
    night_dir = Path(night_dir)
    setup = _read_setup(night_dir / 'photred.setup')
    scripts_dir = setup.get('SCRIPTSDIR', '')
    if not scripts_dir or not Path(scripts_dir).is_dir():
        print(f"run_daophot: SCRIPTSDIR not found: {scripts_dir!r}")
        return [], []

    if workdir is None:
        workdir = setup.get('WORKDIR', '')
    psfstars = setup.get('PSFSTARS', '1') not in ('0', '')
    daopsfva = setup.get('DAOPSFVA', None)
    daofitradfwhm = setup.get('DAOFITRADFWHM', None)
    if not redo:
        redo_str = setup.get('REDO', '0')
        redo = redo_str not in ('0', '', '-1')
    nmulti_cfg = int(setup.get('NMULTI_DAOPHOT',
                               setup.get('NMULTI', str(nmulti))))
    nmulti = max(nmulti_cfg, 1)

    # Input: files from WCS.success (or WCS.inlist)
    input_files = []
    for logname in ('logs/WCS.success', 'logs/WCS.inlist'):
        lines = _read_list(night_dir / logname)
        if lines:
            input_files = lines
            break

    if not input_files:
        print("run_daophot: no input files found")
        return [], []

    # Resolve to absolute paths
    abs_files = []
    for f in input_files:
        p = Path(f)
        if not p.is_absolute():
            p = night_dir / p
        abs_files.append(str(p))

    # Skip already-done files unless redo
    to_process = []
    already_done = []
    for f in abs_files:
        p = Path(f)
        base = p.name
        for ext in ('.fits.fz', '.fits'):
            if base.endswith(ext):
                base = base[:-len(ext)]
                break
        chip_dir = p.parent
        if (not redo and
                _als_ok(chip_dir / f'{base}.als') and
                _als_ok(chip_dir / f'{base}a.als')):
            already_done.append(f)
        else:
            to_process.append(f)

    print(f"run_daophot: {len(to_process)} files to process, "
          f"{len(already_done)} already done")

    args_list = [
        (f, scripts_dir, workdir, psfstars, daopsfva, daofitradfwhm)
        for f in to_process
    ]

    results = []
    if nmulti > 1:
        with ProcessPoolExecutor(max_workers=nmulti) as pool:
            results = list(pool.map(_run_one_chip, args_list))
    else:
        results = [_run_one_chip(a) for a in args_list]

    success = [f for f, ok in results if ok]
    failure = [f for f, ok in results if not ok]
    success += already_done

    # Write log files
    logs_dir = night_dir / 'logs'
    logs_dir.mkdir(exist_ok=True)
    if success:
        _write_list(logs_dir / 'DAOPHOT.success', success)
    if failure:
        _write_list(logs_dir / 'DAOPHOT.failure', failure)

    # ALS outlist
    als_out = [str(Path(f).with_suffix('').with_suffix('') + '.als'
                   if f.endswith('.fits.fz')
                   else str(Path(f).with_suffix('.als')))
               for f in success]
    _write_list(logs_dir / 'DAOPHOT.outlist', als_out)

    print(f"run_daophot: {len(success)} success, {len(failure)} failure")
    return success, failure
