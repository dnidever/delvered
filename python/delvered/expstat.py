"""
DELVERED exposure pipeline status reporting.

Python port of dlvexpstat.pro → exp_stat()
"""

import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from astropy.io import fits


# Pipeline stage names (in order)
STAGES = ['WCS', 'DAOPHOT', 'MATCH', 'APCOR', 'ASTROM',
          'ZEROPOINT', 'CALIB', 'COMBINE', 'DEREDDEN', 'SAVE']

_NIGHT_DTYPE = np.dtype([
    ('night', 'U8'),
    ('dir_mtime', np.float64),
    ('dir_size', np.int64),
    ('logsdir', np.int8),
    ('logsdir_mtime', np.float64),
    ('logsdir_size', np.int64),
    ('fieldsfile', np.int8),
    ('setupfile', np.int8),
    ('exposurefile', np.int8),
    ('nightsumfile', np.int8),
    ('wcs_success', np.int32),
    ('wcs_failure', np.int32),
    ('daophot_success', np.int32),
    ('daophot_failure', np.int32),
    ('match_success', np.int32),
    ('match_failure', np.int32),
    ('apcor_success', np.int32),
    ('apcor_failure', np.int32),
    ('astrom_success', np.int32),
    ('astrom_failure', np.int32),
    ('zeropoint_success', np.int32),
    ('zeropoint_failure', np.int32),
    ('calib_success', np.int32),
    ('calib_failure', np.int32),
    ('combine_success', np.int32),
    ('combine_failure', np.int32),
    ('deredden_success', np.int32),
    ('deredden_failure', np.int32),
    ('save_success', np.int32),
    ('save_failure', np.int32),
    ('success', np.int32, (10,)),
    ('failure', np.int32, (10,)),
    ('nexpall', np.int32),
    ('nexp', np.int32),
    ('done', np.uint8),
])


def _count_lines(path):
    """Return number of non-empty lines in *path*, or 0 if file absent/empty."""
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return 0
    with open(p) as fh:
        return sum(1 for ln in fh if ln.strip())


def _fits_nrows(fitsfile):
    """Return NAXIS2 from the first binary table extension header."""
    try:
        hdr = fits.getheader(str(fitsfile), ext=1)
        return int(hdr.get('NAXIS2', 0))
    except Exception:
        return 0


def exp_stat(delve_dir, redo=False, all_nights=False):
    """
    Report DELVERED pipeline status across all observed nights.

    Python port of dlvexpstat.pro.

    Scans the ``exposures/`` sub-tree for nightly pipeline log files
    (``logs/STAGE.success`` / ``logs/STAGE.failure``), counts successes
    and failures per pipeline stage, and writes a timestamped summary
    FITS file to ``exposures/summary/``.

    Parameters
    ----------
    delve_dir : str or Path
        Root DELVE data directory (contains ``exposures/``).
    redo : bool
        Re-compute every night even if nothing has changed since the
        last summary.
    all_nights : bool
        Print a status row even for nights with no log activity.

    Returns
    -------
    numpy structured array
        Summary table with one row per night.
    """
    delve_dir = Path(delve_dir)
    exp_dir = delve_dir / 'exposures'

    # Discover night directories (8-digit names like 20YYMMDD)
    night_dirs = sorted(
        d for d in exp_dir.iterdir()
        if d.is_dir() and d.name.isdigit() and len(d.name) == 8
    )
    ndirs = len(night_dirs)
    print(f"exp_stat: {ndirs} DELVE nights")
    if ndirs == 0:
        return np.zeros(0, dtype=_NIGHT_DTYPE)

    nights = [d.name for d in night_dirs]

    # Load previous summary if available
    summary_dir = exp_dir / 'summary'
    sumfiles = sorted(summary_dir.glob('delve_expsummary_*.fits')) if summary_dir.exists() else []
    last = np.zeros(ndirs, dtype=_NIGHT_DTYPE)
    for k, night in enumerate(nights):
        last['night'][k] = night
        last['nexpall'][k] = -1
        last['nexp'][k] = -1
        last['success'][k] = -1
        last['failure'][k] = 0

    if sumfiles and not redo:
        try:
            prev = fits.getdata(str(sumfiles[-1]), 1)
            prev_nights = np.char.strip(prev['night'].astype(str))
            for k, night in enumerate(nights):
                match = np.where(prev_nights == night)[0]
                if len(match):
                    for field in _NIGHT_DTYPE.names:
                        if field in prev.dtype.names:
                            last[field][k] = prev[field][match[0]]
        except Exception:
            pass

    str_arr = last.copy()
    if redo:
        str_arr = np.zeros(ndirs, dtype=_NIGHT_DTYPE)
        for k, night in enumerate(nights):
            str_arr['night'][k] = night
            str_arr['nexpall'][k] = -1
            str_arr['nexp'][k] = -1
            str_arr['success'][k] = -1

    ngivesum = 0
    for i, idir in enumerate(night_dirs):
        inight = idir.name
        str_arr['night'][i] = inight

        dir_stat = idir.stat()
        str_arr['dir_mtime'][i] = dir_stat.st_mtime
        str_arr['dir_size'][i] = dir_stat.st_size

        logs_dir = idir / 'logs'
        str_arr['logsdir'][i] = int(logs_dir.is_dir())
        if logs_dir.is_dir():
            ls = logs_dir.stat()
            str_arr['logsdir_mtime'][i] = ls.st_mtime
            str_arr['logsdir_size'][i] = ls.st_size

        # Check if something changed compared to previous summary
        changed = (
            str_arr['dir_mtime'][i] != last['dir_mtime'][i] or
            str_arr['dir_size'][i] != last['dir_size'][i] or
            str_arr['logsdir_mtime'][i] != last['logsdir_mtime'][i] or
            str_arr['logsdir_size'][i] != last['logsdir_size'][i] or
            redo or all_nights
        )

        if changed:
            str_arr['fieldsfile'][i] = int((idir / 'fields').exists())
            str_arr['setupfile'][i] = int((idir / 'photred.setup').exists())
            expfile = idir / f'{inight}_exposures.fits'
            str_arr['exposurefile'][i] = int(expfile.exists())
            nsumfile = idir / f'{inight}_summary.fits'
            str_arr['nightsumfile'][i] = int(nsumfile.exists())

            if expfile.exists():
                str_arr['nexpall'][i] = _fits_nrows(expfile)
            if nsumfile.exists():
                str_arr['nexp'][i] = _fits_nrows(nsumfile)

            if logs_dir.is_dir():
                for j, stage in enumerate(STAGES):
                    sval = _count_lines(logs_dir / f'{stage}.success')
                    fval = _count_lines(logs_dir / f'{stage}.failure')
                    str_arr[f'{stage.lower()}_success'][i] = sval
                    str_arr[f'{stage.lower()}_failure'][i] = fval
                    str_arr['success'][i][j] = sval
                    str_arr['failure'][i][j] = fval

        # Determine "done"
        s = str_arr['success'][i]
        f = str_arr['failure'][i]
        if (np.all(f == 0) and np.all(s > 0) and
                str_arr['fieldsfile'][i] == 1 and
                str_arr['setupfile'][i] == 1 and
                str_arr['exposurefile'][i] == 1 and
                str_arr['nightsumfile'][i] == 1):
            str_arr['done'][i] = 1

        # Print status row
        has_activity = (str_arr['logsdir'][i] == 1 and
                        np.max(str_arr['success'][i]) > 0)
        if has_activity or all_nights:
            if ngivesum % 40 == 0:
                hdr = ('NIGHT      NEXP     WCS     DAOPHOT    MATCH     APCOR'
                       '   ASTROM    ZEROPT  CALIB   COMBINE  DERED   SAVE  NIGHTSUM    DONE')
                print(hdr)
            nexp = (str_arr['nexp'][i] if str_arr['nexp'][i] > 0
                    else str_arr['nexpall'][i])
            comment = ''
            if str_arr['done'][i]:
                comment = '       FINISHED'
            elif np.sum(str_arr['failure'][i]) > 0:
                comment = '       !!!!'

            parts = [f'{inight:<11s}', f'{nexp:4d}']
            for j, stage in enumerate(STAGES):
                s = str_arr[f'{stage.lower()}_success'][i]
                f = str_arr[f'{stage.lower()}_failure'][i]
                parts.append(f'{s:5d}/{f:<4d}')
            parts.append(f'{str_arr["nightsumfile"][i]:5d}')
            parts.append(comment)
            print(' '.join(parts))
            ngivesum += 1

    done_mask = str_arr['done'] == 1
    ndone = int(done_mask.sum())
    print(f"{ndone}/{ndirs} finished")
    nexp_done = int(str_arr['nexp'][done_mask].clip(min=0).sum())
    nexp_all = int(str_arr['nexpall'].clip(min=0).sum())
    print(f"{nexp_done}/{nexp_all} exposures finished")

    # Write timestamped summary
    summary_dir.mkdir(parents=True, exist_ok=True)
    now = datetime.now(tz=timezone.utc)
    tag = now.strftime('%Y%m%d%H%M%S')
    outfile = summary_dir / f'delve_expsummary_{tag}.fits'
    print(f"exp_stat: writing summary to {outfile}")
    hdu = fits.BinTableHDU(str_arr)
    hdu.writeto(str(outfile), overwrite=True)

    return str_arr
