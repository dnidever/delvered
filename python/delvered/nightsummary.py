"""
DELVERED per-night summary file construction.

Python port of delvered_nightsummary.pro → night_summary()
"""

from pathlib import Path

import numpy as np
from astropy.io import fits
from numpy.lib.recfunctions import stack_arrays, append_fields


def _load_fields_file(fields_file):
    """
    Parse the ``fields`` ASCII file with columns ``shname name``.

    Returns a list of (shname, name) tuples.
    """
    entries = []
    path = Path(fields_file)
    if not path.exists():
        return entries
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 2:
                entries.append((parts[0], parts[1]))
    return entries


def _add_str_field(arr, field, values):
    """Append a fixed-length string field to a structured array."""
    return append_fields(arr, field, values, usemask=False)


def night_summary(night, delve_dir, redo=False):
    """
    Build a per-night summary FITS from individual field summary files.

    Python port of delvered_nightsummary.pro.

    Reads the ``fields`` ASCII file in the night directory, then for
    each field loads ``{field}_summary.fits`` (extensions 1 and 2),
    appends ``fieldname`` and ``field`` columns, and concatenates all
    fields into a single nightly summary FITS file written as
    ``{night}_summary.fits``.

    Parameters
    ----------
    night : str
        Night identifier, e.g. ``'20160101'``.
    delve_dir : str or Path
        Root DELVE data directory (contains ``exposures/``).
    redo : bool
        Overwrite the output file if it already exists.

    Returns
    -------
    tuple of (numpy structured array or None, numpy structured array or None)
        ``(expstr, chipstr)`` — exposure-level and chip-level summary
        tables, or ``None`` if no data were found.
    """
    delve_dir = Path(delve_dir)
    exp_dir = delve_dir / 'exposures'
    night_dir = exp_dir / str(night)

    if not night_dir.exists():
        print(f"night_summary: {night_dir} NOT FOUND")
        return None, None

    nightsumfile = night_dir / f'{night}_summary.fits'
    if nightsumfile.exists() and not redo:
        try:
            nrows = fits.getheader(str(nightsumfile), ext=1).get('NAXIS2', 0)
        except Exception:
            nrows = 0
        if nrows > 0:
            print(f"night_summary: {nightsumfile} EXISTS and redo=False")
            expstr = fits.getdata(str(nightsumfile), 1)
            try:
                chipstr = fits.getdata(str(nightsumfile), 2)
            except Exception:
                chipstr = None
            return expstr, chipstr

    fields = _load_fields_file(night_dir / 'fields')
    print(f"night_summary: {len(fields)} fields")

    exp_parts = []
    chip_parts = []

    for shname, name in fields:
        sumfile = night_dir / f'{name}_summary.fits'
        if not sumfile.exists():
            print(f"night_summary: {sumfile} NOT FOUND")
            continue
        try:
            expstr1 = fits.getdata(str(sumfile), 1)
            chipstr1 = fits.getdata(str(sumfile), 2)
        except Exception as exc:
            print(f"night_summary: PROBLEM loading {sumfile}: {exc}")
            continue

        n1 = len(expstr1)
        n2 = len(chipstr1) if chipstr1 is not None else 0

        fname_arr = np.array([name] * n1)
        sname_arr = np.array([shname] * n1)
        expstr1 = append_fields(expstr1, ['fieldname', 'field'],
                                [fname_arr, sname_arr], usemask=False)
        exp_parts.append(expstr1)

        if n2 > 0:
            fname_arr2 = np.array([name] * n2)
            sname_arr2 = np.array([shname] * n2)
            chipstr1 = append_fields(chipstr1, ['fieldname', 'field'],
                                     [fname_arr2, sname_arr2], usemask=False)
            chip_parts.append(chipstr1)

    if not exp_parts:
        print(f"night_summary: no exposures for night={night}")
        return None, None

    expstr = stack_arrays(exp_parts, asrecarray=False, usemask=False,
                          autoconvert=True)
    chipstr = (stack_arrays(chip_parts, asrecarray=False, usemask=False,
                            autoconvert=True)
               if chip_parts else None)

    print(f"night_summary: writing {nightsumfile}")
    primary = fits.PrimaryHDU()
    hdus = fits.HDUList([primary, fits.BinTableHDU(expstr)])
    if chipstr is not None:
        hdus.append(fits.BinTableHDU(chipstr))
    hdus.writeto(str(nightsumfile), overwrite=True)

    return expstr, chipstr
