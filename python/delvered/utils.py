"""
DELVERED utility functions.

Python port of add_elements.pro and remove_tags.pro from the IDL DELVERED pipeline.
"""

import numpy as np
from numpy.lib.recfunctions import drop_fields, stack_arrays


def add_elements(arr, n_new):
    """
    Expand a numpy structured array by appending n_new empty elements.

    Python port of add_elements.pro.

    Parameters
    ----------
    arr : numpy structured array
        Input structured array to expand.
    n_new : int
        Number of new (zero-initialised) elements to append.

    Returns
    -------
    numpy structured array
        Expanded array with len(arr) + n_new elements.  The first
        len(arr) elements are copies of the original; the last n_new
        elements are zero/empty-initialised.
    """
    n_orig = len(arr)
    new_arr = np.zeros(n_orig + n_new, dtype=arr.dtype)
    new_arr[:n_orig] = arr
    return new_arr


def remove_tags(arr, tagnames):
    """
    Remove one or more fields from a numpy structured array.

    Python port of remove_tags.pro.

    Parameters
    ----------
    arr : numpy structured array
        Input structured array.
    tagnames : str or list of str
        Name(s) of fields to remove (case-insensitive).

    Returns
    -------
    numpy structured array or -1
        New structured array without the specified fields.
        Returns -1 if all fields would be removed.
        Returns the original array unchanged if no tags matched.
    """
    if isinstance(tagnames, str):
        tagnames = [tagnames]

    # normalise to upper case for case-insensitive matching
    tagnames_upper = [t.upper() for t in tagnames]
    existing = [f.upper() for f in arr.dtype.names]

    to_remove = [name for name in existing if name in tagnames_upper]

    if not to_remove:
        print("remove_tags: no tags matched, structure unchanged")
        return arr

    remaining = [name for name in arr.dtype.names if name.upper() not in tagnames_upper]
    if not remaining:
        print("remove_tags: all tags removed, returning -1")
        return -1

    # Use numpy recfunctions to drop the fields
    result = drop_fields(arr, [name for name in arr.dtype.names if name.upper() in tagnames_upper])
    return result


def create_index(values):
    """
    Create an index structure similar to IDL's CREATE_INDEX.

    Groups array indices by unique value, returning a dict with:
      'value' : sorted unique values
      'index' : flat array of all original indices, sorted by value
      'lo'    : start position in 'index' for each unique value
      'hi'    : end position (inclusive) in 'index' for each unique value

    Parameters
    ----------
    values : array-like
        1-D array of values to index.

    Returns
    -------
    dict
        Index structure.
    """
    values = np.asarray(values)
    sort_idx = np.argsort(values, kind='stable')
    sorted_vals = values[sort_idx]
    unique_vals, first_occ, counts = np.unique(sorted_vals, return_index=True, return_counts=True)

    lo = first_occ
    hi = first_occ + counts - 1

    return {
        'value': unique_vals,
        'index': sort_idx,
        'lo': lo,
        'hi': hi,
        'count': counts,
    }


def struct_assign(src, dst):
    """
    Copy matching fields from src to dst (in-place), like IDL STRUCT_ASSIGN.

    Parameters
    ----------
    src : numpy structured array or dict
        Source data.
    dst : numpy structured array
        Destination array.  Modified in-place.
    """
    src_names = set(src.dtype.names) if hasattr(src, 'dtype') else set(src.keys())
    dst_names = set(dst.dtype.names)
    common = src_names & dst_names
    for name in common:
        try:
            dst[name] = src[name]
        except (ValueError, TypeError):
            pass
