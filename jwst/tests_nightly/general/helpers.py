"""Helpers for the tests"""

from os.path import splitext


def add_suffix(fname, suffix, range=None):
    """Add suffix to file name

    Parameters
    ----------
    fname: str
        The file name to add the suffix to

    suffix: str
        The suffix to add_suffix

    range: range
        If specified, the set of indexes will be added to the
        outputs.

    Returns
    -------
    fname, fname_with_suffix
        2-tuple of the original file name and name with suffix.
        If `range` is defined, `fname_with_suffix` will be a list.

    """
    fname_root, fname_ext = splitext(fname)
    if range is None:
        with_suffix = ''.join([
            fname_root,
            '_',
            suffix,
            fname_ext
        ])
    else:
        with_suffix = []
        for idx in range:
            with_suffix.append(''.join([
                fname_root,
                '_',
                str(idx),
                '_',
                suffix,
                fname_ext
            ]))

    return fname, with_suffix
