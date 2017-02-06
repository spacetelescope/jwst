"""
Tools for pool creation
"""

from astropy.io.fits import getheader
from . import AssociationPool

IGNORE_KEYS = ('', 'COMMENT', 'HISTORY')


def mkpool(data, *kwargs):
    """Make a pool from data

    Parameters
    ----------
    data: [dataum[, ...]]
        The data to get the pool parameters from.
        Can be pathnames or `astropy.io.fits.HDUL`
        or `astropy.io.fits.ImageHDU

    kwargs: dict
        Other keyword arguments to pass to the
        `astropy.io.fits.getheader` call.
    """
    params = set()
    for datum in data:
        params.update(getheader(datum))

    params = params.difference(IGNORE_KEYS)

    pool = AssociationPool(names=params, dtype=['S20'] * len(params))

    for datum in data:
        valid_params = {
            keyword: value
            for keyword, value in getheader(datum).items()
            if keyword not in IGNORE_KEYS
        }
        pool.add_row(valid_params)

    return pool
