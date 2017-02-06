"""
Tools for pool creation
"""
from os.path import basename

from astropy.io.fits import getheader
import numpy as np

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

    pool = AssociationPool(names=params, dtype=[(np.str_, 20)] * len(params))

    for datum in data:
        valid_params = {
            keyword: value
            for keyword, value in getheader(datum).items()
            if keyword not in IGNORE_KEYS
        }
        valid_params['FILENAME'] = basename(datum)
        pool.add_row(valid_params)

    return pool
