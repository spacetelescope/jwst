"""
Tools for pool creation
"""
from astropy.io.fits import getheader as fits_getheader
import numpy as np

from . import AssociationPool

IGNORE_KEYS = ('', 'COMMENT', 'HISTORY')


def mkpool(data, **kwargs):
    """Make a pool from data

    Parameters
    ----------
    data : [dataum[, ...]]
        The data to get the pool parameters from.
        Can be pathnames or `astropy.io.fits.HDUL`
        or `astropy.io.fits.ImageHDU

    kwargs : dict
        Other keyword arguments to pass to the
        `astropy.io.fits.getheader` call.
    """
    params = set()
    for datum in data:
        params.update(getheader(datum))

    params = params.difference(IGNORE_KEYS)

    pool = AssociationPool(names=params, dtype=[object] * len(params))

    for datum in data:
        valid_params = {
            keyword: value
            for keyword, value in getheader(datum, **kwargs).items()
            if keyword not in IGNORE_KEYS
        }
        pool.add_row(valid_params)

    return pool


def getheader(datum, **kwargs):
    """Get header from the data item

    Parameters
    ----------
    datum : str or HDUList or HDU
        Source of the header information

    kwargs : dict
        Keyword arguments passed to `astropy.io.fits.getheader`.
        Relevant ones are `ext`, `extname`, or `extver`
    """

    # Parse out HDU key
    try:
        key = kwargs['ext']
    except KeyError:
        try:
            key = (kwargs['extname'], kwargs.get('extver', 0))
        except KeyError:
            key = 0

    # Attempt to get an HDU
    try:
        hdu = datum[key]
    except TypeError:
        hdu = datum

    try:
        header = hdu.header
    except AttributeError:
        pass
    else:
        header['FILENAME'] = hdu.fileinfo()['file'].name
        return header

    header = fits_getheader(datum, **kwargs)
    header['FILENAME'] = datum
    return header
