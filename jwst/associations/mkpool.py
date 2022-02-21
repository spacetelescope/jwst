"""
Tools for pool creation
"""
from astropy.io.fits import getheader as fits_getheader
import numpy as np

from . import AssociationPool

IGNORE_KEYS = ('', 'COMMENT', 'HISTORY')


def mkpool(data, asn_candidate=None, **kwargs):
    """Make a pool from data

    Parameters
    ----------
    data : [datum[, ...]]
        The data to get the pool parameters from.
        Can be pathnames or `astropy.io.fits.HDUL`
        or `astropy.io.fits.ImageHDU

    asn_candidate : ['(type, id)'[,...]] or None
        Association candidates to add to each exposure.
        These are added to the default ('observation', 'oXXX') candidate
        created from header information.

    kwargs : dict
        Other keyword arguments to pass to the
        `astropy.io.fits.getheader` call.
    """
    params = set()
    for datum in data:
        params.update(getheader(datum))
    params.add('asn_candidate')

    params = params.difference(IGNORE_KEYS)
    params = [item.lower() for item in params]

    pool = AssociationPool(names=params, dtype=[object] * len(params))

    for datum in data:
        header = getheader(datum, **kwargs)
        valid_params = {
            keyword.lower(): value
            for keyword, value in header.items()
            if keyword not in IGNORE_KEYS
        }

        # Setup association candidates
        combinded_asn_candidates = [f"('observation', 'o{header['observtn']}')"]
        if asn_candidate is not None:
            combinded_asn_candidates += asn_candidate
        valid_params['asn_candidate'] = '[' + ','.join(combinded_asn_candidates) + ']'
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
    header['FILENAME'] = str(datum)
    return header
