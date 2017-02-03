"""
Tools for pool creation
"""

from . import AssociationPool


def mkpool(data):
    """Make a pool from data

    Parameters
    ----------
    data: [dataum[, ...]]
        The data to get the pool parameters from.
        Can be pathnames or `astropy.io.fits.HDUL`
        or `astropy.io.fits.ImageHDU
    """
    pool = AssociationPool()

    return pool
