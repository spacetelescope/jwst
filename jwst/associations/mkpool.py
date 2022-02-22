"""
Tools for pool creation
"""
from astropy.io.fits import getheader as fits_getheader
import numpy as np

from . import AssociationPool

# Header keywords to ignore
IGNORE_KEYS = ('', 'COMMENT', 'HISTORY')

# Non-header columns that need to be defined
NON_HEADER_COLS = {
    'dms_note': '',
    'is_imprt': 'f',
    'is_psf': 'f',
    'pntgtype': 'science',
}


def mkpool(data, asn_candidate=None, dms_note='', is_imprt='f', is_psf='f', pntgtype='science', **kwargs):
    """Make a pool from data

    A number of columns used by the Association rules cannot be derived from the header
    keywords. The columns, and typical other values, are as follows:

    - asn_candidate
      The observation candidate is always defined in table creation. However, higher level
      associations can be created by specifying a list of candidate definitions. An example
      of adding both background and coronographic candidates would be:
      ["('c1000', 'background')", "('c1001', 'coronographic')"]

    - dms_note
      Notes from upstream processing of the downlinked data that may be pertinent
      to the quality of the data. Currently the value "wfsc_los_jitter" is used
      by the Level 2 wavefront sensing rule, Asn_Lv2WFSC, to ignore exposures.

    - is_imprt
      A 't' indicates the exposure is an imprint exposure.

    - is_psf
      A 't' indicate a PSF exposure.

    - pntgtype
      The general class of exposure. The default value is "science".
      For target acquisition, the value is "target_acquisition".

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

    dms_note : str
        Value for the dms_note column.

    is_imprt : 't' or 'f'
        Indicator whether exposures are imprint/leakcal exposures.

    is_psf : 't' ro 'f'
        Indicator whether exposures are PSF's.

    pntgtype : 'science', 'target_acquisition'
        General exposure type.

    kwargs : dict
        Other keyword arguments to pass to the
        `astropy.io.fits.getheader` call.
    """
    params = set()
    for datum in data:
        params.update(getheader(datum))
    params.update(NON_HEADER_COLS)
    params.add('asn_candidate')

    params = params.difference(IGNORE_KEYS)
    params = [item.lower() for item in params]

    pool = AssociationPool(names=params, dtype=[object] * len(params))

    non_header_params = {'dms_note': dms_note, 'is_imprt': is_imprt, 'is_psf': is_psf, 'pntgtype': pntgtype}
    for datum in data:
        header = getheader(datum, **kwargs)
        valid_params = {
            keyword.lower(): value
            for keyword, value in header.items()
            if keyword not in IGNORE_KEYS
        }

        # Update non-header parameters
        valid_params.update(non_header_params)

        # Setup association candidates
        combinded_asn_candidates = [f"('o{header['observtn']}', 'observation')"]
        if asn_candidate is not None:
            combinded_asn_candidates += asn_candidate
        valid_params['asn_candidate'] = '[' + ','.join(combinded_asn_candidates) + ']'

        # Add the exposure
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
