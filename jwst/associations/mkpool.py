"""
Tools for pool creation
"""
import logging
from copy import copy

from astropy.io.fits import getheader as fits_getheader

from . import AssociationPool

__all__ = ['mkpool']

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())
LogLevels = [logging.WARNING, logging.INFO, logging.DEBUG]

# Header keywords to ignore
IGNORE_KEYS = ('', 'COMMENT', 'HISTORY')

# Non-header columns that need to be defined
NON_HEADER_COLS = {
    'asn_candidate': None,
    'dms_note': '',
    'is_imprt': 'f',
    'pntgtype': 'science',
    'targetid': '1',
}


def mkpool(data,
           asn_candidate=NON_HEADER_COLS['asn_candidate'], dms_note=NON_HEADER_COLS['dms_note'],
           is_imprt=NON_HEADER_COLS['is_imprt'], pntgtype=NON_HEADER_COLS['pntgtype'],
           **kwargs):
    """Create an association pool from a list of FITS files.

    Normally, association pools and the associations generated from those pools
    are created by the automatic ground system process. Users should download
    and modify those pools if need be. If desired, this function can be used to
    create pools from scratch using a list of FITS files. Once created, the
    :py:func:`~jwst.associations.generate` can be used to create associations
    from these pools.

    A number of pool columns used by the Association rules cannot be derived
    from the header keywords. The columns, and typical other values, are as
    follows:

    - asn_candidate
        The observation candidate is always defined in table
        creation, based on the observation id of each exposure.

        However, higher level associations can be created by specifying a list
        of candidate definitions. An example of adding both background and
        coronographic candidates would be: [('c1000', 'background'),
        ('c1001', 'coronographic')]

        The specification can be either as a list of 2-tuples, as presented above, or
        as a single string representation of the list. Using the previous example, the
        following is also a valid input:
        "[('c1000', 'background'), ('c1001', 'coronographic')]"
    - dms_note
          Notes from upstream processing of the downlinked data that may be pertinent
          to the quality of the data. Currently the value "wfsc_los_jitter" is used
          by the Level 2 wavefront sensing rule, Asn_Lv2WFSC, to ignore exposures.
    - is_imprt
          A 't' indicates the exposure is a NIRSpec imprint exposure.
    - pntgtype
          The general class of exposure. The default value is "science".
          For target acquisition, the value is "target_acquisition".

    Parameters
    ----------
    data : int
        The data to get the pool parameters from.
        Can be pathnames or `astropy.io.fits.HDUL`
        or `astropy.io.fits.ImageHDU`.

    asn_candidate : [(id, type)[,...]] or None
        Association candidates to add to each exposure.
        These are added to the default ('oXXX', 'observation') candidate
        created from header information.

    dms_note : str
        Value for the dms_note column.

    is_imprt : 't' or 'f'
        Indicator whether exposures are imprint/leakcal exposures.

    pntgtype : 'science', 'target_acquisition'
        General exposure type.

    kwargs : dict
        Other keyword arguments to pass to the
        `astropy.io.fits.getheader` call.

    Returns
    -------
    pool : `jwst.associations.AssociationPool`
        The association pool.

    """
    params = set()
    for datum in data:
        params.update(getheader(datum))
    params.update(NON_HEADER_COLS)

    params = params.difference(IGNORE_KEYS)
    params = [item.lower() for item in params]
    params.sort()
    defaults = {param: 'null' for param in params}

    pool = AssociationPool(names=params, dtype=[object] * len(params))

    # Set default values for user-settable non-header parameters
    non_header_params = {'dms_note': dms_note, 'is_imprt': is_imprt, 'pntgtype': pntgtype}

    # Setup for target id calculation
    targetid = 0  # Start off with no target id.
    target_names = set()

    # Create the table.
    for datum in data:
        header = getheader(datum, **kwargs)
        valid_params = {
            keyword.lower(): str(value)
            for keyword, value in header.items()
            if keyword not in IGNORE_KEYS
        }

        # Update non-header parameters
        valid_params.update(non_header_params)

        # Setup association candidates
        combined_asn_candidates = [(f"o{header['observtn']}", "observation")]
        if isinstance(asn_candidate, str):
            combined_asn_candidates = f'[{combined_asn_candidates[0]}, {asn_candidate[1:]}'
        else:
            if asn_candidate is not None:
                combined_asn_candidates += asn_candidate
            combined_asn_candidates = str(combined_asn_candidates)
        valid_params['asn_candidate'] = combined_asn_candidates

        # Calculate target id.
        if valid_params['targname'] not in target_names:
            target_names.add(valid_params['targname'])
            targetid += 1
        valid_params['targetid'] = str(targetid)

        # Add the exposure
        final_params = copy(defaults)
        final_params.update(valid_params)
        pool.add_row(final_params)

    return pool


def from_cmdline(args=None):
    """Collect command-line options and run mkpool

    Parameters
    ----------
    args : [str[,...]]
        List of arguments to parse

    Returns : dict
        Dict of the arguments and their values.
    """
    import argparse

    parser = argparse.ArgumentParser(
        description='Create an Association Pool file from a list of exposures'
    )

    parser.add_argument(
        'pool',
        help='Name of the pool file to save to.'
    )
    parser.add_argument(
        'data', nargs='+',
        help='List of exposures to create the Association Pool with.'
    )

    parser.add_argument(
        '--asn-candidate', default=NON_HEADER_COLS['asn_candidate'],
        help='Additional candidate information.'
    )
    parser.add_argument(
        '--dms-note', default=NON_HEADER_COLS['dms_note'],
        help='Added notes that may be relevant to association creation.'
    )
    parser.add_argument(
        '--is-imprt', default=NON_HEADER_COLS['is_imprt'],
        help='A "t" indicates the exposure is an imprint exposure.'
    )
    parser.add_argument(
        '--pntgtype', default=NON_HEADER_COLS['pntgtype'],
        help='The general class of exposure.'
    )
    parser.add_argument(
        '-v', '--verbose', action='count', default=0,
        help='Increase verbosity. Specifying multiple times adds more output.'
    )

    parsed = parser.parse_args(args)

    # Set output detail.
    level = LogLevels[min(len(LogLevels) - 1, parsed.verbose)]
    logger.setLevel(level)

    # That's all folks.
    mkpool_args = vars(parsed)
    del mkpool_args['verbose']
    return mkpool_args


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
