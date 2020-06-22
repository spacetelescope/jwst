#!/usr/bin/env python
"""Pointing verification

Check for some consistency in exposure files for pointing accuracy.
Specifically, compare the V1 and Reference RA/DEC to what the intended
TARG RA/DEC were specified.

Examples
--------
>>> pointing_verification exp1.fits

>>> pointing_verification *.fits
"""
import argparse
from collections import namedtuple
import logging

from astropy.coordinates import SkyCoord
from astropy.io import ascii
from astropy.table import Table
import astropy.units as u

import jwst.datamodels as dm

logger = logging.getLogger('jwst')
logger.addHandler(logging.StreamHandler())
LogLevels = [logging.WARNING, logging.INFO, logging.DEBUG]

# Basic delta structure
Delta = namedtuple('Delta', 'target, v1, refpoint, delta_v1, delta_refpoint')


def calc_pointing_deltas(model):
    """Calculate pointing deltas

    Parameters
    ----------
    model : jwst.datamodel.DataModel
        The model to check pointing information in.

    Returns
    -------
    delta : Delta
        A `namedtuple` with the following keys. If for some reason
        a value could not be determined, `None` is given.

        - 'target':         `SkyCoord` of the target
        - 'v1':             `SkyCoord` of V1
        - 'refpoint':       `SkyCoord` of the reference point
        - 'delta_v1':       Difference between V1 and proposed TARGET
        - 'delta_refpoint': Difference between reference pixel pointing and TARGET
    """
    # Retrieve the info from the model
    target = SkyCoord(model.meta.target.ra * u.degree, model.meta.target.dec * u.degree)
    v1 = SkyCoord(model.meta.pointing.ra_v1 * u.degree, model.meta.pointing.dec_v1 * u.degree)
    refpoint = SkyCoord(model.meta.wcsinfo.ra_ref * u.degree, model.meta.wcsinfo.dec_ref * u.degree)

    # Calculate separations
    delta = Delta(
        target=target,
        v1=v1,
        refpoint=refpoint,
        delta_v1=target.separation(v1),
        delta_refpoint=target.separation(refpoint),
    )
    return delta


def calc_deltas(exposures):
    """Create table of pointing deltas

    Parameters
    ----------
    exposures : [file-like[,...]] or [DataModel[,...]]
        List of file-like objects or `jwst.datamodels.DataModel` to retrieve
        pointing information from.

    Returns
    -------
    deltas : astropy.table.Table
        Table of results with the following columns:

        - exposure: The exposure the pointing information is from
        - target:         `SkyCoord` of the proposed target
        - v1:             `SkyCoord` of v1
        - refpoint:       `SkyCoord` of the reference point
        - delta_v1:       target - V1 separation
        - delta_refpoint: target - refpoint separation
    """
    # Initialize structures
    targets = list()
    v1s = list()
    refpoints = list()
    delta_v1s = list()
    delta_refpoints = list()

    # Calculate deltas for all input.
    for exposure in exposures:
        with dm.open(exposure) as model:
            delta = calc_pointing_deltas(model)
            logger.debug(f'{model}: delta v1={delta.delta_v1} delta refpoint={delta.delta_refpoint}')

            targets.append(delta.target)
            v1s.append(delta.v1)
            refpoints.append(delta.refpoint)
            delta_v1s.append(delta.delta_v1.degree)
            delta_refpoints.append(delta.delta_refpoint.degree)

    # Places results into a Table.
    deltas_dict = {
        'exposure':       exposures,
        'target':         targets,
        'v1':             v1s,
        'refpoint':       refpoints,
        'delta_v1':       delta_v1s,
        'delta_refpoint': delta_refpoints,
    }
    deltas = Table(deltas_dict)
    return deltas


# Begin execution
if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Compare various pointing information for consistency.'
    )

    parser.add_argument(
        'exposures', type=str, nargs='+',
        help='List of JWST data files to examine.'
    )
    parser.add_argument(
        '-v', '--verbose', action='count', default=0,
        help='Increase verbosity. Specifying multiple times adds more output.'
    )

    args = parser.parse_args()

    # Set output detail.
    level = LogLevels[min(len(LogLevels)-1, args.verbose)]
    logger.setLevel(level)

    # Process the file list.
    deltas = calc_deltas(args.exposures)
    ascii.write(deltas, format='ecsv')
