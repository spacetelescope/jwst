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
import logging

from astropy.coordinates import SkyCoord
import astropy.units as u

import jwst.datamodels as dm

logger = logging.getLogger('jwst')
handler = logging.StreamHandler()
logger.addHandler(handler)
logger.setLevel(logging.INFO)


def calc_pointing_deltas(model):
    """Calculate pointing deltas

    Parameters
    ----------
    model : jwst.datamodel.DataModel
        The model to check pointing information in.

    Returns
    -------
    deltas : dict
        A dictionary with the following keys. If for some reason
        a value could not be determined, `None` is given.

        - 'v1': Difference between V1 and proposed TARGET
        - 'refpoint': Difference between reference pixel pointing and TARGET
    """
    # Retrieve the info from the model
    targ = SkyCoord(model.meta.target.ra * u.degree, model.meta.target.dec * u.degree)
    v1 = SkyCoord(model.meta.pointing.ra_v1 * u.degree, model.meta.pointing.dec_v1 * u.degree)
    refpoint = SkyCoord(model.meta.wcsinfo.ra_ref * u.degree, model.meta.wcsinfo.dec_ref * u.degree)

    # Calculate separations
    deltas = dict()
    deltas['v1'] = targ.separation(v1)
    deltas['refpoint'] = targ.separation(refpoint)

    return deltas


# Begin execution
if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Compare various pointing information for consistency'
    )

    parser.add_argument(
        'exposures', type=str, nargs='+',
        help='List of JWST data files to examine.'
    )

    args = parser.parse_args()

    # Process the file list
    for path in args.exposures:
        with dm.open(path) as model:
            deltas = calc_pointing_deltas(model)
            logger.info(f'{path}: delta v1={deltas["v1"]} delta refpoint={deltas["refpoint"]}')
