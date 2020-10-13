"""Pointing Summary

Review contents of a set of given models for pointing information.
Compare the calculated V1 and REFPOINT pointing with the proposed
TARGET pointing.

Examples
--------
>>> from jwst.datamodels import ImageModel
>>> im = ImageModel()
>>> im.meta.target.ra       =  90.75541666666666
>>> im.meta.target.dec      = -66.56055555555554
>>> im.meta.pointing.ra_v1  =  91.08142004561715
>>> im.meta.pointing.dec_v1 = -66.60547868904696
>>> im.meta.wcsinfo.ra_ref  =  90.70377653291781
>>> im.meta.wcsinfo.dec_ref = -66.59540223936895
>>> calc_pointing_deltas(im)
    Delta(target=<SkyCoord (ICRS): (ra, dec) in deg
        (90.75541667, -66.56055556)>, v1=<SkyCoord (ICRS): (ra, dec) in deg
        (91.08142005, -66.60547869)>, refpoint=<SkyCoord (ICRS): (ra, dec) in deg
        (90.70377653, -66.59540224)>, delta_v1=<Angle 0.13712727 deg>, delta_refpoint=<Angle 0.04044315 deg>)

>>> calc_deltas([im])
    <Table length=1>
      exposure                  target                ...    delta_refpoint
                               deg,deg                ...
       object                   object                ...       float64
    ------------ ------------------------------------ ... -------------------
    <ImageModel> 90.75541666666666,-66.56055555555554 ... 0.04044314761499765
"""
from collections import defaultdict, namedtuple
import logging

from astropy.coordinates import SkyCoord
from astropy.table import Table
import astropy.units as u

import jwst.datamodels as dm

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

__all__ = ['Delta', 'calc_pointing_deltas', 'calc_deltas']

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


def calc_deltas(exposures, extra_meta=None):
    """Create table of pointing deltas

    Parameters
    ----------
    exposures : [file-like[,...]] or [DataModel[,...]]
        List of file-like objects or `jwst.datamodels.DataModel` to retrieve
        pointing information from.

    extra_meta: [str[,...]] or None
        List of model meta attributes to add to the table.

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
    extra_meta_values = defaultdict(list)

    extra_meta = extra_meta if extra_meta is not None else list()

    # Calculate deltas for all input.
    for exposure in exposures:
        with dm.open(exposure) as model:
            delta = calc_pointing_deltas(model)
            logger.info(f'{model}: delta v1={delta.delta_v1} delta refpoint={delta.delta_refpoint}')

            targets.append(delta.target)
            v1s.append(delta.v1)
            refpoints.append(delta.refpoint)
            delta_v1s.append(delta.delta_v1.degree)
            delta_refpoints.append(delta.delta_refpoint.degree)

            for meta in extra_meta:
                extra_meta_values[meta].append(model[meta])


    # Places results into a Table.
    deltas_dict = {
        'exposure':       exposures,
        'target':         targets,
        'v1':             v1s,
        'refpoint':       refpoints,
        'delta_v1':       delta_v1s,
        'delta_refpoint': delta_refpoints,
    }
    deltas_dict.update(extra_meta_values)
    deltas = Table(deltas_dict)
    return deltas
