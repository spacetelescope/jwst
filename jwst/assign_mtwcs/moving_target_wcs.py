"""
Adjust the WCS of a moving target exposure.

Computes the average RA and DEC of a moving
target in all exposures in an association and adds a step to
each of the WCS pipelines to allow aligning the exposures to the average
location of the target.

"""
import logging
from copy import deepcopy
import numpy as np
from astropy.modeling.models import Shift, Identity
from gwcs import WCS
from gwcs import coordinate_frames as cf

from stdatamodels.jwst import datamodels

from jwst.datamodels import ModelContainer

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

__all__ = ["assign_moving_target_wcs"]


def assign_moving_target_wcs(input_model):

    if not isinstance(input_model, ModelContainer):
        raise ValueError("Expected a ModelContainer object")

    skip = False
    ra_sum = 0.0
    dec_sum = 0.0
    n = 0
    for model in input_model:
        if model.meta.asn.exptype.lower() != 'science':
            continue
        ra = model.meta.wcsinfo.mt_ra
        dec = model.meta.wcsinfo.mt_dec
        if ra is None or dec is None:
            skip = True
        if skip:
            break
        ra_sum += ra
        dec_sum += dec
        n += 1

    if skip:
        log.warning("One or more MT RA/Dec values missing in input images")
        log.warning("Step will be skipped, resulting in target misalignment")
    else:
        if n == 0:
            # TODO should this be an error? It was not previously
            # raise ValueError("assign_moving_target_wcs called with no science exposures")
            return input_model
        mt_avra = ra_sum / n
        mt_avdec = dec_sum / n

    for model in input_model:
        if model.meta.asn.exptype.lower() != 'science':
            continue
        if skip:
            model.meta.cal_step.assign_mtwcs = 'SKIPPED'
            continue

        model.meta.wcsinfo.mt_avra = mt_avra
        model.meta.wcsinfo.mt_avdec = mt_avdec
        if isinstance(model, datamodels.MultiSlitModel):
            for ind, slit in enumerate(model.slits):
                new_wcs = add_mt_frame(slit.meta.wcs,
                                       mt_avra, mt_avdec,
                                       slit.meta.wcsinfo.mt_ra, slit.meta.wcsinfo.mt_dec)
                del model.slits[ind].meta.wcs
                model.slits[ind].meta.wcs = new_wcs
        else:

            new_wcs = add_mt_frame(model.meta.wcs, mt_avra, mt_avdec,
                                   model.meta.wcsinfo.mt_ra, model.meta.wcsinfo.mt_dec)
            del model.meta.wcs
            model.meta.wcs = new_wcs

        model.meta.cal_step.assign_mtwcs = 'COMPLETE'
    return input_model


def add_mt_frame(wcs, ra_average, dec_average, mt_ra, mt_dec):
    """ Add a "moving_target" frame to the WCS pipeline.

    Parameters
    ----------
    wcs : `~gwcs.WCS`
        WCS object for the observation or slit.
    ra_average : float
        The average RA of all observations.
    dec_average : float
        The average DEC of all observations.
    mt_ra, mt_dec : float
        The RA, DEC of the moving target in the observation.

    Returns
    -------
    new_wcs : `~gwcs.WCS`
        The WCS for the moving target observation.
    """
    pipeline = wcs._pipeline[:-1]

    mt = deepcopy(wcs.output_frame)
    mt.name = 'moving_target'

    rdel = ra_average - mt_ra
    ddel = dec_average - mt_dec

    if isinstance(mt, cf.CelestialFrame):
        transform_to_mt = Shift(rdel) & Shift(ddel)
    elif isinstance(mt, cf.CompositeFrame):
        transform_to_mt = Shift(rdel) & Shift(ddel) & Identity(1)
    else:
        raise ValueError("Unrecognized coordinate frame.")

    pipeline.append((
        wcs.output_frame, transform_to_mt))
    pipeline.append((mt, None))
    new_wcs = WCS(pipeline)
    return new_wcs
