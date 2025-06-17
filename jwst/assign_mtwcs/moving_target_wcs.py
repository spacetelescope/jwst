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

from jwst.datamodels import ModelLibrary
from jwst.stpipe.utilities import record_step_status
from jwst.assign_wcs.util import update_s_region_imaging
from jwst.lib.exposure_types import IMAGING_TYPES

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

__all__ = ["assign_moving_target_wcs"]


def assign_moving_target_wcs(input_models):
    """
    Adjust the WCS of a moving target exposure.

    Parameters
    ----------
    input_models : `~jwst.datamodels.ModelLibrary`
        A collection of data models.

    Returns
    -------
    `~jwst.datamodels.ModelLibrary`
        The modified data models.
    """
    if not isinstance(input_models, ModelLibrary):
        raise TypeError(f"Expected a ModelLibrary object, not {type(input_models)}")

    # loop over only science exposures in the ModelLibrary
    ind = input_models.indices_for_exptype("science")
    mt_ra = np.full(len(ind), np.nan)
    mt_dec = np.full(len(ind), np.nan)
    mt_valid = True
    with input_models:
        for i in ind:
            model = input_models.borrow(i)
            if model.meta.wcsinfo.mt_ra is None or model.meta.wcsinfo.mt_dec is None:
                mt_valid = False
            else:
                mt_ra[i] = model.meta.wcsinfo.mt_ra
                mt_dec[i] = model.meta.wcsinfo.mt_dec

            if mt_valid and isinstance(model, datamodels.MultiSlitModel):
                for slit in model.slits:
                    if slit.meta.wcsinfo.mt_ra is None or slit.meta.wcsinfo.mt_dec is None:
                        mt_valid = False

            input_models.shelve(model, i, modify=False)

    # Compute the mean MT RA/Dec over all exposures
    if not mt_valid:
        log.warning("One or more MT RA/Dec values missing in input images.")
        log.warning("Step will be skipped, resulting in target misalignment.")
        record_step_status(input_models, "assign_mtwcs", False)
        return input_models

    mt_avra = mt_ra.mean()
    mt_avdec = mt_dec.mean()

    with input_models:
        for i in ind:
            model = input_models.borrow(i)
            model.meta.wcsinfo.mt_avra = mt_avra
            model.meta.wcsinfo.mt_avdec = mt_avdec
            if isinstance(model, datamodels.MultiSlitModel):
                for ind, slit in enumerate(model.slits):
                    new_wcs = add_mt_frame(
                        slit.meta.wcs,
                        mt_avra,
                        mt_avdec,
                        slit.meta.wcsinfo.mt_ra,
                        slit.meta.wcsinfo.mt_dec,
                    )
                    del model.slits[ind].meta.wcs
                    model.slits[ind].meta.wcs = new_wcs
            else:
                new_wcs = add_mt_frame(
                    model.meta.wcs,
                    mt_avra,
                    mt_avdec,
                    model.meta.wcsinfo.mt_ra,
                    model.meta.wcsinfo.mt_dec,
                )
                del model.meta.wcs
                model.meta.wcs = new_wcs
            if model.meta.exposure.type.lower() in IMAGING_TYPES:
                update_s_region_imaging(model)
            record_step_status(model, "assign_mtwcs", True)
            input_models.shelve(model, i, modify=True)

    return input_models


def add_mt_frame(wcs, ra_average, dec_average, mt_ra, mt_dec):
    """
    Add a "moving_target" frame to the WCS pipeline.

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
    pipeline = wcs.pipeline[:-1]

    mt = deepcopy(wcs.output_frame)
    mt.name = "moving_target"

    rdel = ra_average - mt_ra
    ddel = dec_average - mt_dec

    if isinstance(mt, cf.CelestialFrame):
        transform_to_mt = Shift(rdel) & Shift(ddel)
    elif isinstance(mt, cf.CompositeFrame):
        transform_to_mt = Shift(rdel) & Shift(ddel) & Identity(1)
    else:
        raise TypeError(f"Unrecognized coordinate frame type {type(mt)}.")

    pipeline.append((wcs.output_frame, transform_to_mt))
    pipeline.append((mt, None))
    new_wcs = WCS(pipeline)
    return new_wcs
