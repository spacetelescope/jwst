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

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

__all__ = ["assign_moving_target_wcs"]


def assign_moving_target_wcs(input_models: ModelLibrary) -> ModelLibrary:

    with input_models:
        # get the indices of the science exposures in the ModelLibrary
        indices = input_models.ind_asn_type('science')

        mt_ra = []
        mt_dec = []
        for i in indices:
            sci_model = input_models.borrow(i)
            mt_ra.append(sci_model.meta.wcsinfo.mt_ra)
            mt_dec.append(sci_model.meta.wcsinfo.mt_dec)
            input_models.shelve(sci_model, i, modify=False)

        if None in mt_ra or None in mt_dec:
            log.warning("One or more MT RA/Dec values missing in input images")
            log.warning("Step will be skipped, resulting in target misalignment")
            for i in indices:
                sci_model = input_models.borrow(i)
                sci_model.meta.cal_step.assign_mtwcs = 'SKIPPED'
                input_models.shelve(sci_model, i, modify=True)
            return input_models

        mt_avra = np.mean(mt_ra)
        mt_avdec = np.mean(mt_dec)

        for i in indices:
            sci_model = input_models.borrow(i)
            sci_model.meta.wcsinfo.mt_avra = mt_avra
            sci_model.meta.wcsinfo.mt_avdec = mt_avdec
            if isinstance(sci_model, datamodels.MultiSlitModel):
                for ind, slit in enumerate(sci_model.slits):
                    new_wcs = add_mt_frame(slit.meta.wcs,
                                           mt_avra, mt_avdec,
                                           slit.meta.wcsinfo.mt_ra, slit.meta.wcsinfo.mt_dec)
                    del sci_model.slits[ind].meta.wcs
                    sci_model.slits[ind].meta.wcs = new_wcs
            else:
                new_wcs = add_mt_frame(sci_model.meta.wcs, mt_avra, mt_avdec,
                                       sci_model.meta.wcsinfo.mt_ra, sci_model.meta.wcsinfo.mt_dec)
                del sci_model.meta.wcs
                sci_model.meta.wcs = new_wcs
            sci_model.meta.cal_step.assign_mtwcs = 'COMPLETE'
            input_models.shelve(sci_model, i, modify=True)

    return input_models


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
