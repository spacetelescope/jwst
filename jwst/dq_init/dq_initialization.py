import logging

import numpy as np
from stdatamodels.jwst import datamodels

from jwst.datamodels import dqflags  # type: ignore[attr-defined]
from jwst.lib import reffile_utils

log = logging.getLogger(__name__)


# FGS guide star mode exposure types
guider_list = ["FGS_ID-IMAGE", "FGS_ID-STACK", "FGS_ACQ1", "FGS_ACQ2", "FGS_TRACK", "FGS_FINEGUIDE"]

__all__ = ["do_dqinit", "check_dimensions"]


def do_dqinit(output_model, mask_model):
    """
    Perform the dq_init step on a JWST datamodel.

    Parameters
    ----------
    output_model : `~stdatamodels.jwst.datamodels.RampModel` \
                   or `~stdatamodels.jwst.datamodels.GuiderRawModel`
        The JWST datamodel to be corrected.
    mask_model : `~stdatamodels.jwst.datamodels.MaskModel`
        The mask model to use in the correction.

    Returns
    -------
    output_model : `~stdatamodels.jwst.datamodels.RampModel` \
                   or `~stdatamodels.jwst.datamodels.GuiderRawModel`
        The corrected JWST datamodel, updated in place.
    """
    # Inflate empty DQ array, if necessary
    check_dimensions(output_model)

    # Extract subarray from reference data, if necessary
    if reffile_utils.ref_matches_sci(output_model, mask_model):
        mask_array = mask_model.dq
    else:
        log.info("Extracting mask subarray to match science data")
        mask_sub_model = reffile_utils.get_subarray_model(output_model, mask_model)
        mask_array = mask_sub_model.dq
        del mask_sub_model

    # Set model-specific data quality in output
    if output_model.meta.exposure.type in guider_list:
        output_model.dq |= mask_array
    else:
        output_model.pixeldq |= mask_array

        # Additionally, propagate mask DO_NOT_USE flags to groupdq to
        # ensure no ramps are fit to bad pixels.
        output_model.groupdq |= mask_array & dqflags.group["DO_NOT_USE"]

    output_model.meta.cal_step.dq_init = "COMPLETE"

    return output_model


def check_dimensions(input_model):
    """
    Check the input model pixeldq dimensions.

    The pixeldq attribute should have the same dimensions as
    the image plane of the input model science data
    If it has dimensions (0,0), create an array of zeros with the same shape
    as the image plane of the input model. For the FGS modes, the
    GuiderRawModel has only a regular dq array (no pixeldq or groupdq)

    Parameters
    ----------
    input_model : `~stdatamodels.jwst.datamodels.RampModel` \
                  or `~stdatamodels.jwst.datamodels.GuiderRawModel`
        Input datamodel.
    """
    input_shape = input_model.data.shape

    if isinstance(input_model, datamodels.GuiderRawModel):
        if input_model.dq.shape != input_shape[-2:]:
            # If the shape is different, then the mask model should have
            # a shape of (0,0).
            # If that's the case, create the array
            if input_model.dq.shape == (0, 0):
                input_model.dq = np.zeros(input_shape[-2:], dtype=np.uint32)
            else:
                log.error(f"DQ array has the wrong shape: {input_model.dq.shape}")

    else:  # RampModel
        if input_model.pixeldq.shape != input_shape[-2:]:
            # If the shape is different, then the mask model should have
            # a shape of (0,0).
            # If that's the case, create the array
            if input_model.pixeldq.shape == (0, 0):
                input_model.pixeldq = np.zeros(input_shape[-2:], dtype=np.uint32)
            else:
                log.error(f"Pixeldq array has the wrong shape: {input_model.pixeldq.shape}")

        # Perform the same check for the input model groupdq array
        if input_model.groupdq.shape != input_shape:
            if input_model.groupdq.shape == (0, 0, 0, 0):
                input_model.groupdq = np.zeros(input_shape, dtype=np.uint8)
            else:
                log.error(f"Groupdq array has wrong shape: {input_model.groupdq.shape}")
