import logging

import numpy as np

from stdatamodels.jwst import datamodels

from jwst.lib import reffile_utils
from jwst.datamodels import dqflags  # type: ignore[attr-defined]

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


# FGS guide star mode exposure types
guider_list = ["FGS_ID-IMAGE", "FGS_ID-STACK", "FGS_ACQ1", "FGS_ACQ2", "FGS_TRACK", "FGS_FINEGUIDE"]


def correct_model(input_model, mask_model):
    """
    Perform the dq_init step on a JWST datamodel.

    Parameters
    ----------
    input_model : input JWST datamodel
        The jwst datamodel to be corrected
    mask_model : mask datamodel
        The mask model to use in the correction

    Returns
    -------
    output_model : JWST datamodel
        The corrected JWST datamodel
    """
    output_model = do_dqinit(input_model, mask_model)

    return output_model


def do_dqinit(output_model, mask_model):
    """
    Perform the dq_init step on a JWST datamodel.

    Parameters
    ----------
    output_model : input JWST datamodel
        The jwst datamodel to be corrected
    mask_model : mask datamodel
        The mask model to use in the correction

    Returns
    -------
    output_model : JWST datamodel
        The corrected JWST datamodel
    """
    # Inflate empty DQ array, if necessary
    check_dimensions(output_model)

    # Extract subarray from reference data, if necessary
    if reffile_utils.ref_matches_sci(output_model, mask_model):
        mask_array = mask_model.dq
    else:
        log.info("Extracting mask subarray to match science data")
        mask_sub_model = reffile_utils.get_subarray_model(output_model, mask_model)
        mask_array = mask_sub_model.dq.copy()
        mask_sub_model.close()

    # Set model-specific data quality in output
    if output_model.meta.exposure.type in guider_list:
        dq = np.bitwise_or(output_model.dq, mask_array)
        output_model.dq = dq
    else:
        dq = np.bitwise_or(output_model.pixeldq, mask_array)
        output_model.pixeldq = dq
        # Additionally, propagate mask DO_NOT_USE flags to groupdq to
        # ensure no ramps are fit to bad pixels.
        output_model.groupdq = np.bitwise_or(
            output_model.groupdq, mask_array & dqflags.group["DO_NOT_USE"]
        )

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
    input_model : JWST datamodel
        Input datamodel.
    """
    input_shape = input_model.data.shape

    if isinstance(input_model, datamodels.GuiderRawModel):
        if input_model.dq.shape != input_shape[-2:]:
            # If the shape is different, then the mask model should have
            # a shape of (0,0).
            # If that's the case, create the array
            if input_model.dq.shape == (0, 0):
                input_model.dq = np.zeros(input_shape[-2:]).astype("uint32")
            else:
                log.error(f"DQ array has the wrong shape: {input_model.dq.shape}")

    else:  # RampModel
        if input_model.pixeldq.shape != input_shape[-2:]:
            # If the shape is different, then the mask model should have
            # a shape of (0,0).
            # If that's the case, create the array
            if input_model.pixeldq.shape == (0, 0):
                input_model.pixeldq = np.zeros(input_shape[-2:]).astype("uint32")
            else:
                log.error(f"Pixeldq array has the wrong shape: {input_model.pixeldq.shape}")

        # Perform the same check for the input model groupdq array
        if input_model.groupdq.shape != input_shape:
            if input_model.groupdq.shape == (0, 0, 0, 0):
                input_model.groupdq = np.zeros(input_shape).astype("uint8")
            else:
                log.error(f"Groupdq array has wrong shape: {input_model.groupdq.shape}")
    return
