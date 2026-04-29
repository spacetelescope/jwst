import logging

import numpy as np
from stdatamodels.jwst import datamodels

from jwst.datamodels import dqflags  # type: ignore[attr-defined]
from jwst.lib import reffile_utils
from jwst.lib.exposure_types import FGS_GUIDE_EXP_TYPES

log = logging.getLogger(__name__)


__all__ = ["do_dqinit", "check_dimensions"]


def do_dqinit(output_model, mask_model, user_dq=None):
    """
    Perform the dq_init step on a JWST datamodel.

    Parameters
    ----------
    output_model : `~stdatamodels.jwst.datamodels.RampModel` or \
                   `~stdatamodels.jwst.datamodels.SuperstripeRampModel` or \
                   `~stdatamodels.jwst.datamodels.GuiderRawModel`
        The JWST datamodel to be corrected.
    mask_model : `~stdatamodels.jwst.datamodels.MaskModel`
        The mask model to use in the correction.
    user_dq : ndarray or None
        User-supplied DQ int array, if any.

    Returns
    -------
    output_model : `~stdatamodels.jwst.datamodels.RampModel` or \
                   `~stdatamodels.jwst.datamodels.SuperstripeRampModel` or \
                   `~stdatamodels.jwst.datamodels.GuiderRawModel`
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

    if user_dq is not None:
        if user_dq.shape != mask_array.shape:
            errmsg = f"user_dq has shape={user_dq.shape} but expecting {mask_array.shape}"
            log.error(errmsg)
            raise ValueError(errmsg)

        user_dq = user_dq.astype(mask_array.dtype)
        mask_array |= user_dq

    # Set model-specific data quality in output
    num_superstripe = getattr(output_model.meta.subarray, "num_superstripe", None)
    if str(output_model.meta.exposure.type).lower() in FGS_GUIDE_EXP_TYPES:
        output_model.dq |= mask_array

    elif num_superstripe is not None and num_superstripe > 0:
        # Store 3-D DQ array in pixeldq.
        output_model.pixeldq |= mask_array

        # Generate 4-D groupdq mask_array from pixeldq array, given output groupdq shape
        nints, ngroups, _, _ = output_model.groupdq.shape
        nsci_ints = nints // num_superstripe
        mask_array = mask_array[:, np.newaxis, :, :].repeat(ngroups, axis=1)
        mask_array = np.tile(mask_array, reps=(nsci_ints, 1, 1, 1))
        output_model.groupdq |= mask_array & dqflags.group["DO_NOT_USE"]
    else:
        output_model.pixeldq |= mask_array

        # Additionally, propagate mask DO_NOT_USE flags to groupdq to
        # ensure no ramps are fit to bad pixels.
        output_model.groupdq |= mask_array & dqflags.group["DO_NOT_USE"]

    output_model.meta.cal_step.dq_init = "COMPLETE"

    return output_model


def check_dimensions(input_model):
    """
    Check the input model DQ array dimensions.

    For `~stdatamodels.jwst.datamodels.GuiderRawModel`, compare
    the ``dq`` array to the expected default dimensions. For
    `~stdatamodels.jwst.datamodels.RampModel` or
    `~stdatamodels.jwst.datamodels.SuperstripeRampModel`,
    compare ``pixeldq`` and ``groupdq``.

    If the dimensions do not match, replace the current array with
    a default zero-filled one of the correct dimensions.

    Parameters
    ----------
    input_model : `~stdatamodels.jwst.datamodels.RampModel` or \
                  `~stdatamodels.jwst.datamodels.SuperstripeRampModel`, or \
                  `~stdatamodels.jwst.datamodels.GuiderRawModel`
        Input datamodel; updated in place.
    """
    input_shape = input_model.data.shape
    if isinstance(input_model, datamodels.GuiderRawModel):
        default_dq = input_model.get_default("dq")
        if input_model.dq is None or input_model.dq.shape != default_dq.shape:
            log.warning("Setting input DQ to default array")
            input_model.dq = default_dq

    else:  # RampModel or SuperstripeRampModel
        default_pixeldq = input_model.get_default("pixeldq")
        if input_model.pixeldq is None or input_model.pixeldq.shape != default_pixeldq.shape:
            log.warning("Setting input pixel DQ to default array")
            input_model.pixeldq = default_pixeldq

        # Perform the same check for the input model groupdq array
        # This array should always match the input data
        if input_model.groupdq is None or input_model.groupdq.shape != input_shape:
            log.warning("Setting input group DQ to default array")
            input_model.groupdq = input_model.get_default("groupdq")
