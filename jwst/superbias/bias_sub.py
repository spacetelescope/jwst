"""Module for subtracting a super-bias image from science data sets."""

import numpy as np
import logging
from jwst.lib import reffile_utils

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def do_correction(input_model, bias_model):
    """
    Execute all tasks for Super-Bias Subtraction.

    Parameters
    ----------
    input_model : RampModel
        Science data to be corrected
    bias_model : SuperBiasModel
        Bias data

    Returns
    -------
    output_model : stdatamodels.jwst.datamodels.ramp.RampModel
        Bias-subtracted science data
    """
    # Check for subarray mode and extract subarray from the
    # bias reference data if necessary
    if not reffile_utils.ref_matches_sci(input_model, bias_model):
        bias_model = reffile_utils.get_subarray_model(input_model, bias_model)

    # Replace NaN's in the superbias with zeros
    bias_model.data[np.isnan(bias_model.data)] = 0.0

    # Subtract the bias ref image from the science data
    output_model = subtract_bias(input_model, bias_model)

    output_model.meta.cal_step.superbias = "COMPLETE"

    return output_model


def subtract_bias(output, bias):
    """
    Subtract the superbias image from each group of each integration in the science data.

    The DQ flags in the bias reference image are propagated into the science
    data pixel DQ array. The error array is unchanged.

    Parameters
    ----------
    output : stdatamodels.jwst.datamodels.ramp.RampModel
        Input science data
    bias : stdatamodels.jwst.datamodels.superbias.SuperBiasModel
        Superbias image data

    Returns
    -------
    output : stdatamodels.jwst.datamodels.ramp.RampModel
        Bias-subtracted science data
    """
    # combine the science and superbias DQ arrays
    output.pixeldq = np.bitwise_or(output.pixeldq, bias.dq)

    # Subtract the superbias image from all groups and integrations
    # of the science data
    if not isinstance(type(output.data), float):
        output.data = (output.data).astype(float)
    output.data -= bias.data

    # If ZEROFRAME is present, subtract the super bias.  Zero values
    # indicate bad data, so should be kept zero.
    if output.meta.exposure.zero_frame:
        wh_zero = np.where(output.zeroframe == 0.0)
        output.zeroframe -= bias.data
        output.zeroframe[wh_zero] = 0.0  # Zero values indicate unusable data

    return output
