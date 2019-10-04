#
#  Module for subtracting a super-bias image from science data sets
#

import numpy as np
import logging
from ..lib import reffile_utils

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def do_correction(input_model, bias_model):
    """
    Short Summary
    -------------
    Execute all tasks for Super-Bias Subtraction

    Parameters
    ----------
    input_model: data model object
        science data to be corrected

    bias_model: super-bias model object
        bias data

    Returns
    -------
    output_model: data model object
        bias-subtracted science data

    """

    # Check for subarray mode and extract subarray from the
    # bias reference data if necessary
    if not reffile_utils.ref_matches_sci(input_model, bias_model):
        bias_model = reffile_utils.get_subarray_model(input_model, bias_model)

    # Replace NaN's in the superbias with zeros
    bias_model.data[np.isnan(bias_model.data)] = 0.0

    # Subtract the bias ref image from the science data
    output_model = subtract_bias(input_model, bias_model)

    output_model.meta.cal_step.superbias = 'COMPLETE'

    return output_model


def subtract_bias(input, bias):
    """
    Subtracts a superbias image from a science data set, subtracting the
    superbias from each group of each integration in the science data.
    The DQ flags in the bias reference image are propagated into the science
    data pixeldq array. The error array is unchanged.

    Parameters
    ----------
    input: data model object
        the input science data

    bias: superbias model object
        the superbias image data

    Returns
    -------
    output: data model object
        bias-subtracted science data

    """

    # Create output as a copy of the input science data model
    output = input.copy()

    # combine the science and superbias DQ arrays
    output.pixeldq = np.bitwise_or(input.pixeldq, bias.dq)

    # Subtract the superbias image from all groups and integrations
    # of the science data
    output.data -= bias.data

    return output
