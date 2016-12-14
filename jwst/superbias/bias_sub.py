from __future__ import division

#
#  Module for subtracting a super-bias image from science data sets
#

import numpy as np
import logging
from .. import datamodels

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

    # Replace NaN's in the superbias with zeros
    bias_model.data[np.isnan(bias_model.data)] = 0.0

    # Check for subarray mode and extract subarray from the
    # bias reference data if necessary
    if not ref_matches_sci(bias_model, input_model):
        bias_model = get_subarray(bias_model, input_model)

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


def ref_matches_sci(ref_model, sci_model):

    # See if the reference and science model subarray parameters match
    if (ref_model.meta.subarray.xstart == sci_model.meta.subarray.xstart and
        ref_model.meta.subarray.xsize == sci_model.meta.subarray.xsize and
        ref_model.meta.subarray.ystart == sci_model.meta.subarray.ystart and
        ref_model.meta.subarray.ysize == sci_model.meta.subarray.ysize):
        return True
    else:
        return False


def get_subarray(ref_model, sci_model):

    sci_x1 = sci_model.meta.subarray.xstart
    sci_x2 = sci_model.meta.subarray.xsize
    sci_y1 = sci_model.meta.subarray.ystart
    sci_y2 = sci_model.meta.subarray.ysize
    log.debug("sci xstart=%d, xsize=%d, ystart=%d, ysize=%d" % (sci_x1, sci_x2,
               sci_y1, sci_y2))
    ref_x1 = ref_model.meta.subarray.xstart
    ref_x2 = ref_model.meta.subarray.xsize
    ref_y1 = ref_model.meta.subarray.ystart
    ref_y2 = ref_model.meta.subarray.ysize
    if ref_x1 is None or ref_x2 is None or ref_y1 is None or ref_y2 is None:
        if ref_model.data.shape[-1] == 2048 and ref_model.data.shape[-2] == 2048:
            ref_x1 = 1
            ref_x2 = 2048
            ref_y1 = 1
            ref_y2 = 2048
        else:
            log.error("Missing subarray corner/size keywords in ref file")
            raise ValueError("Can't determine ref file readout properties")
    log.debug("ref xstart=%d, xsize=%d, ystart=%d, ysize=%d" % (ref_x1, ref_x2,
               ref_y1, ref_y2))

    # Compute the slicing indexes
    xstart = sci_x1 - ref_x1
    ystart = sci_y1 - ref_y1
    xstop = xstart + sci_model.meta.subarray.xsize
    ystop = ystart + sci_model.meta.subarray.ysize
    log.debug("slice xstart=%d, xstop=%d, ystart=%d, ystop=%d" % (xstart, xstop, ystart, ystop))

    # Check for errors in the slice indexes
    if (xstart < 0) or (ystart < 0) or \
       ((xstop-1) > ref_model.data.shape[-1]) or \
       ((ystop-1) > ref_model.data.shape[-2]):
        log.error("Science and reference file arrays not compatible")
        raise ValueError("Can't extract matching subarray from reference data")

    # Slice the reference model arrays
    sub_data = ref_model.data[ystart:ystop, xstart:xstop]
    sub_err = ref_model.err[ystart:ystop, xstart:xstop]
    sub_dq = ref_model.dq[ystart:ystop, xstart:xstop]

    # Create the sliced model
    sub_model = datamodels.SuperBiasModel(data=sub_data, err=sub_err, dq=sub_dq)

    # Return the sliced reference model
    return sub_model
