"""Pipeline utilities objects."""

import logging

import numpy as np
from stdatamodels.jwst.datamodels import JwstDataModel, dqflags
from stdatamodels.properties import ObjectNode

from jwst.associations.lib.dms_base import TSO_EXP_TYPES

log = logging.getLogger(__name__)

__all__ = ["is_tso", "is_irs2", "match_nans_and_flags"]


def is_tso(model):
    """
    Check if data is a Time Series Observation (TSO).

    Parameters
    ----------
    model : `~stdatamodels.jwst.datamodels.JwstDataModel`
        Data to check.

    Returns
    -------
    is_tso : bool
       `True` if the model represents TSO data.
    """
    is_tso = False

    # Check on JWST-specific TSOVISIT flag
    try:
        is_tso = model.meta.visit.tsovisit
    except AttributeError:
        pass

    # Check on exposure types
    try:
        is_tso = is_tso or model.meta.exposure.type.lower() in TSO_EXP_TYPES
    except AttributeError:
        pass

    # Check on number of integrations
    try:
        if model.meta.exposure.nints is not None and model.meta.exposure.nints < 2:
            is_tso = False
    except AttributeError:
        pass

    # We've checked everything.
    return is_tso


def is_irs2(model):
    """
    Check whether the data are in IRS2 format.

    This currently assumes that only full-frame, near-infrared data can be
    taken using the IRS2 readout pattern.

    Parameters
    ----------
    model : `~stdatamodels.jwst.datamodels.JwstDataModel` or ndarray
        Data to check.

    Returns
    -------
    status : bool
       `True` if the data are in IRS2 format.
    """
    if isinstance(model, np.ndarray):
        shape = model.shape
    else:
        try:
            shape = model.data.shape
        except AttributeError:
            return False

    max_length = 2048

    irs2_axis_length = max(shape[-1], shape[-2])

    if irs2_axis_length > max_length:
        return True
    else:
        return False


def match_nans_and_flags(input_model):
    """
    Ensure data, error, variance, and DQ are marked consistently for invalid data.

    Invalid data is assumed to be any pixel set to NaN in any one of the
    data, error, or variance arrays, or else set to the DO_NOT_USE flag
    in the DQ array.

    The input model is updated in place with NaNs or DO_NOT_USE flags, as
    appropriate, at all invalid data locations.

    Parameters
    ----------
    input_model : DataModel
        Input model containing some combination of data, dq, err, var_rnoise,
        var_poisson, and var_flat extensions. These extensions must all have
        matching dimensions if present.
    """
    # Check for datamodel input or slit instance
    if not isinstance(input_model, JwstDataModel) and not isinstance(input_model, ObjectNode):
        raise TypeError(f"Input {type(input_model)} is not a datamodel.")

    # Build up the invalid data flags from each available data extension.
    is_invalid = None
    data_shape = None
    nan_extensions = ["data", "err", "var_rnoise", "var_poisson", "var_flat"]
    for extension in nan_extensions:
        data = getattr(input_model, extension, None)
        if data is None:
            continue
        if is_invalid is None:
            is_invalid = np.isnan(data)
            data_shape = data.shape
        else:
            if data.shape != data_shape:
                log.warning(
                    "Mismatched data shapes; skipping invalid data updates for extension '%s'",
                    extension,
                )
                continue
            is_invalid |= np.isnan(data)

    # Nothing to do if no extensions were found to update
    if is_invalid is None:
        return

    # Add in invalid flags from the DQ extension if present
    if getattr(input_model, "dq", None) is not None:
        do_not_use = (input_model.dq & dqflags.pixel["DO_NOT_USE"]).astype(bool)
        if input_model.dq.shape != data_shape:
            log.warning("Mismatched data shapes; skipping invalid data updates for extension 'dq'")
        else:
            is_invalid |= do_not_use
            input_model.dq[is_invalid] |= dqflags.pixel["DO_NOT_USE"]

    # Update all the data extensions
    for extension in nan_extensions:
        data = getattr(input_model, extension, None)
        if data is None:
            continue
        if data.shape != data_shape:
            continue
        data[is_invalid] = np.nan


def generate_substripe_ranges(sci_model, subarray_ranges=False):
    """
    Determine the bounds of each substripe based on the input model multistripe metadata.

    Parameters
    ----------
    sci_model : `~stdatamodels.jwst.datamodels.JwstDataModel`
        The input datamodel with multistripe params defined.
    subarray_ranges : bool
        If true, dict containing ranges in subarray frame will also be provided.

    Returns
    -------
    dict or tuple(dict)
        Dictionary with keys as int counter, values as ranges of array index in slowaxis
        corresponding to stripe shapes. If subarray_ranges was True, a second dictionary
        will be provided with equivalent ranges in the subarray frame.
    """
    nreads1 = sci_model.meta.subarray.multistripe_reads1
    nskips1 = sci_model.meta.subarray.multistripe_skips1
    nreads2 = sci_model.meta.subarray.multistripe_reads2
    nskips2 = sci_model.meta.subarray.multistripe_skips2
    repeat_stripe = sci_model.meta.subarray.repeat_stripe
    interleave_reads1 = sci_model.meta.subarray.interleave_reads1
    slowaxis = sci_model.meta.subarray.slowaxis
    if np.abs(slowaxis) > 1:
        slowsize = sci_model.meta.subarray.ysize
    else:
        slowsize = sci_model.meta.subarray.xsize

    sub_ranges = {}
    ranges = {}
    range_counter = 0
    # Track the read position in the full frame with linecount, and number of lines
    # read into subarray with sub_lines
    linecount = 0
    sub_lines = 0

    # Start at 0, make nreads1 row reads
    linecount += nreads1
    sub_lines += nreads1
    # Now skip nskips1
    linecount += nskips1
    # Nreads2
    ranges[range_counter] = [linecount, linecount + nreads2]
    sub_ranges[range_counter] = [sub_lines, sub_lines + nreads2]
    range_counter += 1
    linecount += nreads2
    sub_lines += nreads2

    # Now, while the output size is less than the science array size:
    # 1a. If repeat_stripe, reset linecount (HEAD) to initial position
    #     after every nreads2.
    # 1b. Else, do nskips2 followed by nreads2 until subarray complete.
    # 2.  Following 1a., repeat sequence of nreads1, skips*, nreads2
    #     until complete. For skips*:
    # 3a. If interleave_reads1, value of skips increments by nreads2 +
    #     nskips2 for each stripe read.
    # 3b. If not interleave, each loop after linecount reset is simply
    #     nreads1 + nskips1 + nreads2.
    interleave_skips = nskips1
    if nreads2 <= 0:
        raise ValueError(
            "Invalid value for multistripe_reads2 - "
            "cutout for reference file could not be "
            "generated!"
        )
    while sub_lines < slowsize:
        # If repeat_stripe, add interleaved rows to output and increment sub_lines
        if repeat_stripe > 0:
            linecount = 0
            linecount += nreads1
            sub_lines += nreads1
            if interleave_reads1:
                interleave_skips += nskips2 + nreads2
                linecount += interleave_skips
            else:
                linecount += nskips1
        else:
            linecount += nskips2
        ranges[range_counter] = [linecount, linecount + nreads2]
        sub_ranges[range_counter] = [sub_lines, sub_lines + nreads2]
        range_counter += 1
        linecount += nreads2
        sub_lines += nreads2

    if sub_lines != slowsize:
        raise ValueError(
            "Stripe readout resulted in mismatched reference array shape "
            "with respect to science array!"
        )

    if slowaxis < 0:
        for rge in ranges.keys():
            for i, row in enumerate(ranges[rge]):
                ranges[rge][i] = 2048 - row
                sub_ranges[rge][i] = 2048 - row

    if subarray_ranges:
        return ranges, sub_ranges
    else:
        return ranges


def generate_superstripe_ranges(sci_model):
    """
    Return a dict of slowaxis ranges read into stripes.

    Given an input model with multistripe parameters in its metadata,
    return the slowaxis pixel ranges corresponding
    to the nreads2 positions read into each stripe of a set
    of superstripes. Dictionary keys are 0-indexed stripe
    labels.

    Parameters
    ----------
    sci_model : `~stdatamodels.jwst.datamodels.JwstDataModel`
        The input datamodel with multistripe params defined.

    Returns
    -------
    dict
        Dictionary with keys as int counter, values as ranges of array index in slowaxis
        corresponding to stripe shapes.
    """
    nreads1 = sci_model.meta.subarray.multistripe_reads1
    nskips1 = sci_model.meta.subarray.multistripe_skips1
    nreads2 = sci_model.meta.subarray.multistripe_reads2
    nskips2 = sci_model.meta.subarray.multistripe_skips2
    repeat_stripe = sci_model.meta.subarray.repeat_stripe
    interleave_reads1 = sci_model.meta.subarray.interleave_reads1
    superstripe_step = sci_model.meta.subarray.superstripe_step
    num_superstripe = sci_model.meta.subarray.num_superstripe
    xsize_sci = sci_model.meta.subarray.xsize
    ysize_sci = sci_model.meta.subarray.ysize
    fastaxis = sci_model.meta.subarray.fastaxis

    if np.abs(fastaxis) == 1:
        slow_size = ysize_sci
    else:
        slow_size = xsize_sci

    ranges = {}
    range_counter = 0

    for stripe in range(num_superstripe):
        # Track the read position in the full frame with linecount, and number of lines
        # read into subarray with sub_lines
        linecount = 0
        sub_lines = 0

        # Start at 0, make nreads1 row reads
        linecount += nreads1
        sub_lines += nreads1
        # Now skip nskips1 + superstripe_step * stripe
        linecount += nskips1 + superstripe_step * stripe
        # Nreads2
        ranges[range_counter] = [(linecount, linecount + nreads2)]
        linecount += nreads2
        sub_lines += nreads2

        # Now, while the output size is less than the science array size:
        # 1a. If repeat_stripe, reset linecount (HEAD) to initial position
        #     after every nreads2.
        # 1b. Else, do nskips2 followed by nreads2 until subarray complete.
        # 2.  Following 1a., repeat sequence of nreads1, skips*, nreads2
        #     until complete. For skips*:
        # 3a. If interleave_reads1, value of skips increments by nreads2 +
        #     nskips2 for each stripe read.
        # 3b. If not interleave, each loop after linecount reset is simply
        #     nreads1 + nskips1 + nreads2.
        interleave_skips = nskips1

        while sub_lines < slow_size:
            # If repeat_stripe, add interleaved rows to output and increment sub_lines
            if repeat_stripe > 0:
                linecount = 0

                linecount += nreads1
                sub_lines += nreads1

                if interleave_reads1:
                    interleave_skips += nskips2 + nreads2
                    linecount += interleave_skips
                else:
                    linecount += nskips1
            else:
                linecount += nskips2
            ranges[range_counter].append((linecount, linecount + nreads2))
            linecount += nreads2
            sub_lines += nreads2
        range_counter += 1

        if sub_lines != slow_size:
            raise ValueError(
                "Stripe readout resulted in mismatched reference array shape "
                "with respect to science array!"
            )
    return ranges
