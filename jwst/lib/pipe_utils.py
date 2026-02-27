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
        if not input_model.hasattr(extension):
            continue
        data = getattr(input_model, extension)
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
    if input_model.hasattr("dq"):
        do_not_use = (input_model.dq & dqflags.pixel["DO_NOT_USE"]).astype(bool)
        if input_model.dq.shape != data_shape:
            log.warning("Mismatched data shapes; skipping invalid data updates for extension 'dq'")
        else:
            is_invalid |= do_not_use

    # Update all the data extensions
    for extension in nan_extensions:
        if not input_model.hasattr(extension):
            continue
        data = getattr(input_model, extension)
        if data.shape != data_shape:
            continue
        data[is_invalid] = np.nan

    # Update the DQ extension
    if input_model.dq.shape == data_shape:
        input_model.dq[is_invalid] |= dqflags.pixel["DO_NOT_USE"]


def generate_substripe_ranges(sci_model):
    """
    TBD.

    Parameters
    ----------
    sci_model : JwstDataModel
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
    ysize_sci = sci_model.meta.subarray.ysize

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
    ranges[range_counter] = (linecount, linecount + nreads2)
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
    while sub_lines < ysize_sci:
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
        ranges[range_counter] = (linecount, linecount + nreads2)
        range_counter += 1
        linecount += nreads2
        sub_lines += nreads2

    if sub_lines != ysize_sci:
        raise ValueError(
            "Stripe readout resulted in mismatched reference array shape "
            "with respect to science array!"
        )

    return ranges
