import logging

import numpy as np
from stdatamodels.jwst import datamodels

from jwst.lib.reffile_utils import (
    detector_science_frame_transform,
    science_detector_frame_transform,
)

log = logging.getLogger(__name__)

__all__ = [
    "STRIPE_TO_STANDARD_SUBARRAY",
    "generate_substripe_ranges",
    "generate_superstripe_ranges",
    "stripe_read",
    "generate_stripe_reference",
    "collate_superstripes",
    "clean_superstripe_metadata",
    "generate_stripe_int_times",
    "stripe_read_times",
]

DETECTOR_FULL_SIZE = 2048

STRIPE_TO_STANDARD_SUBARRAY = {
    # NIRISS SOSS superstripes: reassemble to full size
    "SUB17STRIPE_SOSS": {"slow_size": DETECTOR_FULL_SIZE, "slow_start": 1},
    "SUB60STRIPE_SOSS": {"slow_size": DETECTOR_FULL_SIZE, "slow_start": 1},
    "SUB204STRIPE_SOSS": {"slow_size": DETECTOR_FULL_SIZE, "slow_start": 1},
    "SUB680STRIPE_SOSS": {"slow_size": DETECTOR_FULL_SIZE, "slow_start": 1},
}
"""
Map superstripe subarrays to standard subarray sizes.

Superstripe subarrays must be included here if the size along the slow axis
should be reassembled to include empty placeholder rows/columns.
"""


def _detector_coord_slow_start(slow_axis, slow_start):
    """
    Calculate the slow start value in detector coordinates for striped data.

    For striped data, the start value must be equal to the actual
    starting point for readouts: this is different from the usual definition,
    for which the starting point may be the start (slowaxis > 0) or the end
    (slowaxis < 0) of the readouts.

    If slow_start is 1 and the slow axis is negative, a warning is issued
    and 0 is returned, since it's assumed this is an error in the input
    data.

    Parameters
    ----------
    slow_axis : int
        The slow axis.
    slow_start : int
        The subarray start value in science coordinates, 1-indexed.

    Returns
    -------
    slow_start_idx : int
        The subarray start value in detector coordinates, 0-indexed.
    """
    if slow_axis < 0:
        # Handle some early products that might have the wrong value
        # for the slow start (i.e. start=1 when it should have been 2048).
        # Other incorrect cases are harder to detect and not handled here.
        if slow_start == 1:
            log.warning(
                f"Slow start is set to 1 for slowaxis={slow_axis}. "
                f"Setting detector start index to 0."
            )
            return 0

        # e.g. start = 2046 is pixel index 2
        slow_start_idx = DETECTOR_FULL_SIZE - slow_start
    else:
        # e.g. start = 3 is pixel index 2
        slow_start_idx = slow_start - 1

    return slow_start_idx


def generate_substripe_ranges(sci_model, science_frame=False):
    """
    Determine the bounds of each substripe based on the input model multistripe metadata.

    All ranges are indexed in detector orientation by default.

    Parameters
    ----------
    sci_model : `~stdatamodels.jwst.datamodels.JwstDataModel`
        The input datamodel with multistripe params defined.
    science_frame : bool, optional
        If True, the slow axis ranges are returned in science orientation
        instead of detector orientation.

    Returns
    -------
    dict
        Dictionary with keys "full", "subarray", "reference_full", and
        "reference_subarray". Each key has another dictionary as the value, with
        integer counter keys. Values are the ranges of array index in slowaxis
        corresponding to stripe shapes in the full array, striped subarray for the
        science and reference regions respectively.
    """
    nreads1 = sci_model.meta.subarray.multistripe_reads1
    nskips1 = sci_model.meta.subarray.multistripe_skips1
    nreads2 = sci_model.meta.subarray.multistripe_reads2
    nskips2 = sci_model.meta.subarray.multistripe_skips2
    repeat_stripe = sci_model.meta.subarray.repeat_stripe
    interleave_reads1 = sci_model.meta.subarray.interleave_reads1
    slowaxis = sci_model.meta.subarray.slowaxis
    subname = sci_model.meta.subarray.name

    # Get the slow axis size and start position
    if np.abs(slowaxis) > 1:
        slow_size = sci_model.meta.subarray.ysize
        slow_start = sci_model.meta.subarray.ystart
    else:
        slow_size = sci_model.meta.subarray.xsize
        slow_start = sci_model.meta.subarray.xstart

    # Start values are in science coord, 1-indexed: convert to detector start index
    slow_start = _detector_coord_slow_start(slowaxis, slow_start)
    if "DHS" in str(subname).upper() and slow_start != 0:
        # Accommodate early DHS data with incorrect slow start values
        log.warning("Setting detector start index to 0 for DHS")
        slow_start = 0

    # Check for valid input
    if nreads2 <= 0:
        raise ValueError(
            "Invalid value for multistripe_reads2 - stripe ranges could not be generated."
        )

    sub_ranges = {}
    full_ranges = {}
    reference_sub_ranges = {}
    reference_full_ranges = {}
    range_counter = 0
    reference_range_counter = 0

    # Track the read position in the full frame with linecount, and number of lines
    # read into subarray with sub_lines
    # Start at the subarray start
    linecount = slow_start
    sub_lines = 0

    # Make nreads1 row reads
    if nreads1 > 0:
        reference_full_ranges[reference_range_counter] = [linecount, linecount + nreads1]
        reference_sub_ranges[reference_range_counter] = [sub_lines, sub_lines + nreads1]
        reference_range_counter += 1
        linecount += nreads1
        sub_lines += nreads1

    # Now skip nskips1
    linecount += nskips1

    # Nreads2
    full_ranges[range_counter] = [linecount, linecount + nreads2]
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
    while sub_lines < slow_size:
        # If repeat_stripe, add interleaved rows to output and increment sub_lines
        if repeat_stripe > 0:
            linecount = slow_start

            if nreads1 > 0:
                reference_full_ranges[reference_range_counter] = [linecount, linecount + nreads1]
                reference_sub_ranges[reference_range_counter] = [sub_lines, sub_lines + nreads1]
                reference_range_counter += 1
                linecount += nreads1
                sub_lines += nreads1

            if interleave_reads1:
                interleave_skips += nskips2 + nreads2
                linecount += interleave_skips
            else:
                linecount += nskips1
        else:
            linecount += nskips2
        full_ranges[range_counter] = [linecount, linecount + nreads2]
        sub_ranges[range_counter] = [sub_lines, sub_lines + nreads2]
        range_counter += 1
        linecount += nreads2
        sub_lines += nreads2

    if sub_lines != slow_size:
        raise ValueError("Stripe readout does not match science array shape.")

    all_ranges = {
        "full": full_ranges,
        "subarray": sub_ranges,
        "reference_full": reference_full_ranges,
        "reference_subarray": reference_sub_ranges,
    }

    # Swap slow axis indices to science frame if needed
    if science_frame and slowaxis < 0:
        for range_set in all_ranges:
            for key, rge in all_ranges[range_set].items():
                if "full" in range_set:
                    all_ranges[range_set][key] = [DETECTOR_FULL_SIZE - x for x in rge[::-1]]
                else:
                    all_ranges[range_set][key] = [slow_size - x for x in rge[::-1]]

    return all_ranges


def generate_superstripe_ranges(sci_model, science_frame=False):
    """
    Return a dict of slowaxis ranges read into stripes.

    Given an input model with multistripe parameters in its metadata,
    return the slowaxis pixel ranges corresponding
    to the nreads2 positions read into each stripe of a set
    of superstripes.

    All ranges are indexed in detector orientation.

    Parameters
    ----------
    sci_model : `~stdatamodels.jwst.datamodels.JwstDataModel`
        The input datamodel with multistripe params defined.
    science_frame : bool, optional
        If True, the slow axis ranges are returned in science orientation
        instead of detector orientation.

    Returns
    -------
    dict
        Dictionary with keys "full", "subarray", "reference_full", and
        "reference_subarray". Each key has another dictionary as the value, with
        0-indexed stripe labels as keys. Values are the ranges of array index
        in slowaxis corresponding to stripe shapes in the full array, striped
        subarray for the science and reference regions respectively.
    """
    nreads1 = sci_model.meta.subarray.multistripe_reads1
    nskips1 = sci_model.meta.subarray.multistripe_skips1
    nreads2 = sci_model.meta.subarray.multistripe_reads2
    nskips2 = sci_model.meta.subarray.multistripe_skips2
    repeat_stripe = sci_model.meta.subarray.repeat_stripe
    interleave_reads1 = sci_model.meta.subarray.interleave_reads1
    superstripe_step = sci_model.meta.subarray.superstripe_step
    num_superstripe = sci_model.meta.subarray.num_superstripe
    fastaxis = sci_model.meta.subarray.fastaxis
    slowaxis = sci_model.meta.subarray.slowaxis
    subname = sci_model.meta.subarray.name

    # Get slow axis size and start position
    if np.abs(fastaxis) == 1:
        slow_size = sci_model.meta.subarray.ysize
        slow_start = sci_model.meta.subarray.ystart
    else:
        slow_size = sci_model.meta.subarray.xsize
        slow_start = sci_model.meta.subarray.xstart

    # Start values are in science coord, 1-indexed: convert to detector start index
    slow_start = _detector_coord_slow_start(slowaxis, slow_start)
    if "DHS" in str(subname).upper() and slow_start != 0:
        # Accommodate early DHS data with incorrect slow start values
        log.warning("Setting detector start index to 0 for DHS")
        slow_start = 0

    # Check for the number of pixels to read per stripe
    if nreads2 <= 0:
        # For pure superstripe, nreads2 is 0 and the pixels to
        # read is equal to the step size
        nreads2 = superstripe_step

    # Check for a pure substripe case: set nstripes to 1
    if num_superstripe == 0:
        num_superstripe = 1

    sub_ranges = {}
    full_ranges = {}
    reference_sub_ranges = {}
    reference_full_ranges = {}

    for stripe in range(num_superstripe):
        # Track the read position in the full frame with linecount, and number of lines
        # read into subarray with sub_lines
        # Start at the subarray start.
        linecount = slow_start
        sub_lines = 0

        # Make nreads1 row reads
        if nreads1 > 0:
            reference_full_ranges[stripe] = [(linecount, linecount + nreads1)]
            reference_sub_ranges[stripe] = [(sub_lines, sub_lines + nreads1)]
            linecount += nreads1
            sub_lines += nreads1
        else:
            reference_full_ranges[stripe] = []
            reference_sub_ranges[stripe] = []

        # Now skip nskips1 + superstripe_step * stripe, to get to the
        # beginning of the current stripe
        linecount += nskips1 + superstripe_step * stripe

        # Read Nreads2 for the stripe
        full_ranges[stripe] = [(linecount, linecount + nreads2)]
        sub_ranges[stripe] = [(sub_lines, sub_lines + nreads2)]
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
                linecount = slow_start

                if nreads1 > 0:
                    reference_full_ranges[stripe].append((linecount, linecount + nreads1))
                    reference_sub_ranges[stripe].append((sub_lines, sub_lines + nreads1))
                    linecount += nreads1
                    sub_lines += nreads1

                if interleave_reads1:
                    interleave_skips += nskips2 + nreads2
                    linecount += interleave_skips
                else:
                    linecount += nskips1 + superstripe_step * stripe
            else:
                linecount += nskips2
            full_ranges[stripe].append((linecount, linecount + nreads2))
            sub_ranges[stripe].append((sub_lines, sub_lines + nreads2))
            linecount += nreads2
            sub_lines += nreads2

        if sub_lines != slow_size:
            raise ValueError("Stripe readout does not match science array shape.")

    all_ranges = {
        "full": full_ranges,
        "subarray": sub_ranges,
        "reference_full": reference_full_ranges,
        "reference_subarray": reference_sub_ranges,
    }

    # Swap slow axis indices to science frame if needed
    if science_frame and slowaxis < 0:
        for range_set in all_ranges:
            for key, range_list in all_ranges[range_set].items():
                for i, rge in enumerate(range_list):
                    if "full" in range_set:
                        all_ranges[range_set][key][i] = tuple(
                            [DETECTOR_FULL_SIZE - x for x in rge[::-1]]
                        )
                    else:
                        all_ranges[range_set][key][i] = tuple([slow_size - x for x in rge[::-1]])
    return all_ranges


def stripe_read(sci_model, ref_model, attribs):
    """
    Generate a reference sub-model from science model multistripe params.

    Parameters
    ----------
    sci_model : `~stdatamodels.jwst.datamodels.JwstDataModel`
        Science model.
    ref_model : `~stdatamodels.jwst.datamodels.ReferenceFileModel`
        Reference model.
    attribs : list of str
        Attributes in the model to process.

    Returns
    -------
    sub_model : `~stdatamodels.jwst.datamodels.ReferenceFileModel`
        Generated reference sub-model.
    """
    # Get the science model multistripe params
    sci_meta = sci_model.meta

    # Get the reference model subarray params
    num_superstripe = getattr(sci_meta.subarray, "num_superstripe", 0)
    if num_superstripe > 0:
        sub_model = datamodels.ReferenceFileModel()
    else:
        sub_model = type(ref_model)()

    for attrib in attribs:
        ref_array = getattr(ref_model, attrib)

        # Apply subarray shape in fastaxis; slowaxis cutouts determined
        # in generate_stripe_reference
        # DHS reference files will likely be in FULL frame and will not have these subarray
        # values defined - they should pass over this if/elif block.
        if (
            np.abs(sci_meta.subarray.fastaxis) == 1
            and getattr(ref_model.meta.subarray, "xstart", None) is not None
        ):
            faststart_sci = sci_meta.subarray.xstart
            fastsize_sci = sci_meta.subarray.xsize

            # Get the reference model subarray params
            faststart_ref = ref_model.meta.subarray.xstart

            # Compute the slice indexes, in 0-indexed python frame
            faststart = faststart_sci - faststart_ref
            faststop = faststart + fastsize_sci
            ref_array = ref_array[..., faststart:faststop]
        elif (
            np.abs(sci_meta.subarray.fastaxis) == 2
            and getattr(ref_model.meta.subarray, "ystart", None) is not None
        ):
            faststart_sci = sci_meta.subarray.ystart
            fastsize_sci = sci_meta.subarray.ysize

            # Get the reference model subarray params
            faststart_ref = ref_model.meta.subarray.ystart

            # Compute the slice indexes, in 0-indexed python frame
            faststart = faststart_sci - faststart_ref
            faststop = faststart + fastsize_sci
            ref_array = ref_array[..., faststart:faststop, :]

        sub_model[attrib] = generate_stripe_reference(ref_array, sci_model)
    sub_model.update(ref_model)
    return sub_model


def generate_stripe_reference(ref_array, sci_model):
    """
    Generate stripe array from a reference array matching the full subarray size.

    Parameters
    ----------
    ref_array : ndarray
        The reference array to be sliced. Assumed to be full size along the slow axis.
    sci_model : `~stdatamodels.jwst.datamodels.JwstDataModel`
        The science datamodel.

    Returns
    -------
    stripe_out : ndarray
        Generated stripe array.
    """
    sci_meta = sci_model.meta

    # Extract science metadata
    num_superstripe = sci_meta.subarray.num_superstripe
    xsize_sci = sci_meta.subarray.xsize
    ysize_sci = sci_meta.subarray.ysize
    fastaxis = sci_meta.subarray.fastaxis
    slowaxis = sci_meta.subarray.slowaxis

    # Transform science data to detector frame
    ref_array = science_detector_frame_transform(ref_array, fastaxis, slowaxis)
    ref_shape = ref_array.shape
    if np.abs(fastaxis) == 1:
        slow_size = ysize_sci
        fast_size = xsize_sci
    else:
        slow_size = xsize_sci
        fast_size = ysize_sci

    if num_superstripe == 0:
        # SUBSTRIPE MODE
        stripe_out = np.zeros((*ref_shape[:-2], slow_size, fast_size), dtype=ref_array.dtype)
        all_ranges = generate_substripe_ranges(sci_model)
        for idx in all_ranges["full"]:
            f_reg = all_ranges["full"][idx]
            s_reg = all_ranges["subarray"][idx]
            stripe_out[..., s_reg[0] : s_reg[1], :] = ref_array[..., f_reg[0] : f_reg[1], :]
        for idx in all_ranges["reference_full"]:
            f_reg = all_ranges["reference_full"][idx]
            s_reg = all_ranges["reference_subarray"][idx]
            stripe_out[..., s_reg[0] : s_reg[1], :] = ref_array[..., f_reg[0] : f_reg[1], :]

    else:
        # SUPERSTRIPE MODE
        # All expected reference types for superstripe ramps are 2D
        if len(ref_shape) != 2:
            raise ValueError(f"Unsupported shape: len(ref_shape) == {len(ref_array.shape)}")

        # Rather than broadcast a 2-D array into many wasted dims, just provide the minimal
        # 3-D array, with one slice per superstripe in the third dimension. Expect steps
        # to handle the extra dimension.
        stripe_out = np.zeros((num_superstripe, slow_size, fast_size), dtype=ref_array.dtype)
        all_ranges = generate_superstripe_ranges(sci_model)
        for stripe in range(num_superstripe):
            # read in pixels over science areas
            full_range = all_ranges["full"]
            sub_range = all_ranges["subarray"]
            for s_reg, f_reg in zip(sub_range[stripe], full_range[stripe], strict=True):
                stripe_out[stripe, s_reg[0] : s_reg[1], :] = ref_array[f_reg[0] : f_reg[1], :]

            # read in reference pixels
            full_range = all_ranges["reference_full"]
            sub_range = all_ranges["reference_subarray"]
            for s_reg, f_reg in zip(sub_range[stripe], full_range[stripe], strict=True):
                stripe_out[stripe, s_reg[0] : s_reg[1], :] = ref_array[f_reg[0] : f_reg[1], :]

    # Transform from detector frame back to science frame
    stripe_out = detector_science_frame_transform(stripe_out, fastaxis, slowaxis)

    return stripe_out


def _slow_sub_from_full_range(slowaxis, full_range):
    """
    Get the subarray start and size along the slow axis from the derived ranges in the full frame.

    Parameters
    ----------
    slowaxis : int
        The slow axis.
    full_range : dict
        Keys are stripe numbers, values are lists of range tuples.

    Returns
    -------
    sci_slow_start : int
        The start value in science coordinates, 1-indexed, for recording
        in the output metadata.
    det_slow_start : int
        The start value in detector coordinates, 0-indexed, for determining
        output subarray indices.
    slow_size : int
        The size of the output subarray, to include all science values.
    """
    all_stripe_ranges = np.array(list(full_range.values()))
    if slowaxis < 0:
        # negative axis, convert to science coord
        sci_stripe_ranges = DETECTOR_FULL_SIZE - all_stripe_ranges
    else:
        sci_stripe_ranges = all_stripe_ranges
    sci_slow_start = np.min(sci_stripe_ranges) + 1
    det_slow_start = np.min(all_stripe_ranges)
    slow_size = np.max(all_stripe_ranges) - det_slow_start
    return sci_slow_start, det_slow_start, slow_size


def collate_superstripes(input_model):
    """
    Collate superstripes into arrays resembling the full detector/subarray shape.

    Parameters
    ----------
    input_model : `~stdatamodels.jwst.datamodels.SuperstripeRampModel`
        The datamodel containing a reference-pixel corrected ramp
        for superstripe data.

    Returns
    -------
    `~stdatamodels.jwst.datamodels.RampModel`
        The datamodel with superstripes collated into a single frame per set of stripes.
    """
    # First define the parent array shape
    fastaxis = input_model.meta.subarray.fastaxis
    slowaxis = input_model.meta.subarray.slowaxis
    if np.abs(fastaxis) == 1:
        fast_size = input_model.meta.subarray.xsize
    else:
        fast_size = input_model.meta.subarray.ysize

    # Determine integration/stripe numbers
    nints, ngroups, ny, nx = input_model.data.shape
    n_stripe = input_model.meta.subarray.num_superstripe
    if n_stripe == 0:
        n_stripe = 1
    nints_sci = nints // n_stripe

    # Generate slowaxis ranges to place science regions from stripes into parent frame
    all_ranges = generate_superstripe_ranges(input_model)
    full_range = all_ranges["full"]
    sub_range = all_ranges["subarray"]

    # Get the expected output size and start for the mode, if it needs to have row/column padding
    subarray_name = input_model.meta.subarray.name
    standard_subarray = STRIPE_TO_STANDARD_SUBARRAY.get(subarray_name)
    if standard_subarray is not None:
        slow_size = standard_subarray["slow_size"]
        sci_slow_start = standard_subarray["slow_start"]
        if slowaxis < 1:
            det_slow_start = DETECTOR_FULL_SIZE - (sci_slow_start + slow_size - 1)
        else:
            det_slow_start = sci_slow_start - 1
    else:
        # Otherwise, get the subarray start and extend from the derived ranges
        sci_slow_start, det_slow_start, slow_size = _slow_sub_from_full_range(slowaxis, full_range)
    n_repeat = max(len(r) for r in sub_range.values())
    ngroups_sci = ngroups * n_repeat

    # Initialize new array shapes in detector space
    newdata = np.full((nints_sci, ngroups_sci, slow_size, fast_size), np.nan, dtype="float32")
    newgdq = np.full((nints_sci, ngroups_sci, slow_size, fast_size), 0, dtype="uint8")
    newpdq = np.full(
        (slow_size, fast_size), datamodels.dqflags.pixel["REFERENCE_PIXEL"], dtype="uint32"
    )
    if input_model.meta.exposure.zero_frame:
        newzf = np.full((nints_sci, slow_size, fast_size), np.nan, dtype="float32")
    else:
        newzf = None

    # Transform old data to detector frame for indexing
    olddata = science_detector_frame_transform(input_model.data, fastaxis, slowaxis)
    if input_model.groupdq is None or input_model.groupdq.shape != (nints, ngroups, ny, nx):
        oldgdq = None
    else:
        oldgdq = science_detector_frame_transform(input_model.groupdq, fastaxis, slowaxis)
    if input_model.pixeldq is None or input_model.pixeldq.shape != (n_stripe, ny, nx):
        oldpdq = None
    else:
        oldpdq = science_detector_frame_transform(input_model.pixeldq, fastaxis, slowaxis)
    if (
        newzf is None
        or input_model.zeroframe is None
        or input_model.zeroframe.shape != (nints, ny, nx)
    ):
        oldzf = None
    else:
        oldzf = science_detector_frame_transform(input_model.zeroframe, fastaxis, slowaxis)

    # Work through each set of stripes, pushing them into a common frame
    # Long term, there may be cases where stripes overlap.
    # This could be refactored to store each stripe in a separate
    # plane, and overlaps could be reduced using a function of choice,
    # e.g. np.median.
    for integ in range(nints_sci):
        for group in range(ngroups_sci):
            repeat = group % n_repeat
            old_group_idx = group // n_repeat
            for stripe in range(n_stripe):
                stripe_idx = integ * n_stripe + stripe
                s_reg = sub_range[stripe][repeat]
                f_reg = [int(r - det_slow_start) for r in full_range[stripe][repeat]]

                # Propagate the science data
                newdata[integ, group, f_reg[0] : f_reg[1], :] = olddata[
                    stripe_idx, old_group_idx, s_reg[0] : s_reg[1], :
                ]

                # Propagate input groupdq only if appropriate
                if oldgdq is not None:
                    newgdq[integ, group, f_reg[0] : f_reg[1], :] = oldgdq[
                        stripe_idx, old_group_idx, s_reg[0] : s_reg[1], :
                    ]

                # Propagate pixeldq only for the first integration if appropriate
                if integ == 0:
                    if oldpdq is None:
                        # If not present, set a default DQ value in science regions
                        newpdq[f_reg[0] : f_reg[1], :] = 0
                    else:
                        newpdq[f_reg[0] : f_reg[1], :] = oldpdq[stripe, s_reg[0] : s_reg[1]]

                # Propagate zeroframe if present, but keep only the first repeat
                if oldzf is not None and old_group_idx == 0:
                    newzf[integ, f_reg[0] : f_reg[1], :] = oldzf[stripe_idx, s_reg[0] : s_reg[1], :]

    # Swap back to science frame
    newdata = detector_science_frame_transform(newdata, fastaxis, slowaxis)
    newgdq = detector_science_frame_transform(newgdq, fastaxis, slowaxis)
    newpdq = detector_science_frame_transform(newpdq, fastaxis, slowaxis)
    if newzf is not None:
        newzf = detector_science_frame_transform(newzf, fastaxis, slowaxis)

    new_model = datamodels.RampModel(
        data=newdata,
        groupdq=newgdq,
        pixeldq=newpdq,
        zeroframe=newzf,
        int_times=input_model.int_times,
    )
    new_model.update(input_model)

    new_model = generate_stripe_int_times(new_model)
    new_model = clean_superstripe_metadata(new_model, slow_start=sci_slow_start)

    # Add a list of readout times for the new groups if there were repeats
    # If not, the data is evenly sampled and explicit read times are not needed.
    if n_repeat > 1:
        new_model.meta.exposure.read_times = stripe_read_times(input_model)
    return new_model


def clean_superstripe_metadata(input_model, slow_start=1):
    """
    Update model metadata to match changes to arrays.

    Parameters
    ----------
    input_model : `~stdatamodels.jwst.datamodels.RampModel`
        The model with updated data array shapes matching a parent
        frame, e.g. full frame or a subarray, which consist of
        multiple superstripe integrations.
    slow_start : int
        The slowaxis subarray start position, 1-indexed to record.

    Returns
    -------
    `~stdatamodels.jwst.datamodels.RampModel`
        The model cleaned of metadata indicating the presence
        of superstripe data.
    """
    nstripe = input_model.meta.subarray.num_superstripe
    if nstripe == 0:
        nstripe = 1
    intstart = (
        input_model.meta.exposure.integration_start
        if input_model.meta.exposure.integration_start is not None
        else 1
    )
    intend = (
        input_model.meta.exposure.integration_end
        if input_model.meta.exposure.integration_end is not None
        else input_model.meta.exposure.nints
    )
    input_model.meta.exposure.integration_start = np.ceil(intstart / nstripe).astype(int)
    input_model.meta.exposure.integration_end = np.ceil(intend / nstripe).astype(int)
    input_model.meta.exposure.nints = np.ceil(input_model.meta.exposure.nints / nstripe).astype(int)
    input_model.meta.exposure.ngroups = input_model.data.shape[1]
    input_model.meta.subarray.multistripe_reads1 = None
    input_model.meta.subarray.multistripe_reads2 = None
    input_model.meta.subarray.multistripe_skips1 = None
    input_model.meta.subarray.multistripe_skips2 = None
    input_model.meta.subarray.repeat_stripe = None
    input_model.meta.subarray.interleave_reads1 = None
    input_model.meta.subarray.num_superstripe = None
    input_model.meta.subarray.superstripe_step = None
    input_model.meta.subarray.ysize, input_model.meta.subarray.xsize = input_model.data.shape[-2:]

    if np.abs(input_model.meta.subarray.slowaxis) == 1:
        input_model.meta.subarray.xstart = slow_start
    else:
        input_model.meta.subarray.ystart = slow_start

    return input_model


def generate_stripe_int_times(input_model):
    """
    Move input model INT_TIMES to stripe table, then condense ints to parent frame.

    Each output integration in the parent frame will have the start time from
    the first stripe readout and the end time from the last stripe readout included
    in the integration.  The mid time for the integration is assigned to the average
    of the start and end times.

    Parameters
    ----------
    input_model : `~stdatamodels.jwst.datamodels.RampModel`
        The model with updated data array shapes matching a parent
        frame, e.g. full frame or a subarray, which consist of
        multiple superstripe integrations.

    Returns
    -------
    `~stdatamodels.jwst.datamodels.RampModel`
        The model now with two INT_TIMES tables - one reflecting per-stripe
        information, preserved from input, and one reflecting times
        corresponding to the new parent frames.
    """
    if input_model.int_times is None or len(input_model.int_times) == 0:
        # No int_times, just return
        return input_model

    nstr = input_model.meta.subarray.num_superstripe
    if nstr == 0:
        # No superstripes, so int_times is okay as is
        return input_model

    nints_sci = len(input_model.int_times) // nstr
    default_intstart = (
        input_model.meta.exposure.integration_start
        if input_model.meta.exposure.integration_start is not None
        else 1
    )
    int_start = int(np.ceil(default_intstart / nstr))

    otab = np.array(
        list(
            zip(
                np.repeat(np.arange(len(input_model.int_times) // nstr) + int_start, nstr),
                np.arange(len(input_model.int_times)) % nstr + 1,
                [integ["int_start_MJD_UTC"] for integ in input_model.int_times],
                [integ["int_mid_MJD_UTC"] for integ in input_model.int_times],
                [integ["int_end_MJD_UTC"] for integ in input_model.int_times],
                [integ["int_start_BJD_TDB"] for integ in input_model.int_times],
                [integ["int_mid_BJD_TDB"] for integ in input_model.int_times],
                [integ["int_end_BJD_TDB"] for integ in input_model.int_times],
                strict=True,
            )
        ),
        dtype=input_model.get_dtype("int_times_stripe"),
    )

    input_model.int_times_stripe = otab

    cds_mjd_beg = np.full(nints_sci, np.nan, dtype="<f8")
    cds_mjd_mid = np.full(nints_sci, np.nan, dtype="<f8")
    cds_mjd_end = np.full(nints_sci, np.nan, dtype="<f8")
    cds_bjd_beg = np.full(nints_sci, np.nan, dtype="<f8")
    cds_bjd_mid = np.full(nints_sci, np.nan, dtype="<f8")
    cds_bjd_end = np.full(nints_sci, np.nan, dtype="<f8")

    input_inttimes = input_model.int_times
    nints_sci = len(input_inttimes) // nstr
    for i in range(nints_sci):
        cds_mjd_beg[i] = input_inttimes[i * nstr]["int_start_MJD_UTC"]
        cds_mjd_end[i] = input_inttimes[(i + 1) * nstr - 1]["int_end_MJD_UTC"]
        cds_mjd_mid[i] = (cds_mjd_end[i] - cds_mjd_beg[i]) / 2.0 + cds_mjd_beg[i]
        cds_bjd_beg[i] = input_inttimes[i * nstr]["int_start_BJD_TDB"]
        cds_bjd_end[i] = input_inttimes[(i + 1) * nstr - 1]["int_end_BJD_TDB"]
        cds_bjd_mid[i] = (cds_bjd_end[i] - cds_bjd_beg[i]) / 2.0 + cds_bjd_beg[i]

    otab2 = np.array(
        list(
            zip(
                np.arange(len(input_model.int_times) // nstr) + int_start,
                cds_mjd_beg,
                cds_mjd_mid,
                cds_mjd_end,
                cds_bjd_beg,
                cds_bjd_mid,
                cds_bjd_end,
                strict=True,
            )
        ),
        dtype=input_model.int_times.dtype,
    )

    input_model.int_times = otab2

    return input_model


def _frametime(ncols, nrows, noutputs, end_of_row_pad, extra_rows, extra_clocks):
    """
    Calculate a frametime from the number of pixels read and other detector characteristics.

    Parameters
    ----------
    ncols : int
        Number of columns read.
    nrows : int
        Number of rows read.
    noutputs : int
        Number of outputs.
    end_of_row_pad : int
        Extra pad for the end of row.
    extra_rows : int
        Extra row pad for reset timing.
    extra_clocks : int
        Extra clock pad for reset timing.

    Returns
    -------
    frametime : float
        The total frametime in seconds.
    """
    frametime = (
        ((ncols / noutputs) + end_of_row_pad) * (nrows + extra_rows) + extra_clocks
    ) * 10.0e-6  # in seconds

    return frametime


def stripe_read_times(stripe_ramp):
    """
    Calculate read times for each readout, including subframe samples.

    Parameters
    ----------
    stripe_ramp : `~stdatamodels.jwst.datamodels.SuperstripeRampModel` or \
                  `~stdatamodels.jwst.datamodels.RampModel`
        The input ramp model containing multistripe data.

    Returns
    -------
    read_times : ndarray
        Array of read times with length ``ngroup * nrepeat``, where ``ngroup``
        is the number of groups taken per integration, and ``nrepeat`` is
        the number of stripe repeats contained within a ramp image.  Each element
        contains a list with the read time for each frame averaged into the group
        readout (usually 1 for striped modes).
    """
    nframes = stripe_ramp.meta.exposure.nframes
    groupgap = stripe_ramp.meta.exposure.groupgap
    ngroups = stripe_ramp.meta.exposure.ngroups
    noutputs = stripe_ramp.meta.exposure.noutputs
    nreads1 = stripe_ramp.meta.subarray.multistripe_reads1
    nreads2 = stripe_ramp.meta.subarray.multistripe_reads2
    if np.abs(stripe_ramp.meta.subarray.slowaxis) == 1:
        nrows = stripe_ramp.meta.subarray.xsize
        ncols = stripe_ramp.meta.subarray.ysize
    else:
        nrows = stripe_ramp.meta.subarray.ysize
        ncols = stripe_ramp.meta.subarray.xsize

    nreads1 = 0 if nreads1 is None else nreads1
    nreads2 = 0 if nreads2 is None else nreads2

    # Get extra padding values
    end_of_row_pad = 12
    extra_clocks = 0
    if noutputs == 1 and ncols >= 9:
        extra_rows = 2
    else:
        extra_rows = 3
    if noutputs == 4:
        # for full-frame and stripe mode
        extra_rows = 1
        extra_clocks = 1

    frametime = _frametime(ncols, nrows, noutputs, end_of_row_pad, extra_rows, extra_clocks)

    # Reset frame time: same as frametime, except only the subframe is reset.
    # Each "subframe" has reads1+reads2 rows in it
    subframe_nrows = nreads1 + nreads2
    reset_frametime = _frametime(
        ncols, subframe_nrows, noutputs, end_of_row_pad, extra_rows, extra_clocks
    )

    # "subframe" time: the time between successive repeats of reads1+reads2 after
    # each fsync. This is the time between UTR samples within the frame.
    # There are no extra clocks or extra rows in the subframes.
    extra_clocks = 0
    extra_rows = 0
    subframetime = _frametime(
        ncols, subframe_nrows, noutputs, end_of_row_pad, extra_rows, extra_clocks
    )

    log.debug(f"Frame time is: {frametime} s")
    log.debug(f"Reset time is: {reset_frametime} s")
    log.debug(f"Subframe time is: {subframetime} s")

    if subframe_nrows > 0:
        # Get the repeats per frame
        n_repeat = nrows // subframe_nrows

        # Delta t for each group is: [reset_frametime, subframe_time, subframe_time, ...]
        delta_t = [reset_frametime] + [subframetime] * (n_repeat - 1)
    else:
        # Without repeats, the delta t is just the frame time
        delta_t = [frametime]

    # Assemble the list of read times for all frames
    read_times = []
    total_time = 0
    for _ in range(ngroups):
        # Loop over subframes
        for dt in delta_t:
            # One read at this dt for each frame averaged in the group
            frame_reads = [total_time + dt + frame * frametime for frame in range(nframes)]
            read_times.append(frame_reads)
            # Track the total time elapsed
            total_time += dt

        # Add in the total times for any averaged frames within a group or
        # any dropped frames between groups
        total_time += frametime * (nframes - 1 + groupgap)

    return read_times
