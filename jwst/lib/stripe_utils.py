import numpy as np
from stdatamodels.jwst import datamodels

from jwst.lib.reffile_utils import (
    detector_science_frame_transform,
    science_detector_frame_transform,
)

__all__ = [
    "generate_substripe_ranges",
    "generate_superstripe_ranges",
    "stripe_read",
    "generate_stripe_reference",
    "collate_superstripes",
    "clean_superstripe_metadata",
    "generate_stripe_int_times",
]


def generate_substripe_ranges(sci_model):
    """
    Determine the bounds of each substripe based on the input model multistripe metadata.

    All ranges are indexed in detector orientation.

    Parameters
    ----------
    sci_model : `~stdatamodels.jwst.datamodels.JwstDataModel`
        The input datamodel with multistripe params defined.

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
    if np.abs(slowaxis) > 1:
        slowsize = sci_model.meta.subarray.ysize
    else:
        slowsize = sci_model.meta.subarray.xsize

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
    linecount = 0
    sub_lines = 0

    # Start at 0, make nreads1 row reads
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
    while sub_lines < slowsize:
        # If repeat_stripe, add interleaved rows to output and increment sub_lines
        if repeat_stripe > 0:
            linecount = 0
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

    if sub_lines != slowsize:
        raise ValueError("Stripe readout does not match science array shape.")

    return {
        "full": full_ranges,
        "subarray": sub_ranges,
        "reference_full": reference_full_ranges,
        "reference_subarray": reference_sub_ranges,
    }


def generate_superstripe_ranges(sci_model):
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
    xsize_sci = sci_model.meta.subarray.xsize
    ysize_sci = sci_model.meta.subarray.ysize
    fastaxis = sci_model.meta.subarray.fastaxis
    if np.abs(fastaxis) == 1:
        slow_size = ysize_sci
    else:
        slow_size = xsize_sci

    # Check for valid input
    if nreads2 <= 0:
        raise ValueError(
            "Invalid value for multistripe_reads2 - stripe ranges could not be generated."
        )

    sub_ranges = {}
    full_ranges = {}
    reference_sub_ranges = {}
    reference_full_ranges = {}
    for stripe in range(num_superstripe):
        # Track the read position in the full frame with linecount, and number of lines
        # read into subarray with sub_lines
        linecount = 0
        sub_lines = 0

        # Start at 0, make nreads1 row reads for reference pixels
        reference_full_ranges[stripe] = [(linecount, linecount + nreads1)]
        reference_sub_ranges[stripe] = [(sub_lines, sub_lines + nreads1)]
        linecount += nreads1
        sub_lines += nreads1

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
                linecount = 0

                reference_full_ranges[stripe].append((linecount, linecount + nreads1))
                reference_sub_ranges[stripe].append((sub_lines, sub_lines + nreads1))
                linecount += nreads1
                sub_lines += nreads1

                if interleave_reads1:
                    interleave_skips += nskips2 + nreads2
                    linecount += interleave_skips
                else:
                    linecount += nskips1
            else:
                linecount += nskips2
            full_ranges[stripe].append((linecount, linecount + nreads2))
            sub_ranges[stripe].append((sub_lines, sub_lines + nreads2))
            linecount += nreads2
            sub_lines += nreads2

        if sub_lines != slow_size:
            raise ValueError("Stripe readout does not match science array shape.")

    return {
        "full": full_ranges,
        "subarray": sub_ranges,
        "reference_full": reference_full_ranges,
        "reference_subarray": reference_sub_ranges,
    }


def stripe_read(sci_model, ref_model, attribs):
    """
    Generate sub-model from science model multistripe params.

    Parameters
    ----------
    sci_model, ref_model : DataModel
        Science and reference models, respectively.

    attribs : list of str
        Attributes in the model to process.

    Returns
    -------
    sub_model : DataModel
        Generated sub-model.
    """
    # Get the science model multistripe params
    sci_meta = sci_model.meta

    num_superstripe = getattr(sci_meta.subarray, "num_superstripe", 0)
    # We need to extract the number of science integrations in this file, which is tangled up with
    #  - the number of integrations in the exposure
    #  - the number of superstripes each science integration may be divided between
    # Substripe data may have 0 for num_superstripe, so we use 1 as the integration divisor
    int_divisor = max(num_superstripe, 1)
    sci_nints = sci_model.data.shape[0] // int_divisor

    # Get the reference model subarray params
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
            faststart_sci = sci_model.meta.subarray.xstart
            fastsize_sci = sci_meta.subarray.xsize

            # Get the reference model subarray params
            faststart_ref = ref_model.meta.subarray.xstart

            # Compute the slice indexes, in 0-indexed python frame
            faststart = faststart_sci - faststart_ref
            faststop = faststart + fastsize_sci
            ref_array = ref_array[..., faststart:faststop]
        elif getattr(ref_model.meta.subarray, "ystart", None) is not None:
            faststart_sci = sci_model.meta.subarray.ystart
            fastsize_sci = sci_meta.subarray.ysize

            # Get the reference model subarray params
            faststart_ref = ref_model.meta.subarray.ystart

            # Compute the slice indexes, in 0-indexed python frame
            faststart = faststart_sci - faststart_ref
            faststop = faststart + fastsize_sci
            ref_array = ref_array[..., faststart:faststop, :]

        sub_model[attrib] = generate_stripe_reference(ref_array, sci_model, sci_nints)
    sub_model.update(ref_model)
    return sub_model


def generate_stripe_reference(ref_array, sci_model, sci_nints):
    """
    Generate stripe array from a reference array matching the full subarray size.

    Parameters
    ----------
    ref_array : np.array
        The scene to be sliced.
    sci_model : `~stdatamodels.datamodels.JwstDataModel`
        The science datamodel metadata tree.
    sci_nints : int
        The number of science integrations in the science datamodel. Not equivalent to nints when
        the science exposure is segmented.

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
    ngroups = sci_meta.exposure.ngroups

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
        # First alter subarray shape to broadcast stripe size to fill "full" subarray
        ref_native_dims = len(ref_shape)
        if len(ref_shape) == 2:
            ref_array = ref_array[np.newaxis, np.newaxis, :].repeat(ngroups, axis=1)
            ref_shape = ref_array.shape
        if len(ref_shape) == 4:
            ref_nints, _, ysize, xsize = ref_shape
        else:
            raise ValueError(f"Unsupported shape: len(ref_shape) == {len(ref_array.shape)}")

        # If reference file has more integrations than science, only process through
        # necessary integrations.
        # If science has more integrations, we'll broadcast the reference arrays to
        # match the number of arrays in the science frame.
        nints = min(ref_nints, sci_nints)
        stripe_out = np.zeros(
            (nints * num_superstripe, ngroups, slow_size, fast_size), dtype=ref_array.dtype
        )
        all_ranges = generate_superstripe_ranges(sci_model)
        for integ in range(nints):
            for stripe in range(num_superstripe):
                stripe_idx = integ * num_superstripe + stripe
                full_range = all_ranges["full"]
                sub_range = all_ranges["subarray"]
                for s_reg, f_reg in zip(sub_range[stripe], full_range[stripe], strict=True):
                    stripe_out[stripe_idx, :, s_reg[0] : s_reg[1], :] = ref_array[
                        integ, :, f_reg[0] : f_reg[1], :
                    ]
                full_range = all_ranges["reference_full"]
                sub_range = all_ranges["reference_subarray"]
                for s_reg, f_reg in zip(sub_range[stripe], full_range[stripe], strict=True):
                    stripe_out[stripe_idx, :, s_reg[0] : s_reg[1], :] = ref_array[
                        integ, :, f_reg[0] : f_reg[1], :
                    ]

        # If multistripe but ref array is typically 2-D:
        # rather than broadcast a 2-D array into many wasted dims, just provide the minimal
        # 3-D array, with one slice per superstripe in the third dimension. Expect steps
        # to handle the extra dimension.
        if ref_native_dims == 2:
            stripe_out = stripe_out[:, 0, :, :]
        # If multistripe and output currently has one plane per stripe only,
        # and reference file is expected to have 4 dimensions,
        # broadcast arrays into sci_nints copies so that direct application
        # of the reference array to the science array is possible
        elif stripe_out.shape[0] == num_superstripe:
            stripe_out = np.tile(stripe_out, reps=(sci_nints, 1, 1, 1))

    # Transform from detector frame back to science frame
    stripe_out = detector_science_frame_transform(stripe_out, fastaxis, slowaxis)

    return stripe_out


def collate_superstripes(input_model):
    """
    Collate superstripes into arrays resembling the full detector/subarray shape.

    Parameters
    ----------
    input_model : `~stdatamodels.jwst.datamodels.RampModel`
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
    slow_size = 2048
    if np.abs(fastaxis) == 1:
        fast_size = input_model.meta.subarray.xsize
    else:
        fast_size = input_model.meta.subarray.ysize

    # Determine integration/stripe numbers
    nints, ngroups, ny, nx = input_model.data.shape
    n_sstr = input_model.meta.subarray.num_superstripe
    nints_sci = nints // n_sstr

    # Generate slowaxis ranges to place stripes into parent frame
    all_ranges = generate_superstripe_ranges(input_model)
    full_range = all_ranges["full"]
    sub_range = all_ranges["subarray"]

    # Initialize new array shapes in detector space
    newdata = np.full((nints_sci, ngroups, slow_size, fast_size), np.nan, dtype="float32")
    newgdq = np.full((nints_sci, ngroups, slow_size, fast_size), 0, dtype="uint8")
    newpdq = np.full(
        (slow_size, fast_size), datamodels.dqflags.pixel["REFERENCE_PIXEL"], dtype="uint32"
    )

    # Transform old data to detector frame for indexing
    olddata = science_detector_frame_transform(input_model.data, fastaxis, slowaxis)
    oldgdq = science_detector_frame_transform(input_model.groupdq, fastaxis, slowaxis)
    oldpdq = science_detector_frame_transform(input_model.pixeldq, fastaxis, slowaxis)

    # Work through each set of stripes, pushing them into a common frame
    # Long term, there may be cases where stripes overlap.
    # This could be refactored to store each stripe in a separate
    # plane, and overlaps could be reduced using a function of choice,
    # e.g. np.median.
    for integ in range(nints_sci):
        for stripe in range(n_sstr):
            stripe_idx = integ * n_sstr + stripe
            for s_reg, f_reg in zip(sub_range[stripe], full_range[stripe], strict=True):
                newdata[integ, :, f_reg[0] : f_reg[1], :] = olddata[
                    stripe_idx, :, s_reg[0] : s_reg[1], :
                ]
                newgdq[integ, :, f_reg[0] : f_reg[1], :] = oldgdq[
                    stripe_idx, :, s_reg[0] : s_reg[1], :
                ]
                if integ == 0:
                    newpdq[f_reg[0] : f_reg[1], :] = oldpdq[stripe, s_reg[0] : s_reg[1]]

    # Swap back to science frame
    newdata = detector_science_frame_transform(newdata, fastaxis, slowaxis)
    newgdq = detector_science_frame_transform(newgdq, fastaxis, slowaxis)
    newpdq = detector_science_frame_transform(newpdq, fastaxis, slowaxis)

    new_model = datamodels.RampModel(
        data=newdata,
        groupdq=newgdq,
        pixeldq=newpdq,
        int_times=input_model.int_times,
    )
    new_model.update(input_model)

    new_model = generate_stripe_int_times(new_model)
    new_model = clean_superstripe_metadata(new_model)

    return new_model


def clean_superstripe_metadata(input_model):
    """
    Update model metadata to match changes to arrays.

    Parameters
    ----------
    input_model : `~stdatamodels.jwst.datamodels.RampModel`
        The model with updated data array shapes matching a parent
        frame, e.g. full frame or a subarray, which consist of
        multiple superstripe integrations.

    Returns
    -------
    `~stdatamodels.jwst.datamodels.RampModel`
        The model cleaned of metadata indicating the presence
        of superstripe data.
    """
    input_model.meta.exposure.integration_start = np.ceil(
        input_model.meta.exposure.integration_start / input_model.meta.subarray.num_superstripe
    ).astype(int)
    input_model.meta.exposure.integration_end = np.ceil(
        input_model.meta.exposure.integration_end / input_model.meta.subarray.num_superstripe
    ).astype(int)
    input_model.meta.exposure.nints = np.ceil(
        input_model.meta.exposure.nints / input_model.meta.subarray.num_superstripe
    ).astype(int)

    input_model.meta.subarray.multistripe_reads1 = None
    input_model.meta.subarray.multistripe_reads2 = None
    input_model.meta.subarray.multistripe_skips1 = None
    input_model.meta.subarray.multistripe_skips2 = None
    input_model.meta.subarray.repeat_stripe = None
    input_model.meta.subarray.interleave_reads1 = None
    input_model.meta.subarray.num_superstripe = None
    input_model.meta.subarray.superstripe_step = None
    input_model.meta.subarray.ysize, input_model.meta.subarray.xsize = input_model.data.shape[-2:]

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
        return input_model

    nstr = input_model.meta.subarray.num_superstripe
    nints_sci = len(input_model.int_times) // nstr

    otab = np.array(
        list(
            zip(
                np.repeat(np.arange(len(input_model.int_times) // nstr) + 1, nstr),
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
                np.arange(len(input_model.int_times) // nstr) + 1,
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
