import numpy as np
from stdatamodels.jwst import datamodels

__all__ = [
    "generate_substripe_ranges",
    "generate_superstripe_ranges",
    "collate_superstripes",
    "clean_superstripe_metadata",
    "generate_stripe_int_times",
]


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
    fastaxis = np.abs(input_model.meta.subarray.fastaxis)

    slowsize = 2048
    slowdir = np.sign(input_model.meta.subarray.slowaxis)
    nreads1 = input_model.meta.subarray.multistripe_reads1

    # Generate slowaxis ranges to place stripes into parent frame
    stripe_ranges = generate_superstripe_ranges(input_model)
    srlist = []
    for key in stripe_ranges:
        if slowdir < 0:
            # If slowaxis is negative, need to read the stripes in reverse order
            # (but not reverse column order within a stripe, as they've been
            # transformed to science frame)
            srlist.append([slowsize - x for x in stripe_ranges[key][0][::-1]])
        else:
            srlist.append(list(stripe_ranges[key][0]))

    # Determine integration/stripe numbers
    nints, ngroups, ny, nx = input_model.data.shape
    n_sstr = input_model.meta.subarray.num_superstripe
    nints_sci = nints // n_sstr

    # Initialize new array shapes
    if fastaxis == 1:
        newdata = np.full((nints_sci, ngroups, slowsize, nx), np.nan, dtype=">f4")
        newgdq = np.full((nints_sci, ngroups, slowsize, nx), 0, dtype="uint8")
        newpdq = np.full(
            (slowsize, nx), datamodels.dqflags.pixel["REFERENCE_PIXEL"], dtype="uint32"
        )
    else:
        newdata = np.full((nints_sci, ngroups, ny, slowsize), np.nan, dtype=">f4")
        newgdq = np.full((nints_sci, ngroups, ny, slowsize), 0, dtype="uint8")
        newpdq = np.full(
            (ny, slowsize), datamodels.dqflags.pixel["REFERENCE_PIXEL"], dtype="uint32"
        )

    # Work through each set of stripes, pushing them into a common frame
    for integ in range(nints_sci):
        for stripe in range(n_sstr):
            """
            Long term, there may be cases where stripes overlap.
            This could be refactored to store each stripe in a separate
            plane, and overlaps could be reduced using a function of choice,
            e.g. np.median.
            """
            # Determine fast orient
            if fastaxis == 1:
                newslice = np.s_[integ, :, srlist[stripe][0] : srlist[stripe][1], :]
                # Determine end of slowread to drop refpix from
                if slowdir < 0:
                    dataslice = np.s_[integ * n_sstr + stripe, :, :-nreads1, :]
                else:
                    dataslice = np.s_[integ * n_sstr + stripe, :, nreads1:, :]
                newdata[newslice] = input_model.data[dataslice]
                newgdq[newslice] = input_model.groupdq[dataslice]
                if integ == 0:
                    newpdq[newslice[-2], newslice[-1]] = input_model.pixeldq[
                        stripe, dataslice[-2], dataslice[-1]
                    ]
            else:
                newslice = np.s_[integ, :, :, srlist[stripe][0] : srlist[stripe][1]]
                if slowdir < 0:
                    dataslice = np.s_[integ * n_sstr + stripe, :, :, :-nreads1]
                else:
                    dataslice = np.s_[integ * n_sstr + stripe, :, :, nreads1:]
                newdata[newslice] = input_model.data[dataslice]
                newgdq[newslice] = input_model.groupdq[dataslice]
                if integ == 0:
                    newpdq[newslice[-2], newslice[-1]] = input_model.pixeldq[
                        stripe, dataslice[-2], dataslice[-1]
                    ]

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
