import logging

import numpy as np
from astropy.stats import sigma_clipped_stats
from scipy.ndimage import convolve1d

from stdatamodels.jwst.datamodels import dqflags

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def correct_model(
    output_model, irs2_model, scipix_n_default=16, refpix_r_default=4, pad=8, preserve_refpix=False
):
    """
    Correct an input NIRSpec IRS2 datamodel using reference pixels.

    Parameters
    ----------
    output_model : ramp model
        The input science data model.

    irs2_model : IRS2 model
        The reference file model for IRS2 correction.

    scipix_n_default : int
        Number of regular samples before stepping out to collect
        reference samples.

    refpix_r_default : int
        Number of reference samples before stepping back in to collect
        regular samples.

    pad : int
        The effective number of pixels sampled during the pause at the end
        of each row (new-row overhead).  The padding is needed to preserve
        the phase of temporally periodic signals.

    preserve_refpix : bool
        If True, reference pixels will be preserved in the output.
        This is not used in the science pipeline, but is necessary to
        create new bias files for IRS2 mode.

    Returns
    -------
    output_model : ramp model
        The science data with reference output and reference pixels
        subtracted.
    """
    #
    """Readout parameters
    scipix_n   16         Number of regular samples before stepping out
                          to collect reference samples
    refpix_r    4         Number of reference samples before stepping back
                          in to collect regular samples
    NFOH        1 row     New frame overhead (714?)
    NROH        8 pixels  New row overhead (`pad`)
    JOH         1 pixel   Jump overhead for stepping out to or in from
                          reference pixels
    TPIX       10 microseconds  Pixel dwell time
    tframe     14.5889 s  Frame readout time

    The image and reference data will be rearranged into a 1-D array
    containing the values in time order, i.e. the element number * TPIX
    is the relative time at which a pixel was read out.  This array will
    have elements corresponding to the gaps (overheads) when no pixel was
    being read.  This 1-D array has length 1,458,176, which is equal to
    712 * 2048:

    ((scipix_n + refpix_r + 2) * (512 // scipix_n) + NROH) * 2048

    The total frame readout time is:
    (((scipix_n + refpix_r + 2) * (512 // scipix_n) + NROH) * 2048 + NFOH)
        * TPIX
    This agrees with the above value of tframe (14.5889 s) if NFOH = 714.
    """
    # Copy in SCI and PIXELDQ arrays for now; that's all we need. The rest
    # of the input model will be copied to output at the end of the step.
    data = output_model.data.copy()
    pixeldq = output_model.pixeldq.copy()

    # Load the reference file data.
    # The reference file data are complex, but they're stored as float, with
    # alternating real and imaginary parts.  We therefore check for twice
    # as many rows as we actually want, and we'll divide that number by two
    # when allocating the arrays alpha and beta.
    nrows = len(irs2_model.irs2_table.field("alpha_0"))
    expected_nrows = 2 * 712 * 2048
    if nrows != expected_nrows:
        log.error(f"Number of rows in reference file = {nrows}, but it should be {expected_nrows}.")
        output_model.meta.cal_step.refpix = "SKIPPED"
        return output_model
    alpha = np.ones((4, nrows // 2), dtype=np.complex64)
    beta = np.zeros((4, nrows // 2), dtype=np.complex64)

    alpha[0, :] = float_to_complex(irs2_model.irs2_table.field("alpha_0"))
    alpha[1, :] = float_to_complex(irs2_model.irs2_table.field("alpha_1"))
    alpha[2, :] = float_to_complex(irs2_model.irs2_table.field("alpha_2"))
    alpha[3, :] = float_to_complex(irs2_model.irs2_table.field("alpha_3"))

    beta[0, :] = float_to_complex(irs2_model.irs2_table.field("beta_0"))
    beta[1, :] = float_to_complex(irs2_model.irs2_table.field("beta_1"))
    beta[2, :] = float_to_complex(irs2_model.irs2_table.field("beta_2"))
    beta[3, :] = float_to_complex(irs2_model.irs2_table.field("beta_3"))

    scipix_n = output_model.meta.exposure.nrs_normal
    if scipix_n is None:
        log.warning(f"Keyword NRS_NORM not found; using default value {scipix_n_default}")
        scipix_n = scipix_n_default

    refpix_r = output_model.meta.exposure.nrs_reference
    if refpix_r is None:
        log.warning(f"Keyword NRS_REF not found; using default value {refpix_r_default}")
        refpix_r = refpix_r_default

    # Convert from sky (DMS) orientation to detector orientation.
    detector = output_model.meta.instrument.detector
    if detector == "NRS1":
        data = np.swapaxes(data, 2, 3)
        pixeldq = np.swapaxes(pixeldq, 0, 1)
    elif detector == "NRS2":
        data = np.swapaxes(data, 2, 3)[:, :, ::-1, ::-1]
        pixeldq = np.swapaxes(pixeldq, 0, 1)[::-1, ::-1]
    else:
        log.warning(f"Detector {detector}; not changing orientation (sky vs detector)")

    n_int = data.shape[0]  # number of integrations in file
    ny = data.shape[-2]  # 2048
    nx = data.shape[-1]  # 3200

    # Create a mask that indicates the locations of normal vs interspersed
    # reference pixels. True flags normal pixels, False is reference pixels.
    irs2_mask = make_irs2_mask(nx, ny, scipix_n, refpix_r)

    # Get bad ref pixel flags from the pixeldq, collapsed along rows
    ref_flags = pixeldq & dqflags.pixel["BAD_REF_PIXEL"]
    ref_flags = np.any(ref_flags, axis=0)

    # If the IRS2 reference file includes data quality info, use that to
    # set bad reference pixel values to zero.
    if hasattr(irs2_model, "dq_table") and len(irs2_model.dq_table) > 0:
        output = irs2_model.dq_table.field("output")
        odd_even = irs2_model.dq_table.field("odd_even")
        mask = irs2_model.dq_table.field("mask")

        # Set interleaved reference pixel values to zero if they are flagged
        # as bad in the DQ extension of the CRDS reference file and not yet handled
        is_irs2 = ~irs2_mask.copy()
        # treat the refout like the other sections
        amplifier = nx // 5
        is_irs2[:amplifier] = is_irs2[2 * amplifier : 3 * amplifier]
        clobber_ref(
            data, output, odd_even, mask, ref_flags, is_irs2, scipix_n=scipix_n, refpix_r=refpix_r
        )
    else:
        log.warning("DQ extension not found in reference file")

    # Compute and apply the correction to one integration at a time
    for integ in range(n_int):
        log.info(f"Working on integration {integ + 1} out of {n_int}")

        # The input data have a length of 3200 for the last axis (X), while
        # the output data have an X axis with length 2048, the same as the
        # Y axis.  This is the reason for the slice `nx-ny:` that is used
        # below.  The last axis of output_model.data should be 2048.
        data0 = data[integ, :, :, :]
        data0 = subtract_reference(
            data0, alpha, beta, irs2_mask, scipix_n, refpix_r, pad, preserve_refpix=preserve_refpix
        )
        if not preserve_refpix:
            data[integ, :, :, nx - ny :] = data0
        else:
            data[integ, :, :, :] = data0

    # Convert corrected data back to sky orientation
    if not preserve_refpix:
        temp_data = data[:, :, :, nx - ny :]
    else:
        temp_data = data
    if detector == "NRS1":
        output_model.data = np.swapaxes(temp_data, 2, 3)
    elif detector == "NRS2":
        output_model.data = np.swapaxes(temp_data[:, :, ::-1, ::-1], 2, 3)
    else:  # don't change orientation
        output_model.data = temp_data

    # Strip interleaved ref pixels from the PIXELDQ and GROUPDQ extensions.
    if not preserve_refpix:
        strip_ref_pixels(output_model, irs2_mask)

    return output_model


def float_to_complex(data):
    """
    Convert real and imaginary parts to complex.

    Parameters
    ----------
    data : ndarray
        Data array with interleaved real and imaginary parts

    Returns
    -------
    data : ndarray
        Complex array made from real and imaginary parts
    """
    nelem = len(data)

    return data[0:-1:2] + 1j * data[1:nelem:2]


def make_irs2_mask(nx, ny, scipix_n, refpix_r):
    """
    Make IRS2 mask.

    Parameters
    ----------
    nx : int
        Number of columns in input data
    ny : int
        Number of rows in input data
    scipix_n : int
        Number of regular samples before stepping out to collect reference samples
    refpix_r : int
        Number of reference samples before stepping back in to collect regular samples

    Returns
    -------
    irs2_mask : ndarray
        The IRS2 mask
    """
    # Number of (scipix_n + refpix_r) per output, assuming four amplifier
    # outputs and one reference output.
    irs2_nx = max((ny, nx))

    # Length of the reference output section.
    refout = irs2_nx // 5
    part = refout - (scipix_n // 2 + refpix_r)
    k = part // (scipix_n + refpix_r)
    # `part` consists of k * (scipix_n + refpix_r) + stuff_at_end
    stuff_at_end = part - k * (scipix_n + refpix_r)

    # Create the mask that flags normal pixels as True.
    irs2_mask = np.ones(irs2_nx, dtype=bool)
    irs2_mask[0:refout] = False

    # Check whether the interspersed reference pixels are in the same
    # locations regardless of readout direction.
    if stuff_at_end == scipix_n // 2:
        # Yes, they are in the same locations.
        for i in range(refout + scipix_n // 2, irs2_nx + 1, scipix_n + refpix_r):
            irs2_mask[i : i + refpix_r] = False
    else:
        # Set the flags for each readout direction separately.
        nelem = refout  # number of elements per output
        temp = np.ones(nelem, dtype=bool)
        for i in range(scipix_n // 2, nelem + 1, scipix_n + refpix_r):
            temp[i : i + refpix_r] = False
        j = refout
        irs2_mask[j : j + nelem] = temp.copy()
        j += nelem
        irs2_mask[j : j + nelem] = temp[::-1].copy()
        j += nelem
        irs2_mask[j : j + nelem] = temp.copy()
        j += nelem
        irs2_mask[j : j + nelem] = temp[::-1].copy()

    return irs2_mask


def strip_ref_pixels(output_model, irs2_mask):
    """
    Copy out the normal pixels from PIXELDQ and GROUPDQ arrays.

    Parameters
    ----------
    output_model : ramp model
        The output science data model, to be modified in-place

    irs2_mask : ndarray of bool
        1D array of length 3200.  True means the element corresponds to a normal
        pixel in the raw, IRS2-format data.  False corresponds either to a reference
        output pixel or to one of the interspersed reference pixel values.
    """
    detector = output_model.meta.instrument.detector

    if detector == "NRS1":
        # Select rows.
        temp_array = output_model.pixeldq
        output_model.pixeldq = temp_array[..., irs2_mask, :]

        temp_array = output_model.groupdq
        output_model.groupdq = temp_array[..., irs2_mask, :]

    elif detector == "NRS2":
        # Reverse the direction of the mask, and select rows.
        temp_mask = irs2_mask[::-1]

        temp_array = output_model.pixeldq
        output_model.pixeldq = temp_array[..., temp_mask, :]

        temp_array = output_model.groupdq
        output_model.groupdq = temp_array[..., temp_mask, :]

    else:
        # Select columns.
        temp_array = output_model.pixeldq
        output_model.pixeldq = temp_array[..., irs2_mask]

        temp_array = output_model.groupdq
        output_model.groupdq = temp_array[..., irs2_mask]


def clobber_ref(data, output, odd_even, mask, ref_flags, is_irs2, scipix_n=16, refpix_r=4):
    """
    Set some interleaved reference pixel values to zero.

    This is an explanation of the arithmetic for computing `ref` in the loop
    over the list of bit numbers that is returned by `decode_mask`.
    Reads of reference pixels are interleaved with reads of science data.  The
    pattern of science pixels (S) and reference pixels (r) looks like this:

    SSSSSSSSrrrrSSSSSSSSSSSSSSSSrrrrSSSSSSSSSSSSSSSSrrrr ... rrrrSSSSSSSS

    Within each amplifier output, a row starts and ends with 8 (scipix_n / 2)
    science pixels, and the row contains 32 blocks of 4 reference pixels.
    There are 20 (scipix_n + refpix_r) pixels from the start of one block of
    reference pixels to the start of the next.  `k` is an integer between
    0 and 31, inclusive, an index to identify the block of reference pixels
    that we need to modify (we'll set two of the pixels to zero).  `odd_even`
    is either 1 or 2, indicating that we should set either the first or the
    second pair of reference pixels to 0.

    The same set of interleaved reference pixels will be set to 0 regardless
    of integration number, group number, or image line number.

    Parameters
    ----------
    data : 4-D ndarray
        The data array in detector orientation.  This includes both the
        science and interleaved reference pixel values.  `data` will be
        modified in-place to set some of the reference pixel values to zero.
        The science data values will not be modified.

    output : 1-D ndarray of int16
        An array of amplifier output numbers, 1, 2, 3, or 4, read from the
        OUTPUT column in the DQ extension of the CRDS reference file.

    odd_even : 1-D ndarray of int16
        An array of integer values, which may be either 1 or 2, read from the
        ODD_EVEN column in the DQ extension of the CRDS reference file.

    mask : 1-D ndarray of uint32
        The MASK column read from the CRDS reference file.

    ref_flags : 1-D ndarray of bool
        Bad reference pixel flags, matching the data row size in detector
        orientation.  True indicates a bad reference pixel.

    is_irs2 : 1-D ndarray of bool
        Array matching the data row size in detector orientation.
        True indicates an interleaved reference pixel.

    scipix_n : int, optional
        Number of regular (science) samples before stepping out to collect
        reference samples.

    refpix_r : int, optional
        Number of reference samples before stepping back in to collect
        regular samples.
    """
    nx = data.shape[-1]  # 3200
    nrows = len(output)

    for row in range(nrows):
        # `offset` is the offset in pixels from the beginning of the row
        # to the start of the current amp output.  `offset` starts with
        # 640 in order to skip over the reference output.
        offset = output[row] * (nx // 5)  # nx // 5 is 640
        # The readout direction alternates from one amp output to the next.
        if output[row] // 2 * 2 == output[row]:
            odd_even_row = 3 - odd_even[row]  # 1 --> 2;  2 --> 1
        else:
            odd_even_row = odd_even[row]
        bits = decode_mask(mask[row])
        log.debug(
            f"output {output[row]}  odd_even {odd_even[row]}  mask {mask[row]} DQ bits {bits}"
        )
        new_bad_pix = []
        for k in bits:
            ref = offset + scipix_n // 2 + k * (scipix_n + refpix_r) + 2 * (odd_even_row - 1)
            log.debug(f"bad interleaved reference at pixels {ref} {ref + 1}")

            # track new bad pixel if not already handled
            for bad_pix in (ref, ref + 1):
                if not ref_flags[bad_pix]:
                    new_bad_pix.append(bad_pix)
                    ref_flags[bad_pix] = True

        # replace new bad pixels
        for bad_pix in new_bad_pix:
            replace_refpix(
                bad_pix,
                data,
                ref_flags,
                is_irs2,
                offset,
                offset + nx // 5,
                scipix_n,
                refpix_r,
                axis=-1,
            )


def decode_mask(mask):
    """
    Interpret the MASK column of the DQ table.

    As per the ESA CDP3 document:
    "There is also a DQ extension that holds a binary table with three
    columns (OUTPUT, ODD_EVEN, and MASK) and eight rows. In the current
    IRS2 implementation, one jumps 32 times to odd and 32 times to even
    reference pixels, which are then read twice consecutively. Therefore,
    the masks are 32 bit unsigned integers that encode bad interleaved
    reference pixels/columns from left to right (increasing column index)
    in the native detector frame. When a bit is set, the corresponding
    reference data should not be used for the correction."

    Parameters
    ----------
    mask : uint32
        A mask value.

    Returns
    -------
    bits : list
        A list of the indices of bits set in the `mask` value.
    """
    # The bit number corresponds to a count of groups of reads of the
    # interleaved reference pixels. The 32-bit unsigned integer encoding
    # has increasing index, from left to right.

    flags = np.array([2**n for n in range(32)], dtype=np.uint32)
    temp = np.bitwise_and(flags, mask)
    bits = np.where(temp > 0)[0]
    bits = list(bits)
    bits.sort()

    return bits


def replace_refpix(
    bad_pix, data, bad_mask, is_irs2, low_limit, high_limit, scipix_n, refpix_r, axis=-2
):
    """
    Replace a bad reference pixel with its nearest neighboring value.

    The nearest reference group above and below the bad pixel
    are checked for good values in pixels with the same parity as the
    bad pixel.  If both are good, they are averaged to determine the replacement
    value.  If only one is good, it is directly used.  If neither are good,
    but a neighboring readout with opposite parity is good, that value is used.
    If none of these options are available, the value is set to 0.0 and will be
    interpolated over during the IRS2 correction.

    The data array is modified in place.

    Parameters
    ----------
    bad_pix : int
        The bad pixel index to replace.
    data : 4-D ndarray
        The data array containing reference and science pixels.
        If in science orientation, `axis` should be -2. If in detector
        orientation, `axis` should be set to -1.
    bad_mask : 1-D ndarray
        A boolean mask, where True indicates a bad reference value.
        Should match the shape of the data along `axis`.
    is_irs2 : 1-D ndarray
        A boolean mask, where True indicates an interleaved reference pixel.
        Should match the shape of the data along `axis`.
    low_limit : int
        The lower limit of the data indices along `axis` to check
        for replacement values, usually set to the bottom of the amplifier.
    high_limit : int
        The upper limit of the data indices along `axis` to check
        for replacement values, usually set to the bottom of the amplifier.
    scipix_n : int
        Number of regular (science) samples before stepping out to collect
        reference samples.
    refpix_r : int
        Number of reference samples before stepping back in to collect
        regular samples.
    axis : int, optional
        Indicates the axis containing the reference pixel values.
        Set to -2 for science orientation, -1 for detector orientation.
    """
    # nearest reference pixel group, respecting parity
    ref_period = scipix_n + refpix_r
    nearest_low = bad_pix - ref_period
    nearest_high = bad_pix + ref_period

    # check for neighboring good ref pixels
    # to use as a fallback value
    neighbor = None
    lower_neighbor = bad_pix - 2
    upper_neighbor = bad_pix + 2
    if lower_neighbor >= low_limit and is_irs2[lower_neighbor] and ~bad_mask[lower_neighbor]:
        neighbor = lower_neighbor
    elif upper_neighbor < high_limit and is_irs2[upper_neighbor] and ~bad_mask[upper_neighbor]:
        neighbor = upper_neighbor

    if neighbor is not None:
        if axis == -1:
            replace_value = data[:, :, :, neighbor]
        else:
            replace_value = data[:, :, neighbor, :]
    else:
        # last resort: set to zero and allow cosine interpolation
        replace_value = 0.0

    # try to average upper and lower
    v1, v2 = None, None
    if nearest_low >= low_limit and ~bad_mask[nearest_low]:
        if axis == -1:
            v1 = data[:, :, :, nearest_low]
        else:
            v1 = data[:, :, nearest_low, :]
    if nearest_high < high_limit and ~bad_mask[nearest_high]:
        if axis == -1:
            v2 = data[:, :, :, nearest_high]
        else:
            v2 = data[:, :, nearest_high, :]

    if v1 is not None and v2 is not None:
        log.debug(
            f"   Pixel {bad_pix} replaced with value averaged from {nearest_low},{nearest_high}"
        )
        replace_value = np.mean([v1, v2], axis=0)
    elif v1 is not None:
        log.debug(f"   Pixel {bad_pix} replaced with value at {nearest_low}")
        replace_value = v1
    elif v2 is not None:
        log.debug(f"   Pixel {bad_pix} replaced with value at {nearest_high}")
        replace_value = v2
    elif neighbor is not None:
        log.debug(f"   Pixel {bad_pix} replaced with value at neighbor {neighbor}")
    else:
        log.debug(f"   Pixel {bad_pix} replaced with 0.0")

    if axis == -1:
        data[:, :, :, bad_pix] = replace_value
    else:
        data[:, :, bad_pix, :] = replace_value


def flag_bad_refpix(datamodel, n_sigma=3.0, flag_only=False, replace_only=False):
    """
    Flag bad reference pixels and replace with nearest good values.

    Parameters
    ----------
    datamodel : DataModel
        The data, in science orientation.  This includes both the
        science and interleaved reference pixel values.  Data and pixeldq
        will be modified in-place. The science data values will not be
        modified.
    n_sigma : float, optional
        Flagging threshold, expressed as a factor times the standard deviation.
    flag_only : bool, optional
        If set, bad values will be flagged in the pixeldq extension but
        not replaced.
    replace_only : bool, optional
        If set, previously flagged bad values will be replaced, but new outliers
        will not be flagged.
    """
    data = datamodel.data
    pixeldq = datamodel.pixeldq
    scipix_n = datamodel.meta.exposure.nrs_normal
    refpix_r = datamodel.meta.exposure.nrs_reference
    log.debug(f"Using flagging threshold n_sigma = {n_sigma}")

    # bad pixels will be replaced for all integrations and all groups
    nints, ngroups, ny, nx = np.shape(data)

    # initialize the mask with any previously marked bad pixels
    ref_flags = pixeldq & dqflags.pixel["BAD_REF_PIXEL"]
    mask_bad = np.any(ref_flags, axis=1)
    is_irs2 = np.full(ny, False)

    # calculate differences of readout pairs per amplifier
    amplifier = ny // 5  # 640
    ref_period = scipix_n + refpix_r
    initial_mask = mask_bad.copy()
    for k in range(5):
        offset = int(k * amplifier)

        # get statistics for each integration individually, but
        # apply flags to all integrations
        for j in range(nints):
            ref_pix, rp_diffs, rp_means, rp_stds = [], [], [], []
            int_bad = initial_mask.copy()

            # jump from the start of the reference pixel sequence to the next
            # starting pixel is from 8 to 640 by 20
            for rpstart in range(scipix_n // 2, amplifier, ref_period):
                # amplifier offset
                rpstart += offset

                # go through the reference pixels by pairs
                for ri in range(0, refpix_r, 2):
                    ri = rpstart + ri
                    rp_d = np.mean(np.abs(data[j, :, ri + 1, :] - data[j, :, ri, :]))
                    rp_m = np.mean(data[j, :, ri : ri + 2, :])
                    rp_s = np.std(data[j, :, ri : ri + 2, :])
                    is_irs2[ri : ri + 2] = True

                    # exclude ref pix already flagged
                    good = ~np.any(int_bad[ri : ri + 2])
                    if good and not replace_only:
                        ref_pix.append(ri)
                        rp_means.append(rp_m)
                        rp_stds.append(rp_s)
                        rp_diffs.append(rp_d)

            if not replace_only:
                ref_pix = np.array(ref_pix, dtype=int)
                rp_diffs = np.array(rp_diffs)
                rp_means = np.array(rp_means)
                rp_stds = np.array(rp_stds)
                pair_pixel = ref_pix + 1

                # clipped stats for all tests
                mean_of_diffs, _, std_of_diffs = sigma_clipped_stats(rp_diffs, sigma=n_sigma)
                mean_of_means, _, std_of_means = sigma_clipped_stats(rp_means, sigma=n_sigma)
                mean_of_stds, _, std_of_stds = sigma_clipped_stats(rp_stds, sigma=n_sigma)

                # find the additional intermittent bad pixels, marking both readouts
                high_diffs = (rp_diffs - mean_of_diffs) > (n_sigma * std_of_diffs)
                high_means = (rp_means - mean_of_means) > (n_sigma * std_of_means)
                high_stds = (rp_stds - mean_of_stds) > (n_sigma * std_of_stds)

                log.debug(
                    f"High diffs={np.sum(high_diffs)}, "
                    f"high means={np.sum(high_means)}, "
                    f"high stds={np.sum(high_stds)}"
                )
                int_bad[ref_pix[high_diffs]] = True
                int_bad[pair_pixel[high_diffs]] = True
                int_bad[ref_pix[high_means]] = True
                int_bad[pair_pixel[high_means]] = True
                int_bad[ref_pix[high_stds]] = True
                int_bad[pair_pixel[high_stds]] = True

                log.debug(
                    f"{np.sum(int_bad[offset : offset + amplifier])} "
                    f"suspicious bad reference pixels in "
                    f"amplifier {k}, integration {j}"
                )
                mask_bad |= int_bad

        # replace any flagged pixels if desired
        if not flag_only:
            # list of all bad pixels
            all_bad = np.arange(offset, offset + amplifier)[mask_bad[offset : offset + amplifier]]
            for bad_pix in all_bad:
                replace_refpix(
                    bad_pix, data, mask_bad, is_irs2, offset, offset + amplifier, scipix_n, refpix_r
                )

    if flag_only:
        log.info(f"Total bad reference pixels flagged: {np.sum(mask_bad)}")
    else:
        log.info(f"Total bad reference pixels replaced: {np.sum(mask_bad)}")
    if pixeldq is not None:
        pixeldq[mask_bad] |= dqflags.pixel["BAD_REF_PIXEL"] | dqflags.pixel["DO_NOT_USE"]


def subtract_reference(
    data0, alpha, beta, irs2_mask, scipix_n, refpix_r, pad, preserve_refpix=False
):
    """
    Subtract reference output and pixels for the current integration.

    Parameters
    ----------
    data0 : ndarray
        The science data for the current integration.  The shape is
        expected to be (ngroups, ny, 3200), where ngroups is the number of
        groups, and ny is the pixel height of the image.  The width 3200
        of the image includes the "normal" pixel data, plus the embedded
        reference pixels, and the reference output.
    alpha : ndarray
        This is a 2-D array of values read from the reference file.  The
        first axis is the sector number (but only for the normal pixel
        data and reference pixels, not the reference output).  The second
        axis has length 2048 * 712, corresponding to the time-ordered
        arrangement of the data.  For each sector, the correction is
        applied as follows:  data * alpha[i] + reference_output * beta[i].
    beta : ndarray
        Data read from the reference file.  See `alpha` for details.
    irs2_mask : Boolean, 1-D array
        True means the element corresponds to a normal pixel in the raw,
        IRS2-format data.  False corresponds either to a reference output
        pixel or to one of the interspersed reference pixel values.
    scipix_n : int
        Number of regular samples before stepping out to collect
        reference samples.
    refpix_r : int
        Number of reference samples before stepping back in to collect
        regular samples.
    pad : int
        The effective number of pixels sampled during the pause at the end
        of each row (new-row overhead).
    preserve_refpix : bool
        If True, reference pixels will be preserved in the output.
        This is not used in the science pipeline, but is necessary to
        create new bias files for IRS2 mode.

    Returns
    -------
    data0 : ndarray
        The science data for the current integration, with reference output
        and embedded reference pixels subtracted and also removed, leaving
        only the normal pixel data (including the reference pixels on each
        edge).  The shape is expected to be (ngroups, ny, nx), where
        nx = ny = 2048.
    """
    shape = data0.shape
    ngroups = shape[0]
    ny = shape[1]
    nx = shape[2]

    # See expression in equation 1 in IRS2_Handoff.pdf.
    # row = 712, if scipix_n = 16, refpix_r = 4, pad = 8.
    row = (scipix_n + refpix_r + 2) * 512 // scipix_n + pad

    # s = size(data0)
    # If data0 is the data for one integration, then:
    # s[0] would be 3
    # s[1] = shape[2] = nx, the length of the X axis
    # s[2] = shape[1] = ny, the length of the Y axis
    # s[3] = shape[0] = ngroups, the number of groups (or frames)

    ind_n = np.arange(512, dtype=np.intp)
    ind_ref = np.arange(512 // scipix_n * refpix_r, dtype=np.intp)

    # hnorm is an array of column indices of normal pixels.
    # len(hnorm) = 512; len(href) = 128
    # len(hnorm1) = 512; len(href1) = 128
    hnorm = ind_n + refpix_r * ((ind_n + scipix_n // 2) // scipix_n)

    # href is an array of column indices of reference pixels.
    href = ind_ref + scipix_n * (ind_ref // refpix_r) + scipix_n // 2

    hnorm1 = ind_n + (refpix_r + 2) * ((ind_n + scipix_n // 2) // scipix_n)
    href1 = ind_ref + (scipix_n + 2) * (ind_ref // refpix_r) + scipix_n // 2 + 1

    unpad = np.sort(np.hstack([hnorm1, href1]))

    # Subtract the average over the ramp for each pixel.
    # b_offset is saved so that it can be added back in at the end.
    b_offset = data0.sum(axis=0, dtype=np.float64) / float(ngroups)
    data0 -= b_offset

    # IDL:  data0 = reform(data0, s[1]/5, 5, s[2], s[3], /over)
    #                             nx/5,   5, ny,   ngroups    (IDL)
    data0 = data0.reshape((ngroups, ny, 5, nx // 5))

    # current order:  nx/5, 5, ny, ngroups    (IDL)
    # current order:  ngroups, ny, 5, nx/5    (numpy)
    #                 0        1   2  3       current numpy indices
    # transpose to:   nx/5, ny, ngroups, 5    (IDL)
    # transpose to:   5, ngroups, ny, nx/5    (numpy)
    #                 2  0        1   3       transpose order for numpy
    # Therefore:      0 1 2 3  -->  2 0 1 3   transpose order for numpy
    # Here is another way to look at it:
    # IDL:    0 1 2 3  -->  0 2 3 1
    #         3 2 1 0       1 3 2 0 (IDL indices, but reversed to numpy order)
    # numpy:  0 1 2 3  -->  2 0 1 3
    # IDL:  data0 = transpose(data0, [0,2,3,1])
    data0 = np.transpose(data0, (2, 0, 1, 3))

    # Flip the direction of the X axis for every other output, so the readout
    # direction in data0 will be the same for every output.
    data0[0, :, :, :] = data0[0, :, :, ::-1]
    data0[2, :, :, :] = data0[2, :, :, ::-1]
    data0[4, :, :, :] = data0[4, :, :, ::-1]

    # convert to time sequences of normal pixels and reference pixels.
    # IDL:  d0 = fltarr(s[1] / 5 + pad + 2 * (512 / scipix_n), s[2], s[3], 5)
    # Note:  nx // 5 + pad + 2 * (512 // scipix_n) = 640 + 64 + 8 = 712.
    # hnorm1[-1] = 703, and hnorm[-1] = 639, so 703 - 639 = 64.
    # 8 is the pad value.

    d0 = np.zeros((5, ngroups, ny, row), dtype=np.float32)  # (5, ngroups, 2048, 712)
    # IDL:  d0[hnorm1,*,*,*] = data0[hnorm,*,*,*]
    # IDL:  d0[href1,*,*,*] = data0[href,*,*,*]
    # IDL:  data0 = temporary(d0)
    d0[:, :, :, hnorm1] = data0[:, :, :, hnorm]
    d0[:, :, :, href1] = data0[:, :, :, href]
    del data0
    data0 = d0.copy()
    del d0

    # Fitting and removal of slopes per frame to remove issues at frame boundaries
    remove_slopes(data0, ngroups, ny, row)

    # Use cosine weighted interpolation to replace 0.0 values and bad
    # pixels and gaps. (initial guess)
    replace_bad_pixels(data0, ngroups, ny, row)

    # Fill in bad pixels, gaps, and reference data locations in the normal
    # data, using Fourier filtering/interpolation
    fill_bad_regions(data0, ngroups, ny, nx, row, scipix_n, refpix_r, pad, hnorm, hnorm1)

    # Setup various lists of indices that will be used in subsequent
    # sections for keeping/shuffling reference pixels in various arrays
    #
    # The comments are for scipix_n = 16, refpix_r = 4
    n0 = 512 // scipix_n
    n1 = scipix_n + refpix_r + 2
    ht = np.arange(n0 * n1, dtype=np.int32).reshape((n0, n1))  # (32, 22)
    ht[:, 0 : (scipix_n - refpix_r) // 2 + 1] = -1
    ht[:, scipix_n // 2 + 1 + 3 * refpix_r // 2 :] = -1
    hs = ht.copy()
    # ht is like href1, but extended over gaps and first and last norm pix
    mask = ht >= 0
    ht = ht[mask]  # 1-D, length = 2 * refpix_r * 512 / scipix_n

    # IDL:  hs[scipix_n/2 + 1-refpix_r/2:scipix_n/2 + refpix_r + refpix_r/2,*] =
    #       hs[reform([transpose(reform(indgen(refpix_r),refpix_r/2,2)),
    #           transpose(reform(indgen(refpix_r),refpix_r/2,2))],refpix_r * 2)
    #           + scipix_n/2 + 1,*]  ; WIRED for R=2^(int)

    indr = np.arange(refpix_r, dtype=np.intp).reshape((2, refpix_r // 2))
    # indr_t =
    # [[0 2]
    #  [1 3]]
    indr_t = indr.transpose()

    # Before flattening, two_indr_t =
    # [[0 2 0 2]
    #  [1 3 1 3]]
    # After flattening, two_indr_t = [0 2 0 2 1 3 1 3].
    two_indr_t = np.concatenate((indr_t, indr_t), axis=1).flatten()
    two_indr_t += scipix_n // 2 + 1  # [9 11 9 11 10 12 10 12]
    hs[:, scipix_n // 2 + 1 - refpix_r // 2 : scipix_n // 2 + 1 + refpix_r // 2 + refpix_r] = hs[
        :, two_indr_t
    ]
    mask = hs >= 0
    hs = hs[mask]  # hs is now 1-D

    if refpix_r % 4 == 2:
        len_hs = len(hs)
        temp_hs = hs.reshape(len_hs // 2, 2)
        temp_hs = temp_hs[:, ::-1]
        hs = temp_hs.flatten()

    # Construct the reference data: this is done in a big loop over the
    # four "sectors" of data in the image, corresponding to the amp regions.
    # Data from each sector is operated on independently and ultimately
    # the corrections are subtracted from each sector independently.
    shape_d = data0.shape
    for k in range(1, 5):
        log.debug(f"processing sector {k}")

        # At this point in the processing data0 has shape (5, ngroups, 2048, 712),
        # assuming normal IRS2 readout settings. r0k contains a subset of the
        # data from 1 sector of data0, with shape (ngroups, 2048, 256)
        r0k = np.zeros((shape_d[1], shape_d[2], shape_d[3]), dtype=np.float32)
        temp = data0[k, :, :, hs].copy()
        temp = np.transpose(temp, (1, 2, 0))
        r0k[:, :, ht] = temp
        del temp

        # data0 has shape (5, ngroups, ny, row).  See the section above where
        # d0 was created, then copied (moved) to data0.
        # sd[1] = shape_d[3]   row (712)
        # sd[2] = shape_d[2]   ny (2048)
        # sd[3] = shape_d[1]   ngroups
        # sd[4] = shape_d[0]   5
        # s is used below, so for convenience, here are the values again:
        # s[1] = shape[2] = nx
        # s[2] = shape[1] = ny
        # s[3] = shape[0] = ngroups

        # IDL and numpy differ in where they apply the normalization for the
        # FFT.  This really shouldn't matter.
        normalization = float(shape_d[2] * shape_d[3])

        # Set up refout if alpha was provided
        refout0 = None
        if alpha is not None:
            # IDL:  refout0 = reform(data0[*,*,*,0], sd[1] * sd[2], sd[3])
            refout0 = data0[0, :, :, :].reshape((shape_d[1], shape_d[2] * shape_d[3]))

            # IDL:  refout0 = fft(refout0, dim=1, /over)
            # Divide by the length of the axis to be consistent with IDL.
            refout0 = np.fft.fft(refout0, axis=1) / normalization

        # IDL:  r0 = reform(r0, sd[1] * sd[2], sd[3], 5, /over)
        r0k = r0k.reshape((shape_d[1], shape_d[2] * shape_d[3]))
        r0k = r0k.astype(np.complex64)
        r0k_fft = np.fft.fft(r0k, axis=1) / normalization

        # Note that where the IDL code uses alpha, we use beta, and vice versa.
        # IDL:  for k=0,3 do oBridge[k]->Execute,
        #           "for i=0, s3-1 do r0[*,i] *= alpha"
        r0k_fft *= beta[k - 1]

        # IDL:  for k=0,3 do oBridge[k]->Execute,
        #           "for i=0, s3-1 do r0[*,i] += beta * refout0[*,i]"
        if alpha is not None:
            r0k_fft += alpha[k - 1] * refout0
        del refout0

        # IDL:  for k=0,3 do oBridge[k]->Execute,
        #           "r0 = fft(r0, 1, dim=1, /overwrite)", /nowait
        r0k = np.fft.ifft(r0k_fft, axis=1) * normalization
        del r0k_fft

        # sd[1] = shape_d[3]   row (712)
        # sd[2] = shape_d[2]   ny (2048)
        # sd[3] = shape_d[1]   ngroups
        # sd[4] = shape_d[0]   5
        # IDL:  r0 = reform(r0, sd[1], sd[2], sd[3], 5, /over)
        r0k = r0k.reshape(shape_d[1], shape_d[2], shape_d[3])
        r0k = r0k.real
        if not preserve_refpix:
            r0k = r0k[:, :, hnorm1]
        else:
            r0k = r0k[:, :, unpad]

        # Subtract the correction from the data in this sector
        if not preserve_refpix:
            data0[k, :, :, hnorm1] -= np.transpose(r0k, (2, 0, 1))
        else:
            data0[k, :, :, unpad] -= np.transpose(r0k, (2, 0, 1))
        del r0k

    # End of loop over 4 sectors

    # Original data0 array has shape (5, ngroups, 2048, 712). Now that
    # correction has been applied, remove the interleaved reference pixels.
    # This leaves data0 with shape (5, ngroups, 2048, 512).
    if not preserve_refpix:
        data0 = data0[:, :, :, hnorm1]
    else:
        data0 = data0[:, :, :, unpad]

    # Unflip the data in the sectors that have opposite readout direction
    if preserve_refpix:
        data0[0, :, :, :] = data0[0, :, :, ::-1]
    data0[2, :, :, :] = data0[2, :, :, ::-1]
    data0[4, :, :, :] = data0[4, :, :, ::-1]

    # IDL:  data0 = transpose(data0, [0,3,1,2])  0, 1, 2, 3 --> 0, 3, 1, 2
    # current order:  512, ny, ngroups, 5     (IDL)
    # current order:  5, ngroups, ny, 512     (numpy)
    #                 0  1        2   3       current numpy indices
    # transpose to:   512, 5, ny, ngroups     (IDL)
    # transpose to:   ngroups, ny, 5, 512     (numpy)
    #                 1        2   0  3       transpose order for numpy
    # Therefore:      0 1 2 3  -->  1 2 0 3   transpose order for numpy
    # After transposing, data0 will have shape (ngroups, 2048, 5, 512).
    data0 = np.transpose(data0, (1, 2, 0, 3))

    # Reshape data0 back to its normal (ngroups, 2048, 2048), which has
    # the interleaved reference pixels stripped out.
    # IDL:  data0 = reform(data0[*, 1:*, *, *], s[2], s[2], s[3], /over)
    # Note:  ny x ny, not ny x nx.
    if not preserve_refpix:
        data0 = data0[:, :, 1:, :].reshape((ngroups, ny, ny))
    else:
        data0 = data0.reshape((ngroups, ny, nx))

    # b_offset is the average over the ramp that we subtracted near the
    # beginning; add it back in.
    # Shape of b_offset is (2048, 3200), but data0 is (ngroups, 2048, 2048),
    # so a mask is applied to b_offset to remove the reference pix locations.
    if not preserve_refpix:
        data0 += b_offset[..., irs2_mask]
    else:
        # add in only data value -
        # reference mean should be subtracted if not stripped,
        # except in reference sector
        data0[..., irs2_mask] += b_offset[..., irs2_mask]
        data0[..., : nx // 5] += b_offset[..., : nx // 5]

    return data0


def fft_interp_norm(dd0, mask0, row, hnorm, hnorm1, ny, ngroups, aa, n_iter_norm):
    """
    Filter iteratively in FFT space of the normal pixels in each group.

    Parameters
    ----------
    dd0 : ndarray
        Data array containing all groups, updated in place
    mask0 : ndarray
        Mask for pixels to filter, with dimensions ny x nrow. 1 means use the pixel,
        0 means do not use it.
    row : int
        Row size. Computed from the number of science pixels, reference pixels
        and padding in an amplifier.
    hnorm : ndarray
        Array of column indices for normal pixels.
    hnorm1 : ndarray
        Shifted index values for normal pixels.
    ny : int
        Y size of data array.
    ngroups : int
        Number of groups.
    aa : ndarray
        Filter to apply.
    n_iter_norm : int
        Number of filtering iterations.
    """
    mm = np.zeros((ny, row), dtype=np.int8)
    mm[:, hnorm1] = mask0[:, hnorm]
    hm = mm != 0  # 2-D boolean mask
    for j in range(ngroups):
        dd = dd0[j, :, :].copy()  # make a copy, not a view
        p = dd.flatten()
        for _it in range(n_iter_norm):
            pp = np.fft.fft(p)
            pp *= aa
            p[:] = np.fft.ifft(pp).real
            p[hm.ravel()] = dd[hm]
        dd0[j, :, :] = p.reshape((ny, row))


def ols_line(x, y):
    """
    Fit a straight line using ordinary least squares.

    Parameters
    ----------
    x : ndarray
        Array of independent variables
    y : ndarray
        Array of dependent variables

    Returns
    -------
    intercept : float
        Intercept of straight line fit
    slope : float
        Slope of straight line fit
    """
    xf = x.ravel()
    yf = y.ravel()
    if len(xf) < 1 or len(yf) < 1:
        return 0.0, 0.0

    groups = float(len(xf))
    mean_x = xf.mean()
    mean_y = yf.mean()
    sum_x2 = (xf**2).sum()
    sum_xy = (xf * yf).sum()

    slope = (sum_xy - groups * mean_x * mean_y) / (sum_x2 - groups * mean_x**2)
    intercept = mean_y - slope * mean_x

    return intercept, slope


def remove_slopes(data0, ngroups, ny, row):
    """
    Remove slopes.

    Fitting and removal of slopes per frame to remove issues at frame boundaries.

    Parameters
    ----------
    data0 : ndarray
        Input data array
    ngroups : int
        Number of groups in input data
    ny : int
        Number of rows in input data
    row : int
        Row size
    """
    time_arr = np.arange(ny * row, dtype=np.float32).reshape((ny, row))
    time_arr -= time_arr.mean(dtype=np.float64)
    row4plus4 = np.array([0, 1, 2, 3, 2044, 2045, 2046, 2047], dtype=np.intp)

    # For ab_3, it should be OK to use the same index order as the IDL code.
    ab_3 = np.zeros((2, ngroups, 5), dtype=np.float32)
    for i in range(5):
        for k in range(ngroups):
            # mask is 2-D, since both row4plus4 and : have more than one element.
            mask = data0[i, k, row4plus4, :] != 0.0
            (intercept, slope) = ols_line(
                time_arr[row4plus4, :][mask], data0[i, k, row4plus4, :][mask]
            )
            ab_3[0, k, i] = intercept
            ab_3[1, k, i] = slope

    for i in range(5):
        for k in range(ngroups):
            # weight is 0 where data0 is 0, else 1.
            weight = (data0[i, k, :, :] != 0.0).astype(np.int8)
            data0[i, k, :, :] -= (ab_3[0, k, i] + time_arr * ab_3[1, k, i]) * weight


def replace_bad_pixels(data0, ngroups, ny, row):
    """
    Replace bad pixels.

    Use cosine weighted interpolation to replace 0.0 values and bad
    pixels and gaps.

    s[1] = nx  s[2] = ny  s[3] = ngroups

    Parameters
    ----------
    data0 : ndarray
        Input data array
    ngroups : int
        Number of groups in input data array
    ny : int
        Number of rows in input data array
    row : int
        Row definition - row = (scipix_n + refpix_r + 2) * 512 // scipix_n + pad
        row = 712, if scipix_n = 16, refpix_r = 4, pad = 8
    """
    w_ind = np.arange(1, 32, dtype=np.float32) / 32.0
    w = np.sin(w_ind * np.pi)
    kk = 0
    for jj in range(ngroups):
        dat = data0[kk, jj, :, :].reshape(row * ny)
        mask = (dat != 0.0).astype(np.float32)
        numerator = convolve1d(dat, w, mode="wrap")
        denominator = convolve1d(mask, w, mode="wrap")
        div_zero = denominator == 0.0  # check for divide by zero
        numerator = np.where(div_zero, 0.0, numerator)
        denominator = np.where(div_zero, 1.0, denominator)
        dat = numerator / denominator
        dat = dat.reshape(ny, row)
        mask = mask.reshape(ny, row)
        data0[kk, jj, :, :] += dat * (1.0 - mask)


def fill_bad_regions(data0, ngroups, ny, nx, row, scipix_n, refpix_r, pad, hnorm, hnorm1):
    """
    Fill the bad regions in the data.

    Use Fourier filter/interpolation to replace
      (a) bad pixel, gaps, and reference data in the time-ordered normal data
      (b) gaps and normal data in the time-ordered reference data
    This "improves" upon the cosine interpolation performed above.

    Parameters
    ----------
    data0 : ndarray
        Input data array.  Modified in place.
    ngroups : int
        Number of groups in input data
    ny : int
        Number of rows in input data
    nx : int
        Number of columns in input data
    row : int
        Row definition: row = (scipix_n + refpix_r + 2) * 512 // scipix_n + pad
        row = 712, if scipix_n = 16, refpix_r = 4, pad = 8
    scipix_n : int
        Number of regular samples before stepping out to collect reference samples
    refpix_r : int
        Number of reference samples before stepping back in to collect regular samples
    pad : int
        The effective number of pixels sampled during the pause at the end
        of each row (new-row overhead).  The padding is needed to preserve
        the phase of temporally periodic signals.
    hnorm : ndarray
        Array of column indices for normal pixels
    hnorm1 : ndarray
        Shifted index values for normal pixels
    """
    # Parameters for the filter to be used:
    # length of apodization cosine filter
    elen = 110000 // (scipix_n + refpix_r + 2)

    # max unfiltered frequency
    blen = (512 + 512 // scipix_n * (refpix_r + 2) + pad) // (
        scipix_n + refpix_r + 2
    ) * ny // 2 - elen // 2

    # Construct the filter [1, cos, 0, cos, 1].
    temp_a1 = (np.cos(np.arange(elen, dtype=np.float64) * np.pi / float(elen)) + 1.0) / 2.0

    # elen = 5000
    # blen = 30268
    # row * ny // 2 - 2 * blen - 2 * elen = 658552
    # len(temp_a2) = 729088
    temp_a2 = np.concatenate(
        (
            np.ones(blen, dtype=np.float64),
            temp_a1.copy(),
            np.zeros(row * ny // 2 - 2 * blen - 2 * elen, dtype=np.float64),
            temp_a1[::-1].copy(),
            np.ones(blen, dtype=np.float64),
        )
    )

    roll_a2 = np.roll(temp_a2, -1)
    aa = np.concatenate((temp_a2, roll_a2[::-1]))
    del temp_a1, temp_a2, roll_a2

    # IDL:  aa = a # replicate(1, s[3]) ; for application to the data
    # In IDL, aa is a 2-D array with one column of `a` for each group.  In
    # Python, numpy broadcasting takes care of this.

    n_iter_norm = 3
    dd0 = data0[0, :, :, :]
    # IDL:  fft_interp_norm, dd0, 2, replicate(1, s[1] / 4, s[2], 4),
    #                        row, hnorm, hnorm1, s, aa , n_iter_norm
    fft_interp_norm(
        dd0,
        np.ones((ny, nx // 4), dtype=np.int64),
        row,
        hnorm,
        hnorm1,
        ny,
        ngroups,
        aa,
        n_iter_norm,
    )

    data0[0, :, :, :] = dd0.copy()
    del aa, dd0
