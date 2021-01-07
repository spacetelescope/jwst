"""
Two-Point Difference method for finding outliers in a 4-D ramp data array.
The scheme used in this variation of the method uses numpy array methods
to compute first-differences and find the max outlier in each pixel while
still working in the full 4-D data array. This makes detection of the first
outlier very fast. We then iterate pixel-by-pixel over only those pixels
that are already known to contain an outlier, to look for any additional
outliers and set the appropriate DQ mask for all outliers in the pixel.
This is MUCH faster than doing all the work on a pixel-by-pixel basis.
"""

import logging
import numpy as np
from ..datamodels import dqflags

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

HUGE_NUM = np.finfo(np.float32).max

JUMP_DET = dqflags.group["JUMP_DET"]
SATURATED = dqflags.group["SATURATED"]
DO_NOT_USE = dqflags.group["DO_NOT_USE"]


def find_crs(data, group_dq, read_noise, rej_threshold, nframes, flag_4_neighbors,
             max_jump_to_flag_neighbors, min_jump_to_flag_neighbors):
    """
    Find CRs/Jumps in each integration within the input data array.
    The input data array is assumed to be in units of electrons, i.e. already
    multiplied by the gain. We also assume that the read noise is in units of
    electrons.
    """
    gdq = group_dq.copy()

    # Get data characteristics
    nints, ngroups, nrows, ncols = data.shape
    ndiffs = ngroups - 1

    # Create arrays for output
    row_above_gdq = np.zeros((nints, ngroups, ncols), dtype=np.uint8)
    row_below_gdq = np.zeros((nints, ngroups, ncols), dtype=np.uint8)

    # Square the read noise values, for use later
    read_noise_2 = read_noise ** 2

    # Set saturated values in the input data array to NaN, so they don't get
    # used in any of the subsequent calculations
    data[np.where(np.bitwise_and(gdq, SATURATED))] = np.nan

    # Set pixels flagged as DO_NOT_USE in the input to NaN, so they don't get
    # used in any of the subsequent calculations. MIRI exposures can sometimes have
    # all pixels in the first and last groups flagged with DO_NOT_USE.
    data[np.where(np.bitwise_and(gdq, DO_NOT_USE))] = np.nan

    # Loop over multiple integrations
    for integ in range(nints):

        log.info(f'Working on integration {integ+1}:')

        # Compute first differences of adjacent groups up the ramp
        # note: roll the ngroups axis of data array to the end, to make
        # memory access to the values for a given pixel faster.
        # New form of the array has dimensions [nrows, ncols, ngroups].
        first_diffs = np.diff(np.rollaxis(data[integ], axis=0, start=3), axis=2)
        positive_first_diffs = np.abs(first_diffs)

        # sat_groups is a 3D array that is true when the group is saturated
        sat_groups = np.isnan(positive_first_diffs)

        # number_sat_groups is a 2D array with the count of saturated groups for each pixel
        number_sat_groups = sat_groups.sum(axis=2)

        # Make all the first diffs for saturated groups be equal to
        # 100,000 to put them above the good values in the sorted index
        first_diffs[np.isnan(first_diffs)] = 100000.

        # Here we sort the 3D array along the last axis, which is the group axis.
        # np.argsort returns a 3D array with the last axis containing the indices
        # that would yield the groups in order.
        sort_index = np.argsort(positive_first_diffs)
        del positive_first_diffs

        # median_diffs is a 2D array with the clipped median of each pixel
        median_diffs = get_clipped_median(ndiffs, number_sat_groups, first_diffs, sort_index)

        # Compute uncertainties as the quadrature sum of the poisson noise
        # in the first difference signal and read noise. Because the first
        # differences can be biased by CRs/jumps, we use the median signal
        # for computing the poisson noise. Here we lower the read noise
        # by the square root of number of frames in the group.
        # Sigma is a 2D array.
        poisson_noise = np.sqrt(np.abs(median_diffs))
        sigma = np.sqrt(poisson_noise * poisson_noise + read_noise_2 / nframes)
        del poisson_noise

        # Reset sigma to exclude pixels with both readnoise and signal=0
        sigma_0_pixels = np.where(sigma == 0.)
        if len(sigma_0_pixels[0] > 0):
            log.debug(f'Twopt found {len(sigma_0_pixels[0])} pixels with sigma=0')
            log.debug('which will be reset so that no jump will be detected')
            sigma[sigma_0_pixels] = HUGE_NUM
        del sigma_0_pixels

        # Compute distance of each sample from the median in units of sigma;
        # note that the use of "abs" means we'll detect positive and negative outliers.
        # ratio is a 2D array with the units of sigma deviation of the difference from the median.
        ratio = np.abs(first_diffs - median_diffs[:, :, np.newaxis]) / sigma[:, :, np.newaxis]
        del sigma
        del median_diffs

        # Get the group index for each pixel of the largest non-saturated group,
        # assuming the indices are sorted. 2 is subtracted from ngroups because we are using
        # differences and there is one less difference than the number of groups.
        # This is a 2-D array.
        max_value_index = ngroups - 2 - number_sat_groups

        # Extract from the sorted group indices the index of the largest non-saturated group.
        row, col = np.where(number_sat_groups >= 0)
        max_index1d = sort_index[row, col, max_value_index[row, col]]
        max_index1 = np.reshape(max_index1d, (nrows, ncols))  # reshape to a 2-D array
        del max_index1d
        del max_value_index

        # Get the row and column indices of pixels whose largest non-saturated ratio is above the threshold
        r, c = np.indices(max_index1.shape)
        row1, col1 = np.where(ratio[r, c, max_index1] > rej_threshold)
        del max_index1
        log.info(f'From highest outlier Two-point found {len(row1)} pixels with at least one CR')

        # Loop over all pixels that we found the first CR in
        number_pixels_with_cr = len(row1)
        for j in range(number_pixels_with_cr):

            # Extract the first diffs for this pixel with at least one CR, yielding a 1D array
            pixel_masked_diffs = first_diffs[row1[j], col1[j]]

            # Get the scalar readnoise^2 and number of saturated groups for this pixel.
            pixel_rn2 = read_noise_2[row1[j], col1[j]]
            pixel_sat_groups = number_sat_groups[row1[j], col1[j]]

            # Create a CR mask and set 1st CR to be found
            # cr_mask=0 designates a CR
            pixel_cr_mask = np.ones(pixel_masked_diffs.shape, dtype=bool)
            number_CRs_found = 1
            pixel_sorted_index = sort_index[row1[j], col1[j], :]
            pixel_cr_mask[pixel_sorted_index[ndiffs - pixel_sat_groups - 1]] = 0  # setting largest diff to be a CR
            new_CR_found = True

            # Loop and see if there is more than one CR, setting the mask as you go
            while new_CR_found and ((ndiffs - number_CRs_found - pixel_sat_groups) > 1):
                new_CR_found = False
                # For this pixel get a new median difference excluding the number of CRs found and
                # the number of saturated groups
                pixel_med_diff = get_clipped_median(ndiffs, number_CRs_found + pixel_sat_groups,
                                                    pixel_masked_diffs, pixel_sorted_index)

                # Recalculate the noise and ratio for this pixel now that we have rejected a CR
                pixel_poisson_noise = np.sqrt(np.abs(pixel_med_diff))
                pixel_sigma = np.sqrt(pixel_poisson_noise * pixel_poisson_noise + pixel_rn2 / nframes)
                pixel_ratio = np.abs(pixel_masked_diffs - pixel_med_diff) / pixel_sigma

                # Check if largest remaining difference is above threshold
                if pixel_ratio[pixel_sorted_index[ndiffs - number_CRs_found - pixel_sat_groups - 1]] > rej_threshold:
                    new_CR_found = True
                    pixel_cr_mask[pixel_sorted_index[ndiffs - number_CRs_found - pixel_sat_groups - 1]] = 0
                    number_CRs_found += 1

            # Found all CRs for this pixel. Set CR flags in input DQ array for this pixel
            gdq[integ, 1:, row1[j], col1[j]] = \
                np.bitwise_or(gdq[integ, 1:, row1[j], col1[j]], JUMP_DET * np.invert(pixel_cr_mask))

        # Flag neighbors of pixels with detected jumps, if requested
        if flag_4_neighbors:
            cr_group, cr_row, cr_col = np.where(np.bitwise_and(gdq[integ], JUMP_DET))
            for j in range(len(cr_group)):

                # Jumps must be within a certain range to have neighbors flagged
                if ratio[cr_row[j], cr_col[j], cr_group[j] - 1] < max_jump_to_flag_neighbors and \
                    ratio[cr_row[j], cr_col[j], cr_group[j] - 1] > min_jump_to_flag_neighbors:

                    # This section saves flagged neighbors that are above or below the current range of row. If this
                    # method is running in a single process, the row above and below are not used. If it is running in
                    # multiprocessing mode, then the rows above and below need to be returned to find_jumps to use when
                    # it reconstructs the full group dq array from the slices.
                    if cr_row[j] != 0:
                        gdq[integ, cr_group[j], cr_row[j] - 1, cr_col[j]] = np.bitwise_or(
                            gdq[integ, cr_group[j], cr_row[j] - 1, cr_col[j]], JUMP_DET)
                    else:
                        row_below_gdq[integ, cr_group[j], cr_col[j]] = JUMP_DET

                    if cr_row[j] != nrows - 1:
                        gdq[integ, cr_group[j], cr_row[j] + 1, cr_col[j]] = np.bitwise_or(
                            gdq[integ, cr_group[j], cr_row[j] + 1, cr_col[j]], JUMP_DET)
                    else:
                        row_above_gdq[integ, cr_group[j], cr_col[j]] = JUMP_DET

                    # Here we are just checking that we don't flag neighbors of jumps that are off the detector.
                    if cr_col[j] != 0:
                        gdq[integ, cr_group[j], cr_row[j], cr_col[j] - 1] = np.bitwise_or(
                            gdq[integ, cr_group[j], cr_row[j], cr_col[j] - 1], JUMP_DET)
                    if cr_col[j] != ncols - 1:
                        gdq[integ, cr_group[j], cr_row[j], cr_col[j] + 1] = np.bitwise_or(
                            gdq[integ, cr_group[j], cr_row[j], cr_col[j] + 1], JUMP_DET)

    # All done
    return gdq, row_below_gdq, row_above_gdq


def get_clipped_median(num_differences, diffs_to_ignore, differences, sorted_index):
    """
    This routine will return the clipped median for the input array or pixel.
    It will ignore the input number of largest differences. At a minimum this
    is at least one plus the number of saturated values, to avoid the median being biased
    by a cosmic ray. As cosmic rays are found, the diffs_to_ignore will increase.
    """

    # ignore largest value and number of CRs found when finding new median
    # Check to see if this is a 2-D array or 1-D
    if sorted_index.ndim > 1:
        # Get the index of the median value always excluding the highest value
        # In addition, decrease the index by 1 for every two diffs_to_ignore,
        # these will be saturated values in this case
        row, col = np.indices(diffs_to_ignore.shape)
        pixel_med_index = sorted_index[row, col, (num_differences - (diffs_to_ignore[row, col] + 1))//2]
        pixel_med_diff = differences[row, col, pixel_med_index]

        # For pixels with an even number of differences the median is the mean of the two central values.
        # So we need to get the value the other central difference one lower in the sorted index that the one found
        # above.
        even_group_rows, even_group_cols = np.where((num_differences - diffs_to_ignore - 1) % 2 == 0)
        pixel_med_index2 = np.zeros_like(pixel_med_index)
        pixel_med_index2[even_group_rows, even_group_cols] = \
            sorted_index[even_group_rows, even_group_cols,
                         (num_differences - (diffs_to_ignore[even_group_rows, even_group_cols] + 3))//2]

        # Average together the two central values
        pixel_med_diff[even_group_rows, even_group_cols] = (
            pixel_med_diff[even_group_rows, even_group_cols] +
            differences[even_group_rows, even_group_cols, pixel_med_index2[even_group_rows, even_group_cols]])/2.0

    # The 1-D array case is a lot simplier.
    else:
        pixel_med_index = sorted_index[int(((num_differences - 1 - diffs_to_ignore) / 2))]
        pixel_med_diff = differences[pixel_med_index]
        if (num_differences - diffs_to_ignore - 1) % 2 == 0:  # even number of differences
            pixel_med_index2 = sorted_index[int((num_differences - 1 - diffs_to_ignore) / 2) - 1]
            pixel_med_diff = (pixel_med_diff + differences[pixel_med_index2]) / 2.0

    return pixel_med_diff
