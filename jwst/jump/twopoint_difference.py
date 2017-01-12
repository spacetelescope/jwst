"""
Two-Point Difference method for finding outliers in a 3-d ramp data cube.

The scheme used in this variation of the method uses numpy array methods
to compute first-differences and find the max outlier in each pixel while
still working in the full 3-d data array. This makes detection of the first
outlier very fast. We then iterate pixel-by-pixel over only those pixels
that are already known to contain an outlier, to look for any additional
outliers and set the appropriate DQ mask for all outliers in the pixel.

This is MUCH faster than doing all the work on a pixel-by-pixel basis.
"""

import logging
import numpy as np

from ..datamodels import dqflags
from .. import datamodels

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

HUGE_NUM = np.finfo(np.float32).max

def find_CRs(data, gdq, read_noise, rej_threshold, nframes):

    """
    Find CRs/Jumps in each integration within the input data array.

    The input data array is assumed to be in units of electrons, i.e. already
    multiplied by the gain. We also assume that the read noise is in units of
    electrons.
    """
    # Get data characteristics
    (nints, ngroups, nrows, ncols) = data.shape
    # Create array for output median slope images
    median_slopes = np.zeros((nints, nrows, ncols), dtype=np.float32)

    # Square the read noise values, for use later
    read_noise_2 = read_noise*read_noise

    # Reset saturated values in input data array to NaN, so they don't get
    # used in any of the subsequent calculations
    data[gdq == dqflags.group['SATURATED']] = np.NaN

    # Loop over multiple integrations
    for integration in range(nints):

        # Roll the ngroups axis of data arrays to the end, to make
        # memory access to the values for a given pixel faster
        rdata = np.rollaxis(data[integration], 0, 3)

        # Compute first differences of adjacent groups up the ramp
        first_diffs = np.diff(rdata, axis=2)
        nans = np.where(np.isnan(first_diffs))
        first_diffs[nans] = 100000.

        positive_first_diffs = np.abs(first_diffs)
        diffsum= positive_first_diffs.sum
        #make all the first diffs for saturated groups be equal to 100,000 to put them above the good values in the
        #sorted index
        matching_array = np.ones(shape=(positive_first_diffs.shape[0],
                                        positive_first_diffs.shape[1],
                                        positive_first_diffs.shape[2]))*100000
        sat_groups =(positive_first_diffs == matching_array)
        number_sat_groups = (sat_groups*1).sum(axis=2)
        ndiffs = ngroups - 1
        sort_index = np.argsort(positive_first_diffs)
        med_diffs = return_clipped_median(ndiffs, number_sat_groups, positive_first_diffs, sort_index)

        # Save initial estimate of the median slope for all pixels
        median_slopes[integration] = med_diffs

        # Compute uncertainties as the quadrature sum of the poisson noise
        # in the first difference signal and read noise. Because the first
        # differences can be biased by CRs/jumps, we use the median signal
        # for computing the poisson noise. Here sigma correctly has the
        # read noise taking into account the fact that multiple frames were
        # averaged into each group./
        poisson_noise = np.sqrt(np.abs(med_diffs))
        sigma = np.sqrt(poisson_noise*poisson_noise + read_noise_2/nframes)

        # Reset sigma to exclude pixels with both readnoise and signal=0
        wh_sig = np.where(sigma == 0.)
        if len(wh_sig[0] > 0):
            log.debug('Twopt found %d pixels with sigma=0' %(len(wh_sig[0])))
            log.debug(' which will be reset so that no jump will be detected')
            sigma[wh_sig] = HUGE_NUM

        # Compute distance of each sample from the median in units of sigma;
        # note that the use of "abs" means we'll detect both positive and
        # negative outliers
        ratio = np.abs(first_diffs - med_diffs[:,:,np.newaxis])/sigma[:,:,np.newaxis]

        # Find the group index of the max outlier in each pixel
        # (the RHS is identical to the previous versions:
        #  max_index = np.nanargmax (ratio, axis=2)
        max_index1 = sort_index[:,:,ngroups-2]

        # Get indices of highest values (may be outliers) that are above the
        # rejection threshold
        r, c = np.indices(max_index1.shape)
        row1, col1 = np.where(ratio[r, c, max_index1 - number_sat_groups] > rej_threshold)
        log.debug('From highest outlier Twopt found %d pixels with at least one CR' % (len(row1)))
        number_pixels_with_cr = len(row1)
        for j in range(number_pixels_with_cr):
            pixel_masked_diffs = first_diffs[row1[j], col1[j]]
            pixel_rn2 = read_noise_2[row1[j], col1[j]]
            pixel_sat_groups = number_sat_groups[row1[j], col1[j]]
            sorted_index_of_cr = max_index1[row1[j],col1[j]] - pixel_sat_groups
            # Create a CR mask and set 1st CR to be found
            # cr_mask=0 designates a CR
            pixel_cr_mask = np.ones(pixel_masked_diffs.shape, dtype=bool)
            number_CRs_found = 1
            pixel_sorted_index = sort_index[row1[j], col1[j], :]
            pixel_cr_mask[pixel_sorted_index[ndiffs - pixel_sat_groups - 1]] = 0
            new_CR_found = True

            #loop over all the found CRs and see if there is more than one CR setting the mask as you go
            while new_CR_found and ((ndiffs - number_CRs_found - pixel_sat_groups) > 1):
                new_CR_found = False
                pixel_med_diff = return_clipped_median(ndiffs, number_CRs_found + pixel_sat_groups,
                                                       pixel_masked_diffs, pixel_sorted_index)
                poisson_noise = np.sqrt(np.abs(pixel_med_diff))
                sigma = np.sqrt(poisson_noise * poisson_noise + pixel_rn2 / nframes)
                ratio = np.abs(pixel_masked_diffs - pixel_med_diff) / sigma
                pixel_sorted_ratio = ratio[pixel_sorted_index[:]]

                #check if largest remaining difference is above threshold
                if ratio[pixel_sorted_index[ndiffs - number_CRs_found - pixel_sat_groups - 1]] > rej_threshold:
                    new_CR_found = True
                    pixel_cr_mask[pixel_sorted_index[ndiffs - number_CRs_found-pixel_sat_groups-1]]=0
                    number_CRs_found += 1
            # Found all CRs. Set CR flags in input DQ array for this pixel
            gdq[integration, 1:, row1[j], col1[j]] = np.bitwise_or \
                (gdq[integration, 1:, row1[j], col1[j]],
                 dqflags.group['JUMP_DET'] * np.invert(pixel_cr_mask))

            # Save the CR-cleaned median slope for this pixel
            if not new_CR_found: # the loop ran at least one time
                median_slopes[integration, row1[j], col1[j]] = pixel_med_diff
        # Next pixel with an outlier (j loop)
    # Next integration (integration loop)
    return median_slopes


#This routine will return the clipped median for the input array or pixel. It will ignore the input number of largest
#differences. At a minimum this is at least one plus the number of saturated values to avoid the median being biased
# by a cosmic ray. As cosmic rays are found the diffs_to_ignore will increase.
def return_clipped_median(num_differences, diffs_to_ignore, differences, sorted_index):
    # ignore largest value and number of CRs found when finding new median

    if sorted_index.ndim > 1:
        # always exclude the highest value
        pixel_med_index = sorted_index[:, :, int((num_differences - 1) / 2)] # always exclude the highest value
        row, col = np.indices(pixel_med_index.shape)
        # in addition decrease the index by 1 for every two diffs_to_ignore, these will be saturated values in this case
        pixel_med_diff = differences[row, col, pixel_med_index - ((diffs_to_ignore) / 2).astype(int)]
        if (num_differences - 1) % 2 == 0:  # even
            pixel_med_index2 = sorted_index[:, :, int((num_differences - 1) / 2 ) - 1]
            pixel_med_diff = (pixel_med_diff + differences[row, col, pixel_med_index2- ((diffs_to_ignore) / 2).astype(int)]) / 2.0
    else:
        pixel_med_index = sorted_index[int(((num_differences - 1  - diffs_to_ignore )/ 2))]
        pixel_med_diff = differences[pixel_med_index]
        if (num_differences - diffs_to_ignore - 1) % 2 == 0:  # even
            pixel_med_index2 = sorted_index[int((num_differences - 1 - diffs_to_ignore) / 2) - 1]
            pixel_med_diff = (pixel_med_diff + differences[pixel_med_index2]) / 2.0
    return pixel_med_diff
