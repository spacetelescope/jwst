'''
Two-Point Difference method for finding outliers in a 3-d ramp data cube.

The scheme used in this variation of the method uses numpy array methods
to compute first-differences and find the max outlier in each pixel while
still working in the full 3-d data array. This makes detection of the first
outlier very fast. We then iterate pixel-by-pixel over only those pixels
that are already known to contain an outlier, to look for any additional
outliers and set the appropriate DQ mask for all outliers in the pixel.

This is MUCH faster than doing all the work on a pixel-by-pixel basis.
'''

import logging
import numpy as np
from ..datamodels import dqflags

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

HUGE_NUM = np.finfo(np.float32).max

def find_CRs(data, gdq, read_noise, rej_threshold):

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
    read_noise_2 = read_noise * read_noise

    # Reset saturated values in input data array to NaN, so they don't get
    # used in any of the subsequent calculations
    data[gdq == dqflags.group['SATURATED']] = np.NaN

    # Loop over multiple integrations
    for integration in range(nints):

        # Roll the ngroups axis of data arrays to the end, to make
        # memory access to the values for a given pixel faster
        rdata = np.rollaxis(data[integration], 0, 3)
        rdq = np.rollaxis(gdq[integration], 0, 3)

        # Compute first differences of adjacent groups up the ramp
        first_diffs = np.diff(rdata, axis=2)

        # Compute median of the first differences for all pixels
        med_diffs = np.nanmedian(first_diffs, axis=2)

        # Zero-out results for pixels that have NaN's in all groups so they
        # don't cause trouble in later calculations
        nans = np.where(np.isnan(med_diffs))
        med_diffs[nans] = 0.
        first_diffs[nans] = 0.

        # Save initial estimate of the median slope for all pixels
        median_slopes[integration] = med_diffs

        # Compute uncertainties as the quadrature sum of the poisson noise
        # in the first difference signal and read noise. Because the first
        # differences can be biased by CRs/jumps, we use the median signal
        # for computing the poisson noise.
        poisson_noise = np.sqrt(np.abs(med_diffs))
        sigma = np.sqrt(poisson_noise * poisson_noise + 2 * read_noise_2)

        # Reset sigma to exclude pixels with both readnoise and signal=0
        wh_sig = np.where(sigma == 0.)
        if len(wh_sig[0] > 0):
            log.debug('Twopt found %d pixels with sigma=0' % (len(wh_sig[0])))
            log.debug(' which will be reset so that no jump will be detected')
            sigma[wh_sig] = HUGE_NUM

        # Compute distance of each sample from the median in units of sigma;
        # note that the use of "abs" means we'll detect both positive and
        # negative outliers
        ratio = np.abs(first_diffs - med_diffs[:, :, np.newaxis]) / sigma[:, :, np.newaxis]

        # Find the group index of the max outlier in each pixel
        max_index = np.nanargmax(ratio, axis=2)

        # Get indices of outliers that are above the rejection threshold
        r, c = np.indices(max_index.shape)
        r1, c1 = np.where(ratio[r, c, max_index] > rej_threshold)
        log.debug('Twopt found %d pixels with at least one CR' % (len(r1)))

        # Loop over pixels that have an outlier, checking to see if they have
        # more than 1 outlier
        for j in range(len(r1)):

            # Row/col indexes of current pixel with an outlier
            row, col = r1[j], c1[j]

            # Copy first_diff group values to an array with a mask
            masked_diffs = first_diffs[row, col]
            rn2 = read_noise_2[row, col]

            # Create a saturation mask based on NaN's in the first_diffs
            sat_mask = np.isfinite(masked_diffs)

            # Create a CR mask and initialize with the max outlier
            cr_mask = np.ones(masked_diffs.shape, dtype=bool)
            cr_mask[max_index[row, col]] = 0

            # Now iteratively search for and reject additional outliers
            iter = 1
            while iter:

                # Recompute the masked median, noise, and ratios for this pixel
                med = np.median(masked_diffs[cr_mask * sat_mask])
                poisson_noise = np.sqrt(np.abs(med))
                sigma = np.sqrt(poisson_noise * poisson_noise + 2 * rn2)
                ratio = np.abs(masked_diffs - med) / sigma

                # Get a list of group indexes sorted from largest to smallest
                # deviation from the median
                sortindx = np.argsort(ratio)[::-1]

                # Check through the list to see if any qualify as an outlier
                for i in sortindx:

                    # If already masked, continue to next index
                    if not cr_mask[i] * sat_mask[i]:
                        continue

                    # If above threshold, set a CR mask and iterate on this pixel
                    elif ratio[i] > rej_threshold:
                        cr_mask[i] = 0
                        iter = 1
                        break

                    # If not above threshold, we're done with this pixel
                    else:
                        iter = 0
                        break

            # Set CR flags in input DQ array for this pixel
            gdq[integration, 1:, row, col] = np.bitwise_or \
                              (gdq[integration, 1:, row, col],
                               dqflags.group['JUMP_DET'] * np.invert(cr_mask))

            # Save the CR-cleaned median slope for this pixel
            median_slopes[integration, row, col] = med

        # Next pixel with an outlier
    # Next integration

    return median_slopes
