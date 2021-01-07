import time
import logging

import numpy as np
from ..datamodels import dqflags
from ..lib import reffile_utils
from . import twopoint_difference as twopt
import multiprocessing

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

def detect_jumps (input_model, gain_model, readnoise_model,
                  rejection_threshold, max_cores,
                  max_jump_to_flag_neighbors, min_jump_to_flag_neighbors,
                  flag_4_neighbors):
    """
    This is the high-level controlling routine for the jump detection process.
    It loads and sets the various input data and parameters needed by each of
    the individual detection methods and then calls the detection methods in
    turn.

    Note that the detection methods are currently setup on the assumption
    that the input science and error data arrays will be in units of
    electrons, hence this routine scales those input arrays by the detector
    gain. The methods assume that the read noise values will be in units
    of DN.

    The gain is applied to the science data and error arrays using the
    appropriate instrument- and detector-dependent values for each pixel of an
    image.  Also, a 2-dimensional read noise array with appropriate values for
    each pixel is passed to the detection methods.
    """
    if max_cores is None:
        numslices = 1
    else:
        num_cores = multiprocessing.cpu_count()
        log.info("Found %d possible cores to use for jump detection " % num_cores)
        if max_cores == 'quarter':
            numslices = num_cores // 4 or 1
        elif max_cores == 'half':
            numslices = num_cores // 2 or 1
        elif max_cores == 'all':
            numslices = num_cores
        else:
            numslices = 1

    # Load the data arrays that we need from the input model
    output_model = input_model.copy()
    data = input_model.data
    err  = input_model.err
    gdq  = input_model.groupdq
    pdq  = input_model.pixeldq


    # Get 2D gain and read noise values from their respective models
    if reffile_utils.ref_matches_sci(input_model, gain_model):
        gain_2d = gain_model.data
    else:
        log.info('Extracting gain subarray to match science data')
        gain_2d = reffile_utils.get_subarray_data(input_model, gain_model)

    if reffile_utils.ref_matches_sci(input_model, readnoise_model):
        readnoise_2d = readnoise_model.data
    else:
        log.info('Extracting readnoise subarray to match science data')
        readnoise_2d = reffile_utils.get_subarray_data(input_model, readnoise_model)

    # Flag the pixeldq where the gain is <=0 or NaN so they will be ignored
    wh_g = np.where( gain_2d <= 0.)
    if len(wh_g[0] > 0):
        pdq[wh_g] = np.bitwise_or( pdq[wh_g], dqflags.pixel['NO_GAIN_VALUE'] )
        pdq[wh_g] = np.bitwise_or( pdq[wh_g], dqflags.pixel['DO_NOT_USE'] )

    wh_g = np.where( np.isnan( gain_2d ))
    if len(wh_g[0] > 0):
        pdq[wh_g] = np.bitwise_or( pdq[wh_g], dqflags.pixel['NO_GAIN_VALUE'] )
        pdq[wh_g] = np.bitwise_or( pdq[wh_g], dqflags.pixel['DO_NOT_USE'] )

    # Apply gain to the SCI, ERR, and readnoise arrays so they're in units
    # of electrons

    data *= gain_2d
    err  *= gain_2d
    readnoise_2d *= gain_2d

    # Apply the 2-point difference method as a first pass
    log.info('Executing two-point difference method')
    start = time.time()
    nrows = data.shape[-2]
    ncols = data.shape[-1]
    num_groups = data.shape[1]
    num_ints = data.shape[0]
    frames_per_group = input_model.meta.exposure.nframes
    row_above_gdq = np.zeros((num_ints, num_groups, ncols), dtype=np.uint8)
    previous_row_above_gdq = np.zeros((num_ints, num_groups, ncols), dtype=np.uint8)
    row_below_gdq = np.zeros((num_ints, num_groups, ncols), dtype=np.uint8)

    yincrement = int(nrows / numslices)
    slices = []
    # Slice up data, gdq, readnoise_2d into slices
    # Each element of slices is a tuple of
    # (data, gdq, readnoise_2d, rejection_threshold, nframes)
    for i in range(numslices - 1):
        slices.insert(i, (data[:, :, i * yincrement:(i + 1) * yincrement, :],
                          gdq[:, :, i * yincrement:(i + 1) * yincrement, :],
                          readnoise_2d[i * yincrement:(i + 1) * yincrement, :],
                          rejection_threshold, frames_per_group, flag_4_neighbors,
                          max_jump_to_flag_neighbors, min_jump_to_flag_neighbors))
    # last slice get the rest
    slices.insert(numslices - 1, (data[:, :, (numslices - 1) * yincrement:nrows, :],
                                 gdq[:, :, (numslices - 1) * yincrement:nrows, :],
                                 readnoise_2d[(numslices - 1) * yincrement:nrows, :],
                                 rejection_threshold, frames_per_group, flag_4_neighbors,
                                 max_jump_to_flag_neighbors, min_jump_to_flag_neighbors))
    if numslices == 1:
        gdq, row_below_dq, row_above_dq = twopt.find_crs(data, gdq, readnoise_2d, rejection_threshold,
                                                                        frames_per_group, flag_4_neighbors,
                                                                        max_jump_to_flag_neighbors,
                                                                        min_jump_to_flag_neighbors)
        elapsed = time.time() - start
    else:
        log.info("Creating %d processes for jump detection " % numslices)
        pool = multiprocessing.Pool(processes=numslices)
        # Starts each slice in it's own process. Starmap allows more than one parameter to be passed.
        real_result = pool.starmap(twopt.find_crs, slices)
        pool.close()
        pool.join()
        k = 0
        # Reconstruct gdq, the row_above_gdq, and the row_below_gdq from the slice result
        for resultslice in real_result:

            if len(real_result) == k + 1:  # last result
                gdq[:, :, k * yincrement:nrows, :] = resultslice[0]
            else:
                gdq[:, :, k * yincrement:(k + 1) * yincrement, :] = resultslice[0]
            row_below_gdq[:, :, :] = resultslice[1]
            row_above_gdq[:, :, :] = resultslice[2]
            if k != 0: # for all but the first slice, flag any CR neighbors in the top row of the previous slice and
                # flag any neighbors in the bottom row of this slice saved from the top of the previous slice
                gdq[:, :, k * yincrement - 1, :] = np.bitwise_or(gdq[:, :, k * yincrement - 1, :],
                                                                     row_below_gdq[:, :, :])
                gdq[:, :, k * yincrement, :] = np.bitwise_or(gdq[:, :, k * yincrement, :],
                                                                     previous_row_above_gdq[:, :, :])
            # save the neighbors to be flagged that will be in the next slice
            previous_row_above_gdq = row_above_gdq.copy()
            k += 1
        elapsed = time.time() - start

    elapsed = time.time() - start
    log.info('Total elapsed time = %g sec' % elapsed)

    # Update the DQ arrays of the output model with the jump detection results
    output_model.groupdq = gdq
    output_model.pixeldq = pdq

    return output_model
