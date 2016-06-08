from __future__ import absolute_import

import time
import logging

import numpy as np
from . import twopoint_difference as twopt
from . import yintercept as yint

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

def detect_jumps (input_model, gain_model, readnoise_model,
                  rejection_threshold, do_yint, signal_threshold):
    """
    This is the high-level controlling routine for the jump detection process.
    It loads and sets the various input data and parameters needed by each of
    the individual detection methods and then calls the detection methods in
    turn.

    Note that the detection methods are currently setup on the assumption
    that the input science and error data arrays will be in units of
    electrons, hence this routine scales those input arrays by the detector
    gain. The methods assume that the read noise values will also be in units
    of electrons.

    The gain is applied to the science data and error arrays using the
    appropriate instrument- and detector-dependent values for each pixel of an
    image.  Also, a 2-dimensional read noise array with appropriate values for
    each pixel is passed to the detection methods.
    """

    # Load the data arrays that we need from the input model
    output_model = input_model.copy()
    data = input_model.data
    err  = input_model.err
    gdq  = input_model.groupdq

    ngroups = data.shape[1]

    # Get subarray limits from metadata of input model
    xstart = input_model.meta.subarray.xstart
    xsize  = input_model.meta.subarray.xsize
    xstop  = xstart + xsize - 1
    ystart = input_model.meta.subarray.ystart
    ysize  = input_model.meta.subarray.ysize
    ystop  = ystart + ysize - 1

    # Get 2D gain and read noise values from their respective models
    gain_2d = gain_model.data[ystart-1:ystop,xstart-1:xstop]

    if (readnoise_model.meta.subarray.xstart==xstart and
        readnoise_model.meta.subarray.xsize==xsize   and
        readnoise_model.meta.subarray.ystart==ystart and
        readnoise_model.meta.subarray.ysize==ysize):

        log.debug('Readnoise subarray matches science data')
        readnoise_2d = readnoise_model.data
    else:
        log.debug('Extracting readnoise subarray to match science data')
        readnoise_2d = readnoise_model.data[ystart-1:ystop,xstart-1:xstop]

    # Apply gain to the SCI and ERR arrays so they're in units of electrons
    data *= gain_2d
    err  *= gain_2d

    # Apply the 2-point difference method as a first pass
    log.info('Executing two-point difference method')
    start = time.time()
    median_slopes = twopt.find_CRs( data, gdq, readnoise_2d,
                                    rejection_threshold)
    elapsed = time.time() - start
    log.debug('Elapsed time = %g sec' %elapsed)

    # Apply the y-intercept method as a second pass, if requested
    if do_yint:

        # Set up the ramp time array for the y-intercept method
        group_time = output_model.meta.exposure.group_integration_time
        times = np.array([(k+1)*group_time for k in range(ngroups)])
        median_slopes /= group_time

        # Now apply the y-intercept method
        log.info('Executing yintercept method')
        start = time.time()
        yint.find_CRs( data, err, gdq, times, readnoise_2d,
                       rejection_threshold, signal_threshold, median_slopes)
        elapsed = time.time() - start
        log.debug('Elapsed time = %g sec' %elapsed)

    # Update the DQ array of the output model with the jump detection results
    output_model.groupdq = gdq

    return output_model



