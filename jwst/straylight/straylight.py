# Module for  applying straylight correction.
#
# The routine correct_xartifact applies a cross-artifact correction to MRS
# science slope images.  After computing this correction based on reference
# model parameters applied to the observed detector image it will be subtracted
# from the detector image.  This effectively removes unpleasant 'detector'
# PSF effects that are non-local on the sky.

import numpy as np
import logging

from stdatamodels.jwst.datamodels import dqflags
from astropy.stats import sigma_clipped_stats as scs
from .calc_xart import xart_wrapper  # c extension

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


# C version of the fitting code
def makemodel_ccode(fimg, xvec, imin, imax, lor_fwhm, lor_amp, g_fwhm, g_dx, g1_amp, g2_amp):

    fuse = fimg.copy()
    badval = np.where(fuse < 0.)
    if len(badval[0]) > 0:
        fuse[badval] = 0.
    fuse1d = fuse.ravel()

    gamma = lor_fwhm / 2.
    g_std = g_fwhm / (2 * np.sqrt(2. * np.log(2)))

    xsize, ysize = 1032, 1024

    result = xart_wrapper(imin, imax, xsize, ysize, xvec, fuse1d, gamma, lor_amp, g_std, g_dx, g1_amp, g2_amp)

    model = np.reshape(result[0], fimg.shape)

    result = None

    return model


# Python version of the fitting code
def makemodel_composite(fimg, xvec, imin, imax, lor_fwhm, lor_amp, g_fwhm, g_dx, g1_amp, g2_amp):

    model = np.zeros_like(fimg)
    model1d = model.ravel()

    fuse = fimg.copy()
    badval = np.where(fuse < 0.)
    if len(badval[0]) > 0:
        fuse[badval] = 0.
    fuse1d = fuse.ravel()

    gamma = lor_fwhm / 2.
    gstd = g_fwhm / (2 * np.sqrt(2. * np.log(2)))

    for yy in range(0, 1024):
        for ii in range(imin, imax):
            model1d[1032 * yy:1032 * (yy + 1)] += \
                (fuse1d[yy * 1032 + ii] * lor_amp[yy] * gamma[yy] * gamma[yy]) \
                / (gamma[yy] * gamma[yy] + (xvec - ii) * (xvec - ii))
            model1d[1032 * yy:1032 * (yy + 1)] += \
                (fuse1d[yy * 1032 + ii] * g1_amp[yy]
                 * np.exp(-((xvec - ii - g_dx[yy]) * (xvec - ii - g_dx[yy])) / (2 * gstd[yy] * gstd[yy])))
            model1d[1032 * yy:1032 * (yy + 1)] += \
                (fuse1d[yy * 1032 + ii] * g1_amp[yy]
                 * np.exp(-((xvec - ii + g_dx[yy]) * (xvec - ii + g_dx[yy])) / (2 * gstd[yy] * gstd[yy])))
            model1d[1032 * yy:1032 * (yy + 1)] += \
                (fuse1d[yy * 1032 + ii] * g2_amp[yy]
                 * np.exp(-((xvec - ii - 2 * g_dx[yy]) * (xvec - ii - 2 * g_dx[yy])) / (8 * gstd[yy] * gstd[yy])))
            model1d[1032 * yy:1032 * (yy + 1)] += \
                (fuse1d[yy * 1032 + ii] * g2_amp[yy]
                 * np.exp(-((xvec - ii + 2 * g_dx[yy]) * (xvec - ii + 2 * g_dx[yy])) / (8 * gstd[yy] * gstd[yy])))

    return model


def correct_xartifact(input_model, modelpars):
    """
    Corrects the MIRI MRS data for 'straylight' produced by the cross-artifact.

    Parameters
    ----------
    input_model : `~jwst.datamodels.IFUImageModel`
        Science data to be corrected.

    modelpars : FITS binary table
        Holds the reference parameters to be used to build the cross-artifact model

    Returns
    -------
    output : `~jwst.datamodels.IFUImageModel`
        Straylight-subtracted science data.

    """

    # Save some data parameters for easy use later
    nrows, ncols = input_model.data.shape

    # Create a copy of the input data array that will be modified
    # for use in the straylight calculations
    usedata = input_model.data.copy()

    # mask is same size as image - set = 1 everywhere to start
    mask = np.ones_like(usedata)

    # if there are nans remove them because they mess up the correction
    # Set them to zero in the copy of input data used for calculation
    index_inf = np.isinf(usedata).nonzero()
    usedata[index_inf] = 0.0
    index_inf = np.isnan(usedata).nonzero()
    usedata[index_inf] = 0.0
    # flag associated mask so we do not  use any
    # slice gaps that are nans, now data=0.
    mask[index_inf] = 0

    # flag pixels that should not be used for computation
    mask_dq = input_model.dq.copy()  # * mask # find DQ flags of the gap values
    all_flags = (dqflags.pixel['DO_NOT_USE'] + dqflags.pixel['NON_SCIENCE']
                 + dqflags.pixel['DEAD'] + dqflags.pixel['HOT'])
    # where are pixels set to any one of the all_flags cases
    testflags = np.bitwise_and(mask_dq, all_flags)
    # where are testflags ne 0 and mask == 1
    bad_flags = np.where(testflags != 0)
    mask[bad_flags] = 0

    # Which channel/band are we using?
    channel = input_model.meta.instrument.channel
    band = input_model.meta.instrument.band

    # Define the 'dark' region to use between the channels
    if (channel == '12'):
        xd1, xd2 = 503, 516
        yd1, yd2 = 50, 1024-50
    else:
        xd1, xd2 = 474, 507
        yd1, yd2 = 50, 1024-50

    # Estimate a first-pass pedestal correction based on sigma-clipped
    # statistics in the dark region. We need to do this here so that
    # the cross-artifact correction doesn't get applied to dark glow, but since
    # the pedestal measurement can also be biased by cross-artifact we'll estimate
    # the final pedestal correction afterwards
    pedestal_guess, _, pedestal_rms = scs(usedata[yd1:yd2, xd1:xd2])

    # Don't apply the cross artifact correction to pixels dominated by detector noise
    indx = np.where(usedata < pedestal_guess + 3 * pedestal_rms)
    mask[indx] = 0

    # Deal with normal cases only, we won't apply to cross-dichroic cases for now
    left, right = 'N/A', 'N/A'

    if channel == '12' and band == 'SHORT':
        left, right = 'ch1a_table', 'ch2a_table'
    if channel == '12' and band == 'MEDIUM':
        left, right = 'ch1b_table', 'ch2b_table'
    if channel == '12' and band == 'LONG':
        left, right = 'ch1c_table', 'ch2c_table'
    if channel == '34' and band == 'SHORT':
        left, right = 'ch4a_table', 'ch3a_table'
    if channel == '34' and band == 'MEDIUM':
        left, right = 'ch4b_table', 'ch3b_table'
    if channel == '34' and band == 'LONG':
        left, right = 'ch4c_table', 'ch3c_table'

    # Catch failure cases with a log warning
    if left == 'N/A' or right == 'N/A':
        log.info("Warning: no parameters found for channel = " + str(channel) + " band = " + str(band))

    xvec = (np.arange(ncols)).astype(float)
    left_model = np.zeros_like(usedata)
    right_model = np.zeros_like(usedata)

    # Left-half of detector
    try:
        param = modelpars[left]
        log.info("Found parameters for left detector half, applying Cross-Artifact correction.")
        istart, istop = 0, 516
        # Subtract off intial guess at pedestal first
        fimg = (usedata - pedestal_guess) * mask
        left_model = makemodel_ccode(fimg, xvec, istart, istop, param['LOR_FWHM'],
                                     param['LOR_SCALE'], param['GAU_FWHM'],
                                     param['GAU_XOFF'], param['GAU_SCALE1'], param['GAU_SCALE2'])
    except Exception:
        left_model[:, :] = 0
        log.info("No parameters for left detector half, not applying Cross-Artifact correction.")

    # Right-half of detector
    try:
        param = modelpars[right]
        log.info("Found parameters for right detector half, applying Cross-Artifact correction.")
        istart, istop = 516, 1024
        fimg = (usedata - pedestal_guess) * mask
        right_model = makemodel_ccode(fimg, xvec, istart, istop, param['LOR_FWHM'],
                                      param['LOR_SCALE'], param['GAU_FWHM'],
                                      param['GAU_XOFF'], param['GAU_SCALE1'], param['GAU_SCALE2'])
    except Exception:
        right_model[:, :] = 0
        log.info("No parameters for right detector half, not applying Cross-Artifact correction.")

    model = left_model + right_model
    # remove the straylight correction for the reference pixels
    model[:, 1028:1032] = 0.0
    model[:, 0:4] = 0.0

    # Create output as a copy of the real data prior to replacement of NaNs with zeros
    output = input_model.copy()
    # Subtract the model from the data
    output.data = output.data - model
    usedata = usedata - model

    # Now measure and remove the pedestal dark count rate measured between the channels
    # Embed in a try/except block to catch unusual failures
    try:
        _, themed, therms = scs(usedata[yd1:yd2, xd1:xd2])
        pedestal = np.zeros_like(output.data) + themed
        # remove the pedestal correction for the reference pixels
        pedestal[:, 1028:1032] = 0.0
        pedestal[:, 0:4] = 0.0

        output.data = output.data - pedestal
        log.info("Derived pedestal correction " + str(themed) + " DN/s")
    except Exception:
        log.info("Straylight pedestal correction failed.")

    # Delete our temporary working copy of the data
    del usedata

    log.info("Cross-artifact model complete.")
    return output
