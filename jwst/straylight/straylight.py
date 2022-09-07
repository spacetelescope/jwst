# Module for  applying straylight correction.
#
# The routine correct_xartifact applies a cross-artifact correction to MRS
# science slope images.  After computing this correction based on reference
# model parameters applied to the observed detector image it will be subtracted
# from the detector image.  This effectively removes unpleasant 'detector'
# PSF effects that are non-local on the sky.

import numpy as np
import logging
from ..datamodels import dqflags
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

    # Create output as a copy of the input science data model
    output = input_model.copy()  # this is used in algorithm to
    # find the straylight correction.

    # mask is same size as image - set = 1 everywhere to start
    mask = np.ones_like(output.data)

    # if there are nans remove them because they mess up the correction
    index_inf = np.isinf(output.data).nonzero()
    output.data[index_inf] = 0.0
    index_inf = np.isnan(output.data).nonzero()
    output.data[index_inf] = 0.0
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
    left_model = np.zeros_like(output.data)
    right_model = np.zeros_like(output.data)

    # Left-half of detector
    try:
        param = modelpars[left]
        log.info("Found parameters for left detector half, applying Cross-Artifact correction.")
        istart, istop = 0, 516
        fimg = output.data * mask
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
        fimg = output.data * mask
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
    output.data = output.data - model

    # Remove any remaining pedestal stray light values based on inter-channel region
    # chosen based on analysis of multiple bands of data from various programs
    if (channel == '12'):
        pedestal = np.nanmedian(output.data[280:740, 503:516])
    else:
        pedestal = np.nanmedian(output.data[280:740, 474:507])
    try:
        output.data = output.data - pedestal
        log.info("Derived pedestal correction " + str(pedestal) + " DN/s")
    except Exception:
        log.info("Straylight pedestal correction failed.")

    log.info("Cross-artifact model complete.")
    return output
