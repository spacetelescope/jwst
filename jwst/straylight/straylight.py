# Module for  applying straylight correction.
#
# The routine correct_xartifact applies a cross-artifact correction to MRS
# science slope images.  This new routine entirely replaces previous routines
# which treated the cross-artifact as if it were traditional stray light
# rather than a manifestation of internal reflections within the detector
# substrate (i.e., a 'detector' PSF).

import numpy as np
import logging
from ..datamodels import dqflags
from astropy.io import fits
import pdb
from .straylight_xartifact import xartifact_wrapper  # c extension

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def makemodel_ccode(fimg,xvec,imin,imax,lor_fwhm,lor_amp,g_fwhm,g_dx,g1_amp,g2_amp):
    model = np.zeros_like(fimg)
    model1d = model.ravel()

    fuse = fimg.copy()
    badval = np.where(fuse < 0.)
    if (len(badval[0]) > 0):
        fuse[badval] = 0.
    fuse1d=fuse.ravel()

    gamma = lor_fwhm / 2.
    gstd = g_fwhm / (2 * np.sqrt(2. * np.log(2)))

    xsize, ysize = 1032, 1024
    print('Entering C code')
    #pdb.set_trace()
    result = xartifact_wrapper(fuse1d,xvec,model1d,xsize,ysize,imin,imax,gamma,lor_amp,
                               gstd,g_dx,g1_amp,g2_amp)
    result = None
    pdb.set_trace()

    return model

def makemodel_composite(fimg,xvec,imin,imax,lor_fwhm,lor_amp,g_fwhm,g_dx,g1_amp,g2_amp):
    model = np.zeros_like(fimg)
    model1d = model.ravel()

    fuse = fimg.copy()
    badval = np.where(fuse < 0.)
    if (len(badval[0]) > 0):
        fuse[badval] = 0.
    fuse1d=fuse.ravel()

    gamma = lor_fwhm / 2.
    gstd = g_fwhm / (2 * np.sqrt(2. * np.log(2)))

    for yy in range(0, 1024):
        for ii in range(imin, imax):
            model1d[1032 * yy:1032 * (yy + 1)] += (fuse1d[yy * 1032 + ii] * lor_amp[yy] * gamma[yy] * gamma[yy]) / (
                        gamma[yy] * gamma[yy] + (xvec - ii) * (xvec - ii))
            model1d[1032 * yy:1032 * (yy + 1)] += (fuse1d[yy * 1032 + ii] * g1_amp[yy] * np.exp(
                -((xvec - ii - g_dx[yy]) * (xvec - ii - g_dx[yy])) / (2 * gstd[yy] * gstd[yy])))
            model1d[1032 * yy:1032 * (yy + 1)] += (fuse1d[yy * 1032 + ii] * g1_amp[yy] * np.exp(
                -((xvec - ii + g_dx[yy]) * (xvec - ii + g_dx[yy])) / (2 * gstd[yy] * gstd[yy])))
            model1d[1032 * yy:1032 * (yy + 1)] += (fuse1d[yy * 1032 + ii] * g2_amp[yy] * np.exp(
                -((xvec - ii - 2 * g_dx[yy]) * (xvec - ii - 2 * g_dx[yy])) / (8 * gstd[yy] * gstd[yy])))
            model1d[1032 * yy:1032 * (yy + 1)] += (fuse1d[yy * 1032 + ii] * g2_amp[yy] * np.exp(
                -((xvec - ii + 2 * g_dx[yy]) * (xvec - ii + 2 * g_dx[yy])) / (8 * gstd[yy] * gstd[yy])))

    return model

def correct_xartifact(input_model, modelpars):
    """
    Corrects the MIRI MRS data for 'straylight' produced by the cross-artifact.

    Parameters
    ----------
    input_model: `~jwst.datamodels.IFUImageModel`
        Science data to be corrected.

    modelpars: FITS binary table
        Holds the reference parameters to be used to build the cross-artifact model

    Returns
    -------
    output: `~jwst.datamodels.IFUImageModel`
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

    # TODO: Deal with crossed-band
    if ((channel == '12') & (band == 'SHORT')):
        left, right = 'CH1A', 'CH2A'
    if ((channel == '12') & (band == 'MEDIUM')):
        left, right = 'CH1B', 'CH2B'
    if ((channel == '12') & (band == 'LONG')):
        left, right = 'CH1C', 'CH2C'
    if ((channel == '34') & (band == 'SHORT')):
        left, right = 'CH4A', 'CH3A'
    if ((channel == '34') & (band == 'MEDIUM')):
        left, right = 'CH4B', 'CH3B'
    if ((channel == '34') & (band == 'LONG')):
        left, right = 'CH4C', 'CH3C'

    xvec = (np.arange(ncols)).astype(float)
    left_model = np.zeros_like(output.data)
    right_model = np.zeros_like(output.data)

    # Left-half of detector
    try:
        param = modelpars[left].data
        log.info("Found parameters for left detector half, applying Cross-Artifact correction.")
        istart, istop = 0, 516
        fimg = output.data * mask
        left_model = makemodel_ccode(fimg, xvec, istart, istop, param['LOR_FWHM'],
                                         param['LOR_SCALE'], param['GAU_FWHM'],
                                         param['GAU_XOFF'], param['GAU_SCALE1'],
                                         param['GAU_SCALE2'])
    except:
        log.info("No parameters for left detector half, not applying Cross-Artifact correction.")

    # Right-half of detector
    try:
        param = modelpars[right].data
        log.info("Found parameters for left detector half, applying Cross-Artifact correction.")
        istart, istop = 516, 1024
        fimg = output.data * mask
        right_model = makemodel_composite(fimg, xvec, istart, istop, param['LOR_FWHM'],
                                         param['LOR_SCALE'], param['GAU_FWHM'],
                                         param['GAU_XOFF'], param['GAU_SCALE1'],
                                         param['GAU_SCALE2'])
    except:
        log.info("No parameters for right detector half, not applying Cross-Artifact correction.")

    model = left_model + right_model
    # remove the straylight correction for the reference pixels
    model[:, 1028:1032] = 0.0
    model[:, 0:4] = 0.0
    output.data = output.data - model

    # Test
    hdu=fits.PrimaryHDU(model)
    hdu.writeto('pipemodel.fits',overwrite=True)

    log.info("Cross-artifact model complete.")
    return output
