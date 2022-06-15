#  Module for applying the LG-PLUS algorithm to an AMI exposure
import logging
import numpy as np

from .find_affine2d_parameters import find_rotation
from . import instrument_data
from . import nrm_core
from .utils import img_median_replace

from astropy import units as u

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def apply_LG_plus(input_model, filter_model, oversample, rotation,
                  psf_offset, rotsearch_parameters):
    """
    Short Summary
    -------------
    Applies the image plane algorithm to an AMI image

    Parameters
    ----------
    input_model : data model object
        AMI science image to be analyzed

    filter_model : filter model object
        filter throughput reference data

    oversample : integer
        Oversampling factor

    rotation : float (degrees)
        Initial guess at rotation of science image relative to model

    Returns
    -------
    output_model : Fringe model object
        Fringe analysis data

    """
    # Create copy of input_model to avoid overwriting input
    input_copy = input_model.copy()

    # If the input data were taken in full-frame mode, extract a region
    # equivalent to the SUB80 subarray mode to make execution time acceptable.
    if input_model.meta.subarray.name.upper() == 'FULL':
        log.info("Extracting 80x80 subarray from full-frame data")
        xstart = 1045
        ystart = 1
        xsize = 80
        ysize = 80
        xstop = xstart + xsize - 1
        ystop = ystart + ysize - 1
        input_copy.data = input_copy.data[ystart - 1:ystop, xstart - 1:xstop]
        input_copy.dq = input_copy.dq[ystart - 1:ystop, xstart - 1:xstop]
        input_copy.err = input_copy.err[ystart - 1:ystop, xstart - 1:xstop]

    # Replace NaN's and DO_NOT_USE pixels in the input image
    # with median of surrounding pixel values in a 3x3 box
    box_size = 3
    input_copy = img_median_replace(input_copy, box_size)

    data = input_copy.data
    dim = data.shape[1]

    # Set transformation parameters:
    #   mx, my: dimensionless magnifications
    #   sx, sy: dimensionless shears
    #   x0, y0: offsets in pupil space
    mx = 1.0
    my = 1.0
    sx = 0.0
    sy = 0.0
    xo = 0.0
    yo = 0.0

    psf_offset_ff = None

    lamc = 4.3e-6
    oversample = 11
    bandpass = np.array([(1.0, lamc), ])
    pixelscale_as = 0.0656
    arcsec2rad = u.arcsec.to(u.rad)
    PIXELSCALE_r = pixelscale_as * arcsec2rad
    holeshape = 'hex'
    filt = "F430M"
    rotsearch_d = np.append(np.arange(rotsearch_parameters[0], rotsearch_parameters[1], rotsearch_parameters[2]),
                            rotsearch_parameters[1])

    log.info(f'Initial values to use for rotation search {rotsearch_d}')

    affine2d = find_rotation(data[:, :], psf_offset, rotsearch_d,
                             mx, my, sx, sy, xo, yo,
                             PIXELSCALE_r, dim, bandpass, oversample, holeshape)

    niriss = instrument_data.NIRISS(filt, bandpass=bandpass, affine2d=affine2d)

    ff_t = nrm_core.FringeFitter(niriss, psf_offset_ff=psf_offset_ff,
                                 oversample=oversample)

    output_model = ff_t.fit_fringes_all(input_copy)

    # Copy header keywords from input to output
    output_model.update(input_model, only="PRIMARY")

    return output_model
