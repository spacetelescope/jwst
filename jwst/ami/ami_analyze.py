#
#  Module for applying the LG-PLUS algorithm to an AMI exposure
#

import logging
import numpy as np

from . import find_affine2d_parameters as FAP
from . import InstrumentData
from . import nrm_core

from astropy import units as u

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def apply_LG_plus(input_model, filter_model, oversample, rotation):
    """
    Short Summary
    -------------
    Applies the image plane algorithm to an AMI image

    Parameters
    ----------
    input_model: data model object
        AMI science image to be analyzed

    filter_model: filter model object
        filter throughput reference data

    oversample: integer
        Oversampling factor

    rotation: float (degrees)
        Initial guess at rotation of science image relative to model

    Returns
    -------
    output_model: Fringe model object
        Fringe analysis data

    """
    data = input_model.data
    dim = data.shape[1]

    # Set transformation parameters:
    #   mx, my: dimensionless magnifications
    #   sx, sy: dimensionless shears
    #   x0, y0: offsets in pupil space
    mx,my,sx,sy,xo,yo, = (1.0,1.0, 0.0,0.0, 0.0,0.0) 

    psf_offset_find_rotation = (0.0,0.0)
    psf_offset_ff = None
    rotsearch_d = None

    lamc = 4.3e-6
    oversample = 11
    bandpass = np.array([(1.0, lamc),])
    pixelscale_as = 0.0656
    arcsec2rad = u.arcsec.to(u.rad)
    PIXELSCALE_r = pixelscale_as * arcsec2rad
    holeshape = 'hex'
    filt = "F430M"
    rotsearch_d = np.arange(-3, 3.1, 1)

    affine2d = FAP.find_rotation(data[:,:], psf_offset_find_rotation, rotsearch_d,
                   mx, my, sx, sy, xo, yo,
                   PIXELSCALE_r, dim, bandpass, oversample, holeshape)

    niriss = InstrumentData.NIRISS(filt, bandpass=bandpass, affine2d=affine2d)

    ff_t = nrm_core.FringeFitter(niriss, psf_offset_ff=psf_offset_ff,
                oversample=oversample)

    output_model = ff_t.fit_fringes_all(input_model)

    # Copy header keywords from input to output
    output_model.update(input_model)

    return output_model
