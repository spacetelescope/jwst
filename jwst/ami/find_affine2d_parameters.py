import logging
import numpy as np

from . import lg_model
from . import utils

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def create_afflist_rot(rotdegs):
    """
    Create a list of affine objects with various rotations.

    Find which affine rotation fits an image plane data best.

    Parameters
    ----------
    rotdegs : float 1D array
        Search window for rotation fine_tuning, in degrees

    Returns
    -------
    alist : list
        Affine2d objects having various rotations
    """
    alist = []
    for rotd in rotdegs:
        rotd_ = utils.avoidhexsingularity(rotd)
        alist.append(utils.Affine2d(rotradccw=np.pi * rotd_ / 180.0, name=f"affrot_{rotd_:+.3f}"))
    return alist


def find_rotation(
    imagedata, nrm_model, psf_offset, rotdegs, pixel, npix, bandpass, over, holeshape
):
    """
    Create an affine2d object using the known rotation and scale.

    Parameters
    ----------
    imagedata : 2D float array
        Image data
    nrm_model : NRMModel datamodel
        Datamodel containing mask geometry information
    psf_offset : 2D float array
        Offset from image center in detector pixels
    rotdegs : list of floats
        Range of rotations to search (degrees)
    pixel : float
        Pixel size in radians
    npix : int
        Number of detector pixels on a side
    bandpass : 2D float array, default=None
        Array of the form: [(weight1, wavl1), (weight2, wavl2), ...]
    over : int
        Oversampling factor
    holeshape : str
        Shape of hole; possible values are 'circ', 'hex', and 'fringe'

    Returns
    -------
    new_affine2d : Affine2d object
        Affine2d object using the known rotation and scale.
    """
    if hasattr(rotdegs, "__iter__") is False:
        rotdegs = (rotdegs,)

    affine2d_list = create_afflist_rot(rotdegs)

    crosscorr_rots = []

    for _rot, aff in zip(rotdegs, affine2d_list, strict=False):
        jw = lg_model.LgModel(
            nrm_model,
            bandpass=bandpass,
            mask="jwst_ami",
            holeshape=holeshape,
            over=over,
            affine2d=aff,
            pixscale=pixel,
        )

        # psf_offset in data coords & pixels.  Does it get rotated?  Second order errors poss.
        #  Some numerical testing needed for big eg 90 degree affine2d rotations.  Later.
        jw.simulate(fov=npix, psf_offset=psf_offset)

        crosscorr_rots.append(utils.rcrosscorrelate(imagedata, jw.psf).max())
        del jw

    rot_measured_d, _max_cor = utils.findpeak_1d(rotdegs, crosscorr_rots)

    # return convenient affine2d
    new_affine2d = utils.Affine2d(
        rotradccw=np.pi * rot_measured_d / 180.0, name=f"{rot_measured_d:.4f}"
    )

    return new_affine2d
