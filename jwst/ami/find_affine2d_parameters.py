#
#  Module for calculation of the optimal rotation for an image plane
#

import logging
import numpy as np

from ..fringefitting.LG_Model import NRM_Model
from ..misctools import utils

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

def create_afflist_rot(rotdegs):
    """
    Short Summary
    -------------
    Create a list of affine objects with various rotations to use in order to
    go through and find which fits an image plane data best.

    Parameters
    ----------
    rotdegs: float 1D array
        Search window for rotation fine_tuning, in degrees

    Returns
    -------
    alist: list
        Affine2d objects having various rotations
    """
    alist = []
    for nrot, rotd in enumerate(rotdegs):
        rotd_ = utils.avoidhexsingularity(rotd)
        alist.append(utils.Affine2d(rotradccw=np.pi*rotd_/180.0,
                                    name="affrot_{0:+.3f}".format(rotd_)))
    return alist


def find_rotation(imagedata, psf_offset, rotdegs, mx, my, sx, sy, xo, yo,
                  pixel, npix, bandpass, over, holeshape):
    """
    Short Summary
    -------------

    Parameters
    ----------
    imagedata:
    psf_offset:
    rotdegs:
    mx, my:
    sx,sy:
    xo,yo:
    pixel:
    npix:
    bandpass:
    over:
    holeshape:

    Returns
    -------

    """

    if hasattr(rotdegs, '__iter__') is False:
        rotdegs = (rotdegs,)

    affine2d_list = create_afflist_rot(rotdegs, mx, my, sx,sy, xo,yo)
    crosscorr_rots = []

    for (rot,aff) in zip(rotdegs,affine2d_list):

        jw = NRM_Model(mask='jwst', holeshape=holeshape,
                       over=over, affine2d=aff)
        jw.set_pixelscale(pixel)
        # psf_offset in data coords & pixels.  Does it get rotated?  Second order errors poss.
        #  Some numerical testing needed for big eg 90 degree affine2d rotations.  Later.
        jw.simulate(fov=npix, bandpass=bandpass, over=over, psf_offset=psf_offset)

        crosscorr_rots.append(utils.rcrosscorrelate(imagedata, jw.psf).max())
        del jw

    rot_measured_d, max_cor = utils.findpeak_1d(crosscorr_rots, rotdegs)

    # return convenient affine2d
    return utils.Affine2d(rotradccw=np.pi*rot_measured_d/180.0,
                          name="{0:.4f}".format(rot_measured_d))
