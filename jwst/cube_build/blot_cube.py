""" Map the detector pixels to the cube coordinate system.
This is where the weight functions are used.
"""
import numpy as np
import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def blot_overlap(ipt, xstart,
                 xcenter, ycenter,
                 x_cube, y_cube,
                 flux_cube,
                 blot_xsize,
                 blot_flux, blot_weight):

    """ Blot the median sky image back to the detector

    ipt is the median element.
    xcenter, ycenter are the detector pixel arrays
    xstart is only valid for MIRI and is the start of the x detector value for channel
          0  for channels on left and ~512 for channels on right
    x_cube, y_cube:  median cube IFU mapped backwards to detector
    flux_cube: median flux
    blot_flux & blot_weight: blotted values of flux_cube to detector
    """

    roi_det = 1.0  # Just large enough that we don't get holes
    xdistance = np.absolute(x_cube[ipt] - xcenter)
    ydistance = np.absolute(y_cube[ipt] - ycenter)

    index_x = np.where(xdistance <= roi_det)
    index_y = np.where(ydistance <= roi_det)

    if len(index_x[0]) > 0 and len(index_y[0]) > 0:

        d1pix = x_cube[ipt] - xcenter[index_x]
        d2pix = y_cube[ipt] - ycenter[index_y]

        dxy = [(dx * dx + dy * dy) for dy in d2pix for dx in d1pix]
        dxy = np.array(dxy)
        dxy = np.sqrt(dxy)
        weight_distance = np.exp(-dxy)
        weighted_flux = weight_distance * flux_cube[ipt]

        index2d = [iy * blot_xsize + ix for iy in index_y[0] for ix in (index_x[0] + xstart)]
        index2d = np.array(index2d)

        blot_flux[index2d] = blot_flux[index2d] + weighted_flux
        blot_weight[index2d] = blot_weight[index2d] + weight_distance
