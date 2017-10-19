from __future__ import (division, print_function, unicode_literals,
    absolute_import)

import numpy as np

from astropy import wcs as fitswcs
from astropy.coordinates import SkyCoord
from astropy.modeling.models import Scale, AffineTransformation2D
from astropy.modeling import Model
from gwcs import WCS, wcstools

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def compute_output_transform(refwcs, filename, fiducial):
    """Compute a simple FITS-type WCS transform
    """
    x0, y0 = refwcs.backward_transform(*fiducial)
    x1 = x0 + 1
    y1 = y0 + 1
    ra0, dec0 = refwcs(x0, y0)
    ra_xdir, dec_xdir = refwcs(x1, y0)
    ra_ydir, dec_ydir = refwcs(x0, y1)

    position0 = SkyCoord(ra=ra0, dec=dec0, unit='deg')
    position_xdir = SkyCoord(ra=ra_xdir, dec=dec_xdir, unit='deg')
    position_ydir = SkyCoord(ra=ra_ydir, dec=dec_ydir, unit='deg')
    offset_xdir = position0.spherical_offsets_to(position_xdir)
    offset_ydir = position0.spherical_offsets_to(position_ydir)

    xscale = np.abs(position0.separation(position_xdir).value)
    yscale = np.abs(position0.separation(position_ydir).value)
    scale = np.sqrt(xscale * yscale)

    c00 = offset_xdir[0].value / scale
    c01 = offset_xdir[1].value / scale
    c10 = offset_ydir[0].value / scale
    c11 = offset_ydir[1].value / scale
    pc_matrix = AffineTransformation2D(matrix=[[c00, c01], [c10, c11]])
    cdelt = Scale(scale) & Scale(scale)

    return pc_matrix | cdelt


def bounding_box_from_shape(shape):
    """ Create bounding_box for WCS based on shape of model data.
    """
    bb = []
    for s in reversed(shape):
        bb.append((-0.5, s - 0.5))
    return tuple(bb)


def shape_from_bounding_box(bounding_box):
    """ Return the size of the frame based on the provided domain
    """
    size = []
    for axs in bounding_box:
        delta = axs[1] - axs[0]
        size.append(int(delta + 0.5))
    return tuple(reversed(size))


def calc_gwcs_pixmap(in_wcs, out_wcs, shape=None):
    """ Return a pixel grid map from input frame to output frame.
    """
    if shape:
        bb = bounding_box_from_shape(shape)
        log.debug("Bounding box from data shape: {}".format(bb))
    else:
        bb = in_wcs.bounding_box
        log.debug("Bounding box from WCS: {}".format(in_wcs.bounding_box))

    grid = wcstools.grid_from_bounding_box(bb, step=(1, 1))
    pixmap = np.dstack(reproject(in_wcs, out_wcs)(grid[0], grid[1]))
    return pixmap


def reproject(wcs1, wcs2):
    """
    Given two WCSs or transforms return a function which takes pixel
    coordinates in the first WCS or transform and computes them in the second
    one. It performs the forward transformation of ``wcs1`` followed by the
    inverse of ``wcs2``.

    Parameters
    ----------
    wcs1, wcs2 : `~astropy.wcs.WCS` or `~gwcs.wcs.WCS` or `~astropy.modeling.Model`
        WCS objects.

    Returns
    -------
    _reproject : func
        Function to compute the transformations.  It takes x, y
        positions in ``wcs1`` and returns x, y positions in ``wcs2``.
    """

    if isinstance(wcs1, fitswcs.WCS):
        forward_transform = wcs1.all_pix2world
    elif isinstance(wcs1, WCS):
        forward_transform = wcs1.forward_transform
    elif issubclass(wcs1, Model):
        forward_transform = wcs1
    else:
        raise TypeError("Expected input to be astropy.wcs.WCS or gwcs.WCS "
            "object or astropy.modeling.Model subclass")

    if isinstance(wcs2, fitswcs.WCS):
        backward_transform = wcs2.all_world2pix
    elif isinstance(wcs2, WCS):
        backward_transform = wcs2.backward_transform
    elif issubclass(wcs2, Model):
        backward_transform = wcs2.inverse
    else:
        raise TypeError("Expected input to be astropy.wcs.WCS or gwcs.WCS "
            "object or astropy.modeling.Model subclass")

    def _reproject(x, y):
        sky = forward_transform(x, y)
        flat_sky = []
        for axis in sky:
            flat_sky.append(axis.flatten())
        det = backward_transform(*tuple(flat_sky))
        det_reshaped = []
        for axis in det:
            det_reshaped.append(axis.reshape(x.shape))
        return tuple(det_reshaped)
    return _reproject
