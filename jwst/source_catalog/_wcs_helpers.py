# Licensed under a 3-clause BSD style license - see LICENSE.rst
# (taken from photutils: should probably migrate into astropy.wcs)
"""
This module provides WCS helper tools.
"""

import astropy.units as u
import numpy as np


def pixel_scale_angle_at_skycoord(skycoord, wcs, offset=1 * u.arcsec):
    """
    Calculate the pixel coordinate and scale and WCS rotation angle at
    the position of a SkyCoord coordinate.

    Parameters
    ----------
    skycoord : `~astropy.coordinates.SkyCoord`
        The SkyCoord coordinate.

    wcs : WCS object
        A world coordinate system (WCS) transformation that
        supports the `astropy shared interface for WCS
        <https://docs.astropy.org/en/stable/wcs/wcsapi.html>`_ (e.g.,
        `astropy.wcs.WCS`, `gwcs.wcs.WCS`).

    offset : `~astropy.units.Quantity`
        A small angular offset to use to compute the pixel scale and
        position angle.

    Returns
    -------
    xypos : tuple of float
        The (x, y) pixel coordinate.

    scale : `~astropy.units.Quantity`
        The pixel scale in arcsec/pixel.

    angle : `~astropy.units.Quantity`
        The angle (in degrees) measured counterclockwise from the
        positive x axis to the "North" axis of the celestial coordinate
        system.

    Notes
    -----
    If distortions are present in the image, the x and y pixel scales
    likely differ.  This function computes a single pixel scale along
    the North/South axis.
    """
    # Convert to pixel coordinates
    xpos, ypos = wcs.world_to_pixel(skycoord)

    # We take a point directly North (i.e., latitude offset) the
    # input sky coordinate and convert it to pixel coordinates,
    # then we use the pixel deltas between the input and offset sky
    # coordinate to calculate the pixel scale and angle.
    skycoord_offset = skycoord.directional_offset_by(0.0, offset)
    x_offset, y_offset = wcs.world_to_pixel(skycoord_offset)

    dx = x_offset - xpos
    dy = y_offset - ypos
    scale = offset.to(u.arcsec) / (np.hypot(dx, dy) * u.pixel)
    angle = (np.arctan2(dy, dx) * u.radian).to(u.deg)

    return (xpos, ypos), scale, angle
