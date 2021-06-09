#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# TODO I don't much like the way the rotation point is handled.
# TODO Theoretically it could be removed entirely, but the way coords and image
# TODO handle it are different by default (lower-left vs center).

from astropy.io import fits
import numpy as np
from scipy.ndimage import shift, rotate
from scipy.optimize import minimize
import warnings

from .soss_centroids import get_soss_centroids


def transform_coords(angle, xshift, yshift, xpix, ypix, cenx=1024, ceny=50):
    """Apply a rotation and shift to the trace centroids positions. This
    assumes that the trace centroids are already in the CV3 coordinate system.

    :param angle: The angle by which to rotate the coordinates, in degrees.
    :param xshift: The shift to apply to the x-coordinates after rotating.
    :param yshift: The shift to apply to the y-coordinates after rotating.
    :param xpix: The x-coordinates to be transformed.
    :param ypix: The y-coordinates to be transformed.
    :param cenx: The x-coordinate around which to rotate.
    :param ceny: The y-coordinate around which to rotate.

    :type angle: float
    :type xshift: float
    :type yshift: float
    :type xpix: array[float]
    :type ypix: array[float]
    :type cenx: The x-coordinate around which to apply the rotation.
    :type ceny: The y-coordinate around which to apply the rotation.

    :returns: xrot, yrot - The rotated and shifted coordinates.
    :rtype: Tuple(array[float], array[float])
    """

    # Convert to numpy arrays.
    xpix = np.atleast_1d(xpix)
    ypix = np.atleast_1d(ypix)

    # Required rotation in the detector frame to match the data.
    angle = np.deg2rad(angle)
    rot_mat = np.array([[np.cos(angle), -np.sin(angle)],
                        [np.sin(angle), np.cos(angle)]])

    # Rotation center set to o1 trace centroid halfway along spectral axis.
    points = np.array([xpix - cenx, ypix - ceny])
    rot_points = rot_mat @ points
    rot_points[0] += cenx
    rot_points[1] += ceny

    # Apply the offsets.
    xrot = rot_points[0] + xshift
    yrot = rot_points[1] + yshift

    return xrot, yrot


def _chi_squared(transform, xref, yref, xdat, ydat):
    """"Compute the chi-squared statistic for fitting the reference positions
    to the true positions.

    :param transform: The transformation parameters.
    :param xref: The reference x-positions.
    :param yref: The reference y-positions.
    :param xdat: The data x-positions.
    :param ydat: The data y-positions.

    :type transform: Tuple, List, Array
    :type xref: array[float]
    :type xref: array[float]
    :type xdat: array[float]
    :type ydat: array[float]

    :returns: chisq - The chi-squared value of the model fit.
    :rtype: float
    """

    angle, xshift, yshift = transform

    # Calculate rotated reference positions.
    xrot, yrot = transform_coords(angle, xshift, yshift, xref, yref)

    # After rotation, need to re-sort the x-positions.
    sort = np.argsort(xrot)
    xrot, yrot = xrot[sort], yrot[sort]

    # Interpolate rotated model onto same x scale as data.
    ymod = np.interp(xdat, xrot, yrot)

    # Compute the chi-square.
    chisq = np.nansum((ydat - ymod)**2)

    return chisq


def solve_transform(scidata, scimask, xref, yref, subarray, verbose=False):
    """Given a science image, determine the centroids and find the simple
    transformation needed to match xcen_ref and ycen_ref to the image.

    :param scidata: the image of the SOSS trace.
    :param scimask: a boolean mask of pixls to be excluded.
    :param xref: a priori expectation of the trace x-positions.
    :param yref: a priori expectation of the trace y-positions.
    :param subarray: the subarray of the observations.
    :param verbose: If set True provide diagnostic information.

    :type scidata: array[float]
    :type scimask: array[bool]
    :type xref: array[float]
    :type yref: array[float]
    :type subarray: str
    :type verbose: bool

    :returns: simple_transform - Array containing the angle, x-shift and y-shift
        needed to match xcen_ref and ycen_ref to the image.
    :rtype: array[float]
    """

    # Remove any NaNs used to pad the xref, yref coordinates.
    mask = np.isfinite(xref) & np.isfinite(yref)
    xref = xref[mask]
    yref = yref[mask]

    # Get centroids from data.
    centroids = get_soss_centroids(scidata, mask=scimask, subarray=subarray,
                                   verbose=verbose)
    xdat = centroids['order 1']['X centroid']
    ydat = centroids['order 1']['Y centroid']

    # Set up the optimization problem.
    guess_transform = np.array([0.15, 1, 1])
    min_args = (xref, yref, xdat, ydat)

    # Find the best-fit transformation.
    result = minimize(_chi_squared, guess_transform, args=min_args)
    simple_transform = result.x

    return simple_transform


def rotate_image(image, angle, origin):
    """Rotate an image around a specific pixel.

    :param image: The image to rotate.
    :param angle: The rotation angle in degrees.
    :param origin: The x, y pixel position around which to rotate.

    :type image: array[float]
    :type angle: float
    :type origin: Tuple, List, Array

    :returns: image_rot - The rotated image.
    :rtype: array[float]
    """

    # Pad image so we can safely rotate around the origin.
    padx = [image.shape[1] - origin[0], origin[0]]
    pady = [image.shape[0] - origin[1], origin[1]]
    image_pad = np.pad(image, [pady, padx], 'constant')

    # Rotate the image.
    image_pad_rot = rotate(image_pad, angle, reshape=False)

    # Remove the padding.
    image_rot = image_pad_rot[pady[0]:-pady[1], padx[0]:-padx[1]]

    return image_rot


def transform_image(angle, xshift, yshift, image, cenx=1024, ceny=50):
    """Apply the rotation and offset to a 2D reference map, and bin
    the map down the native size and resolution.

    :param angle: The angle by which to rotate the file, in degrees.
    :param xshift: The x-shift to apply in native pixels, will be rounded to the
        nearest (oversampled) pixel.
    :param yshift: The y-shift to apply in native pixels, will be rounded to the
        nearest (oversampled) pixel.
    :param image: An image to transform.
    :param cenx: The x-coordinate around which to rotate.
    :param ceny: The y-coordinate around which to rotate.

    :type angle: float
    :type xshift: float
    :type yshift: float
    :type image: array[float]
    :type cenx: The x-coordinate around which to apply the rotation.
    :type ceny: The y-coordinate around which to apply the rotation.

    :returns: image_rot - The image, after applying the shift and rotation.
    :rtype: array[float]
    """

    # Rotate the image.
    image_rot = rotate_image(image, angle, [cenx, ceny])

    # Shift the image.
    image_rot = shift(image_rot, [yshift, xshift])

    return image_rot


def apply_transform(simple_transform, ref_map, oversample, pad, native=True,
                    norm=False):
    """Apply the transformation found by solve_transform() to a 2D reference map.

    :param simple_transform: The transformation parameters returned by
        solve_transform().
    :param ref_map: A reference map, e.g. a 2D Wavelength map or Trace Profile
        map.
    :param oversample: The oversampling factor the reference map.
    :param pad: The padding (in native pixels) on the reference map.
    :param native: If True bin down to native pixel sizes and remove padding.
        Default is True
    :param norm: If True normalize columns to 1, used for trace profile
        reference maps. Default is False.

    :type simple_transform: Tuple, List, Array
    :type ref_map: array[float]
    :type oversample: int
    :type pad: int
    :type native: bool
    :type norm: bool

    :returns: trans_maps - the ref_maps after having the transformation applied.
    :rtype: array[float]
    """

    ovs = oversample

    # Unpack the transformation.
    angle, xshift, yshift = simple_transform

    # Modify the transformation with the oversampling and padding.
    xshift = ovs*xshift
    yshift = ovs*yshift
    cenx = ovs*(pad + 1024)
    ceny = ovs*(pad + 50)

    # Set NaN pixels to zero - the rotation doesn't handle NaNs well.
    ref_map[np.isnan(ref_map)] = 0

    # Apply the transformation to the reference map.
    trans_map = transform_image(-angle, xshift, yshift, ref_map, cenx, ceny)

    if native:
        # Bin the transformed map down to native resolution.
        nrows, ncols = trans_map.shape
        trans_map = trans_map.reshape(nrows//ovs, ovs, ncols//ovs, ovs)
        trans_map = trans_map.mean(1).mean(-1)

        # Remove the padding.
        trans_map = trans_map[pad:-pad, pad:-pad]

    if norm:
        # Normalize so that the columns sum to 1.
        with warnings.catch_warnings():
            warnings.simplefilter(action="ignore", category=RuntimeWarning)
            trans_map = trans_map/np.nansum(trans_map, axis=0)

    return trans_map


def write_to_file(stack, filename):  # TODO function not needed?
    """Utility function to write transformed 2D trace profile or wavelength map
    files to disk. Data will be saved as a multi-extension fits file.

    Parameters
    ----------
    stack : np.ndarray (2xYx2048)
        Array containing transformed 2D trace profile or wavelength map data.
        The first dimension must be the spectral order, the second dimension
        the spatial dimension, and the third the spectral dimension.
    filename : str
        Name of the file to which to write the data.
    """

    hdu_p = fits.PrimaryHDU()
    hdulist = [hdu_p]
    for order in [1, 2]:
        hdu_o = fits.ImageHDU(data=stack[order-1])
        hdu_o.header['ORDER'] = order
        hdu_o.header.comments['ORDER'] = 'Spectral order.'
        hdulist.append(hdu_o)

    hdu = fits.HDUList(hdulist)
    hdu.writeto('{}.fits'.format(filename), overwrite=True)
    hdu.close()

    return


def main():
    """Placeholder for potential multiprocessing."""

    return


if __name__ == '__main__':
    main()
