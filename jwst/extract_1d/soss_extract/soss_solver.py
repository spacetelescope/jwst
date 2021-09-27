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

from .soss_syscor import aperture_mask
from .soss_centroids import get_centroids_com


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


def evaluate_model(xmod, transform, xref, yref):
    """Evaluate the transformed reference coordinates at particular x-values.

    :param xmod: The x-values at which to evaluate the transformed coordinates.
    :param transform: The transformation parameters.
    :param xref: The reference x-positions.
    :param yref: The reference y-positions.

    :type xmod: array[float]
    :type transform: Tuple, List, Array
    :type xref: array[float]
    :type yref: array[float]

    :returns: ymod - The transformed y-coordinates corresponding to xmod.
    :rtype: array[float]
    """

    angle, xshift, yshift = transform

    # Calculate rotated reference positions.
    xrot, yrot = transform_coords(angle, xshift, yshift, xref, yref)

    # After rotation, need to re-sort the x-positions for interpolation.
    sort = np.argsort(xrot)
    xrot, yrot = xrot[sort], yrot[sort]

    # Interpolate rotated model onto same x scale as data.
    ymod = np.interp(xmod, xrot, yrot)

    return ymod


def _chi_squared(transform, xref_o1, yref_o1, xref_o2, yref_o2,
                 xdat_o1, ydat_o1, xdat_o2, ydat_o2):
    """Compute the chi-squared statistic for fitting the reference positions
    to the true positions.

    :param transform: The transformation parameters.
    :param xref_o1: The order 1 reference x-positions.
    :param yref_o1: The order 1 reference y-positions.
    :param xref_o2: The order 2 reference x-positions.
    :param yref_o2: The order 2 reference y-positions.
    :param xdat_o1: The order 1 data x-positions.
    :param ydat_o1: The order 1 data y-positions.
    :param xdat_o2: The order 2 data x-positions.
    :param ydat_o2: The order 2 data y-positions.

    :type transform: Tuple, List, Array
    :type xref_o1: array[float]
    :type xref_o1: array[float]
    :type xref_o2: array[float]
    :type xref_o2: array[float]
    :type xdat_o1: array[float]
    :type ydat_o1: array[float]
    :type xdat_o2: array[float]
    :type ydat_o2: array[float]

    :returns: chisq - The chi-squared value of the model fit.
    :rtype: float
    """

    # Interpolate rotated model onto same x scale as data.
    ymod_o1 = evaluate_model(xdat_o1, transform, xref_o1, yref_o1)
    ymod_o2 = evaluate_model(xdat_o2, transform, xref_o2, yref_o2)

    # Compute the chi-square.
    chisq_o1 = np.nansum((ydat_o1 - ymod_o1)**2)
    chisq_o2 = np.nansum((ydat_o2 - ymod_o2)**2)
    chisq = chisq_o1 + chisq_o2

    return chisq


def solve_transform(scidata, scimask, xref_o1, yref_o1, xref_o2, yref_o2,
                    halfwidth=30.):
    """Given a science image, determine the centroids and find the simple
    transformation needed to match xcen_ref and ycen_ref to the image.

    :param scidata: the image of the SOSS trace.
    :param scimask: a boolean mask of pixls to be excluded.
    :param xref_o1: a priori expectation of the order 1 trace x-positions.
    :param yref_o1: a priori expectation of the order 1 trace y-positions.
    :param xref_o2: a priori expectation of the order 2 trace x-positions.
    :param yref_o2: a priori expectation of the order 2 trace y-positions.
    :param halfwidth: size of the aperture mask used when extracting the trace
        positions from the data.

    :type scidata: array[float]
    :type scimask: array[bool]
    :type xref_o1: array[float]
    :type yref_o1: array[float]
    :type xref_o2: array[float]
    :type yref_o2: array[float]
    :type halfwidth: float

    :returns: simple_transform - Array containing the angle, x-shift and y-shift
        needed to match xcen_ref and ycen_ref to the image.
    :rtype: array[float]
    """

    # Remove any NaNs used to pad the xref, yref coordinates.
    mask = np.isfinite(xref_o1) & np.isfinite(yref_o1)
    xref_o1 = xref_o1[mask]
    yref_o1 = yref_o1[mask]

    mask = np.isfinite(xref_o2) & np.isfinite(yref_o2)
    xref_o2 = xref_o2[mask]
    yref_o2 = yref_o2[mask]

    # Get centroids from data.
    aper_mask_o1 = aperture_mask(xref_o1, yref_o1, halfwidth, scidata.shape)
    mask = aper_mask_o1 | scimask
    xdat_o1, ydat_o1, _ = get_centroids_com(scidata, mask=mask, poly_order=None)

    aper_mask_o2 = aperture_mask(xref_o2, yref_o2, halfwidth, scidata.shape)
    mask = aper_mask_o2 | scimask
    xdat_o2, ydat_o2, _ = get_centroids_com(scidata, mask=mask, poly_order=None)

    # Use only the clean range between x=800 and x=1700.
    mask = (xdat_o1 >= 800) & (xdat_o1 <= 1700)
    xdat_o1 = xdat_o1[mask]
    ydat_o1 = ydat_o1[mask]

    mask = (xdat_o2 >= 800) & (xdat_o2 <= 1700)
    xdat_o2 = xdat_o2[mask]
    ydat_o2 = ydat_o2[mask]

    # Set up the optimization problem.
    guess_transform = np.array([0., 0., 0.])
    min_args = (xref_o1, yref_o1, xref_o2, yref_o2,
                xdat_o1, ydat_o1, xdat_o2, ydat_o2)

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


def apply_transform(simple_transform, ref_map, oversample, pad, native=True):
    """Apply the transformation found by solve_transform() to a 2D reference map.

    :param simple_transform: The transformation parameters returned by
        solve_transform().
    :param ref_map: A reference map, e.g. a 2D Wavelength map or Trace Profile
        map.
    :param oversample: The oversampling factor the reference map.
    :param pad: The padding (in native pixels) on the reference map.
    :param native: If True bin down to native pixel sizes and remove padding.
        Default is True

    :type simple_transform: Tuple, List, Array
    :type ref_map: array[float]
    :type oversample: int
    :type pad: int
    :type native: bool

    :returns: trans_map - the ref_map after having the transformation applied.
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

    # Apply the transformation to the reference map.
    trans_map = transform_image(-angle, xshift, yshift, ref_map, cenx, ceny)

    if native:

        # Bin the transformed map down to native resolution.
        nrows, ncols = trans_map.shape
        trans_map = trans_map.reshape(nrows//ovs, ovs, ncols//ovs, ovs)
        trans_map = trans_map.mean(1).mean(-1)

        # Remove the padding.
        trans_map = trans_map[pad:-pad, pad:-pad]

    return trans_map


def transform_wavemap(simple_transform, wavemap, oversample, pad, native=True):
    """Apply the transformation found by solve_transform() to a 2D reference map.

    :param simple_transform: The transformation parameters returned by
        solve_transform().
    :param wavemap: A reference wavelength map.
    :param oversample: The oversampling factor the reference map.
    :param pad: The padding (in native pixels) on the reference map.
    :param native: If True bin down to native pixel sizes and remove padding.
        Default is True

    :type simple_transform: Tuple, List, Array
    :type wavemap: array[float]
    :type oversample: int
    :type pad: int
    :type native: bool

    :returns: trans_wavemap - the wavemap after having the transformation applied.
    :rtype: array[float]
    """

    # Find the minimum and maximum wavelength of the wavemap.
    minval = np.nanmin(wavemap)
    maxval = np.nanmax(wavemap)

    # Set NaNs to zero to prevent errors when shifting/rotating.
    mask = np.isnan(wavemap)
    wavemap[mask] = 0.

    # Apply the transformation to the wavemap.
    trans_wavemap = apply_transform(simple_transform, wavemap, oversample, pad, native=native)

    # Set pixels with interpolation artefacts zero by enforcing the original min/max.
    mask = (trans_wavemap < minval) | (trans_wavemap > maxval)
    trans_wavemap[mask] = 0

    return trans_wavemap


def transform_profile(simple_transform, profile, oversample, pad, native=True, norm=True):
    """Apply the transformation found by solve_transform() to a 2D reference map.

    :param simple_transform: The transformation parameters returned by
        solve_transform().
    :param profile: A reference trace profile map.
    :param oversample: The oversampling factor the reference map.
    :param pad: The padding (in native pixels) on the reference map.
    :param native: If True bin down to native pixel sizes and remove padding.
        Default is True
    :param norm: If True (re-)normalize columns so they sum to 1, as expected by the engine.

    :type simple_transform: Tuple, List, Array
    :type profile: array[float]
    :type oversample: int
    :type pad: int
    :type native: bool

    :returns: trans_profile - the trace profile after having the transformation applied.
    :rtype: array[float]
    """

    # Apply the transformation to the wavemap.
    trans_profile = apply_transform(simple_transform, profile, oversample, pad, native=native)

    if norm:

        # Normalize so that the columns sum to 1.
        with warnings.catch_warnings():
            warnings.simplefilter(action="ignore", category=RuntimeWarning)
            trans_profile = trans_profile/np.nansum(trans_profile, axis=0)

        trans_profile[~np.isfinite(trans_profile)] = 0.

    return trans_profile


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
