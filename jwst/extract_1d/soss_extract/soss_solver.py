import logging

import numpy as np
from scipy.ndimage import shift, rotate
from scipy.optimize import minimize
import warnings

from .soss_syscor import aperture_mask
from .soss_centroids import get_centroids_com

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def transform_coords(angle, xshift, yshift, xpix, ypix, cenx=1024, ceny=50):
    """Apply a rotation and shift to the trace centroids positions. This
    assumes that the trace centroids are already in the CV3 coordinate system.

    Parameters
    ----------
    angle : float
        The angle by which to rotate the coordinates, in degrees.
    xshift : float
        The shift to apply to the x-coordinates after rotating.
    yshift : float
        The shift to apply to the y-coordinates after rotating.
    xpix : array[float]
        The x-coordinates to be transformed.
    ypix : array[float]
        The y-coordinates to be transformed.
    cenx : float (optional)
        The x-coordinate around which to rotate.
    ceny : float (optional)
        The y-coordinate around which to rotate.

    Returns
    -------
    xrot : array[float]
        The rotated and shifted x coordinate.
    yrot : array[float]
        The rotated and shifted y coordinate.
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

    Parameters
    ----------
    xmod : array[float]
        The x-values at which to evaluate the transformed coordinates.
    transform : array[float]
        The transformation parameters.
    xref : array[float]
        The reference x-positions.
    yref : array[float]
        The reference y-positions.

    Returns
    -------
    ymod : array[float]
        The transformed y-coordinates corresponding to xmod.
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

    Parameters
    ----------
    transform : array[float]
        The transformation parameters.
    xref_o1 : array[float]
        The order 1 reference x-positions.
    yref_o1 : array[float]
        The order 1 reference y-positions.
    xref_o2 : array[float]
        The order 2 reference x-positions.
    yref_o2 : array[float]
        The order 2 reference y-positions.
    xdat_o1 : array[float]
        The order 1 data x-positions.
    ydat_o1 : array[float]
        The order 1 data y-positions.
    xdat_o2 : array[float]
        The order 2 data x-positions.
    ydat_o2 : array[float]
        The order 2 data y-positions.

    Returns
    -------
    chisq : float
        The chi-squared value of the model fit.
    """
    # Interpolate rotated model of first order onto same x scale as data.
    ymod_o1 = evaluate_model(xdat_o1, transform, xref_o1, yref_o1)

    # Compute the chi-square.
    chisq_o1 = np.nansum((ydat_o1 - ymod_o1)**2)

    # If second order centroids are provided, include them in the calculation.
    if xdat_o2 is not None:
        # Interpolate rotated model onto same x scale as data.
        ymod_o2 = evaluate_model(xdat_o2, transform, xref_o2, yref_o2)

        # Compute the chi-square and add to the first order.
        chisq_o2 = np.nansum((ydat_o2 - ymod_o2)**2)
        chisq = chisq_o1 + chisq_o2
    # If not, use only the first order.
    else:
        chisq = chisq_o1

    return chisq


def solve_transform(scidata_bkg, scimask, xref_o1, yref_o1, xref_o2=None,
                    yref_o2=None, halfwidth=30., is_fitted=(True, True, True),
                    soss_filter='CLEAR', guess_transform=(None, None, None),
                    bounds_theta=(-5., 5.), bounds_x=(-3, 3), bounds_y=(-3., 3.)):
    """Given a science image, determine the centroids and find the simple
    transformation (rotation + vertical & horizonal offset, or some combination
    thereof) needed to match xref_o1 and yref_o1 to the image.

    Parameters
    ----------
    scidata_bkg : array[float]
        A background subtracted image of the SOSS trace.
    scimask : array[float]
        A boolean mask of pixels to be excluded.
    xref_o1 : array[float]
        A priori expectation of the order 1 trace x-positions.
    yref_o1 : array[float]
        A priori expectation of the order 1 trace y-positions.
    xref_o2 : array[float] (optional)
        A priori expectation of the order 2 trace x-positions. Providing these
        will improve the accuracy of the solver.
    yref_o2 : array[float] (optional)
        A priori expectation of the order 2 trace y-positions. Providing these
        will improve the accuracy of the solver.
    halfwidth : float (optional)
        Size of the aperture mask used when extracting the trace positions
        from the data.
    is_fitted: tuple, list or array [bool]
        Tuple indicating which parameter is fitted in transform.
        It is ordered as (rotation, x-shift, y-shift).
        Default is (True, True, True), so all parameters are fitted.
        When False, the value is fixed to the value in `guess_transform`.
    guess_transform: tuple, list or array
        Tuple of the initial guess value of each transformation parameters.
        It is ordered as (rotation, x-shift, y-shift). If set to None,
        the value 0. is given. Default is (None, None, None).
    soss_filter : str (optional)
        Designator for the SOSS filter used in the observation. Either CLEAR
        or F277W. Setting F277W here will force shift to False.
    bounds_theta : array[float] (optional)
        Boundaries on the rotation angle to consider in the Chi-squared
        minimization.
    bounds_x : array[float] (optional)
        Boundaries on the horizontal offset to consider in the Chi-squared
        minimization.
    bounds_y : array[float] (optional)
        Boundaries on the vertical offset to consider in the Chi-squared
        minimization.

    Returns
    -------
    simple_transform : array[float]
        Array containing the angle, x-shift and y-shift needed to match
        xref_o1 and yref_o1 to the image.
    """
    # Convert None to 0. in guess_transform and then convert to array
    guess_transform = [val if val is not None else 0. for val in guess_transform]
    guess_transform = np.array(guess_transform)

    # Convert is_fitted to array
    is_fitted = np.array(is_fitted)

    # Start with order 1 centroids as they will be available for all subarrays.
    # Remove any NaNs used to pad the xref, yref coordinates.
    mask_o1 = np.isfinite(xref_o1) & np.isfinite(yref_o1)
    xref_o1 = xref_o1[mask_o1]
    yref_o1 = yref_o1[mask_o1]

    # Get centroids from data.
    aper_mask_o1 = aperture_mask(xref_o1, yref_o1, halfwidth, scidata_bkg.shape)
    mask = aper_mask_o1 | scimask
    xdat_o1, ydat_o1, _ = get_centroids_com(scidata_bkg, mask=mask,
                                            poly_order=None)

    # If order 2 centroids are provided, include them in the analysis. The
    # inclusion of the order 2 centroids will allow for a more accurate
    # determination of the rotation and offset, as the addition of the second
    # order provides an anchor in the spatial direction. However, there are
    # instances (a SUBSTRIP96, or F277W observation for example) where the
    # second order is not available. In this case, work only with order 1.
    if xref_o2 is not None and yref_o2 is not None and (soss_filter == 'CLEAR' or soss_filter == 'FULL'):
        # Remove any NaNs used to pad the xref, yref coordinates.
        log.info('Measuring trace position for orders 1 and 2.')
        mask_o2 = np.isfinite(xref_o2) & np.isfinite(yref_o2)
        xref_o2 = xref_o2[mask_o2]
        yref_o2 = yref_o2[mask_o2]

        # Get centroids from data.
        aper_mask_o2 = aperture_mask(xref_o2, yref_o2, halfwidth, scidata_bkg.shape)
        mask = aper_mask_o2 | scimask
        xdat_o2, ydat_o2, _ = get_centroids_com(scidata_bkg, mask=mask,
                                                poly_order=None)

        # Use only the uncontaminated range between x=800 and x=1700.
        mask = (xdat_o1 >= 800) & (xdat_o1 <= 1700)
        xdat_o1 = xdat_o1[mask]
        ydat_o1 = ydat_o1[mask]

        mask = (xdat_o2 >= 800) & (xdat_o2 <= 1700)
        xdat_o2 = xdat_o2[mask]
        ydat_o2 = ydat_o2[mask]

    elif soss_filter == 'F277W':
        # If the exposure uses the F277W filter, there is no second order, and
        # first order centroids are only useful at wavelengths greater than 2.5µm.
        # Restrict centroids to lie within region lambda>~2.5µm, where the
        # F277W filter response is strong.
        log.info('Measuring trace position for order 1 spanning the F277W pixels.')
        mask = (xdat_o1 >= 25) & (xdat_o1 <= 425)
        xdat_o1 = xdat_o1[mask]
        ydat_o1 = ydat_o1[mask]
        # Force shift to False as there is not enough information to
        # constrain dx, dy and dtheta simultaneously.
        is_fitted[1:] = False
        # Second order centroids are not available.
        xdat_o2, ydat_o2 = None, None

    else:
        # If the exposure is SUBSTRIP96 using the CLEAR filter, there is no
        # order 2. Use the entire first order to enable the maximum possible
        # positional constraint on the centroids.
        log.info('Measuring trace position for order 1 only.')
        xdat_o2, ydat_o2 = None, None

    # Find the simple transformation via a chi-squared minimization of the
    # extracted and reference centroids. This transformation considers by
    # default rotation as well as vertical and horizontal offsets.

    # Setup minimizing function depending on the parameters to fit
    simple_transform = guess_transform.copy()

    def _chi2_to_fit(values, *args):
        simple_transform[is_fitted] = values
        return _chi_squared(simple_transform, *args)

    # Set up the optimization problem.
    min_args = (xref_o1, yref_o1, xref_o2, yref_o2,
                xdat_o1, ydat_o1, xdat_o2, ydat_o2)

    # Define the boundaries
    bounds = [bounds_theta, bounds_x, bounds_y]
    # Only keep the ones that are fitted
    bounds = [bnds for idx, bnds in enumerate(bounds) if is_fitted[idx]]

    # Find the best-fit transformation.
    result = minimize(_chi2_to_fit, guess_transform[is_fitted], bounds=bounds,
                      args=min_args)
    simple_transform[is_fitted] = result.x

    return simple_transform


def rotate_image(image, angle, origin):
    """Rotate an image around a specific pixel.

    Parameters
    ----------
    image : array[float]
        The image to rotate.
    angle : float
        The rotation angle in degrees.
    origin : Tuple, List, Array
        The x and y pixel position around which to rotate.

    Returns
    -------
    image_rot : array[float]
        The rotated image.
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
    """Apply the transformation found by solve_transform() to a 2D reference
     map.

    Parameters
    ----------
    angle : float
        The angle by which to rotate the file, in degrees.
    xshift : float
        The x-shift to apply in native pixels, will be rounded to the
        nearest (oversampled) pixel.
    yshift : float
        The y-shift to apply in native pixels, will be rounded to the
        nearest (oversampled) pixel.
    image : array[float]
        An image to transform.
    cenx : float (optional)
        The x-coordinate around which to rotate.
    ceny : float (optional)
        The y-coordinate around which to rotate.

    Returns
    -------
    image_rot : array[float]
        The image, after applying the shift and rotation.
    """

    # Rotate the image.
    image_rot = rotate_image(image, angle, [cenx, ceny])

    # Shift the image.
    image_rot = shift(image_rot, [yshift, xshift])

    return image_rot


def apply_transform(simple_transform, ref_map, oversample, pad, native=True):
    """Apply the calculated rotation and offset to a 2D reference map, and bin
     the map down the native size and resolution.

    Parameters
    ----------
    simple_transform : Tuple, List, Array
        The transformation parameters returned by solve_transform().
    ref_map : array[float]
        A reference map: e.g., a 2D Wavelength map or Trace Profile map.
    oversample : int
        The oversampling factor the reference map.
    pad : int
        The padding (in native pixels) on the reference map.
    native : bool (optional)
        If True bin down to native pixel sizes and remove padding.

    Returns
    -------
    trans_map : array[float]
        The ref_map after having the transformation applied.
    """

    ovs = oversample

    # Unpack the transformation.
    angle, xshift, yshift = simple_transform

    # Modify the transformation with the oversampling and padding.
    xshift = ovs * xshift
    yshift = ovs * yshift
    cenx = ovs * (pad + 1024)
    ceny = ovs * (pad + 50)

    # Apply the transformation to the reference map.
    if (angle == 0 and xshift == 0 and yshift == 0):
        trans_map = ref_map
    else:
        trans_map = transform_image(-angle, xshift, yshift, ref_map, cenx, ceny)

    if native:
        if ovs > 1:
            # Bin the transformed map down to native resolution.
            nrows, ncols = trans_map.shape
            trans_map = trans_map.reshape((nrows // ovs), ovs, (ncols // ovs), ovs)
            trans_map = trans_map.mean(1).mean(-1)

        # Remove the padding.
        if pad != 0:
            trans_map = trans_map[pad:-pad, pad:-pad]

    return trans_map


def transform_wavemap(simple_transform, wavemap, oversample, pad, native=True):
    """Apply the transformation found by solve_transform() to a 2D reference
     wavelength map.

    Parameters
    ----------
    simple_transform : Tuple, List, Array
        The transformation parameters returned by solve_transform().
    wavemap : array[float]
        A reference 2D wavelength map.
    oversample : int
        The oversampling factor the reference map.
    pad : int
        The padding (in native pixels) on the reference map.
    native : bool (optional)
        If True bin down to native pixel sizes and remove padding.

    Returns
    -------
    trans_wavemap : array[float]
        The ref_map after having the transformation applied.
    """

    # Find the minimum and maximum wavelength of the wavelength map.
    minval = np.nanmin(wavemap)
    maxval = np.nanmax(wavemap)

    # Set NaNs to zero to prevent errors when shifting/rotating.
    mask = np.isnan(wavemap)
    wavemap = np.where(mask, 0., wavemap)

    # Apply the transformation to the wavelength map.
    trans_wavemap = apply_transform(simple_transform, wavemap, oversample,
                                    pad, native=native)

    # Set pixels with interpolation artifacts to zero by enforcing the
    # original min/max.
    mask = (trans_wavemap < minval) | (trans_wavemap > maxval)
    trans_wavemap[mask] = 0

    return trans_wavemap


def transform_profile(simple_transform, profile, oversample, pad, native=True,
                      norm=True):
    """Apply the transformation found by solve_transform() to a 2D reference
     trace profile map.

    Parameters
    ----------
    simple_transform : Tuple, List, Array
        The transformation parameters returned by solve_transform().
    profile : array[float]
        A reference 2D trace profile map.
    oversample : int
        The oversampling factor the reference map.
    pad : int
        The padding (in native pixels) on the reference map.
    native : bool (optional)
        If True bin down to native pixel sizes and remove padding.
    norm : bool (optional)
        If True, normalize each column of the trace profile to sum to one.

    Returns
    -------
    trans_profile : array[float]
        The ref_map after having the transformation applied.
    """

    # Apply the transformation to the 2D trace map.
    trans_profile = apply_transform(simple_transform, profile, oversample,
                                    pad, native=native)

    if norm:
        # Normalize so that the columns sum to 1.
        with warnings.catch_warnings():
            warnings.simplefilter(action="ignore", category=RuntimeWarning)
            trans_profile = trans_profile / np.nansum(trans_profile, axis=0)

        trans_profile[~np.isfinite(trans_profile)] = 0.

    return trans_profile
