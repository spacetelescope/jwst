from functools import partial
import logging
import numpy as np
from scipy.interpolate import interp1d

log = logging.getLogger(__name__)

SOSS_XDIM = 2048
SOSS_YDIM = 300
XTRACE_ORD1_LEN = SOSS_XDIM
XTRACE_ORD2_LEN = 1783
WAVEMAP_WLMIN = 0.5
WAVEMAP_WLMAX = 5.5
WAVEMAP_NWL = 5001
SUBARRAY_YMIN = 2048 - 256

__all__ = ["get_soss_traces", "get_soss_wavemaps"]


def _get_wavelengths(refmodel, x, pwcpos, order):
    """
    Get the associated wavelength values for a given spectral order.

    Parameters
    ----------
    refmodel : PastasossModel
        The reference model holding the wavecal models and scale extents
    x : float or numpy.ndarray
        The input pixel values for which the function will estimate wavelengths
    pwcpos : float
        The position of the pupil wheel; used to determine
        the difference between current and commanded position to rotate the model
    order : int
        The spectral order to find trace and wavecal model indices for.

    Returns
    -------
    wavelengths : numpy.ndarray
        The estimated wavelengths for the given pixel values.
    """
    if order == 1:
        wavelengths = wavecal_model_order1_poly(refmodel, x, pwcpos)
    elif order == 2:
        wavelengths = wavecal_model_order2_poly(refmodel, x, pwcpos)

    return wavelengths


def _min_max_scaler(x, x_min, x_max):
    """
    Apply min-max scaling to input values.

    Parameters
    ----------
    x : float or numpy.ndarray
        The input value(s) to be scaled.
    x_min : float
        The minimum value in the range to which 'x' will be scaled.
    x_max : float
        The maximum value in the range to which 'x' will be scaled.

    Returns
    -------
    float or numpy.ndarray
        The scaled value(s) in the range [0, 1].

    Notes
    -----
    Min-max scaling is a data normalization technique that scales input values
    'x' to the range [0, 1] based on the provided minimum and maximum values,
    'x_min' and 'x_max'. This function is applicable to both individual values
    and arrays of values. This function will use the min/max values from the
    training data of the wavecal model.
    """
    x_scaled = (x - x_min) / (x_max - x_min)
    return x_scaled


def wavecal_model_order1_poly(refmodel, x, pwcpos):
    """
    Compute order 1 wavelengths.

    Parameters
    ----------
    refmodel : PastasossModel
        The reference model holding the wavecal models
        and scale extents
    x : float or numpy.ndarray
        The input pixel values for which the function
        will estimate wavelengths
    pwcpos : float
        The position of the pupil wheel; used to determine
        the difference between current and commanded position
        to rotate the model

    Returns
    -------
    wavelengths : numpy.ndarray
        The estimated wavelengths for the given pixel values.
    """
    x_scaler = partial(
        _min_max_scaler,
        **{
            "x_min": refmodel.wavecal_models[0].scale_extents[0][0],
            "x_max": refmodel.wavecal_models[0].scale_extents[1][0],
        },
    )

    pwcpos_offset_scaler = partial(
        _min_max_scaler,
        **{
            "x_min": refmodel.wavecal_models[0].scale_extents[0][1],
            "x_max": refmodel.wavecal_models[0].scale_extents[1][1],
        },
    )

    def get_poly_features(x, offset):
        """
        Polynomial features for the order 1 wavecal model.

        Parameters
        ----------
        x : float or numpy.ndarray
            The input pixel values for which the function will estimate wavelengths
        offset : float or numpy.ndarray
            The offset values for the pupil wheel position

        Returns
        -------
        numpy.ndarray
            The polynomial features for the order 1 wavecal model.
        """
        poly_features = np.array(
            [
                x,
                offset,
                x**2,
                x * offset,
                offset**2,
                x**3,
                x**2 * offset,
                x * offset**2,
                offset**3,
                x**4,
                x**3 * offset,
                x**2 * offset**2,
                x * offset**3,
                offset**4,
                x**5,
                x**4 * offset,
                x**3 * offset**2,
                x**2 * offset**3,
                x * offset**4,
                offset**5,
            ]
        )
        return poly_features

    # extract model weights and intercept
    coef = refmodel.wavecal_models[0].coefficients

    # get pixel columns and then scaled
    x_scaled = x_scaler(x)

    # offset
    offset = np.ones_like(x) * (pwcpos - refmodel.meta.pwcpos_cmd)
    offset_scaled = pwcpos_offset_scaler(offset)

    # polynomial features
    poly_features = get_poly_features(x_scaled, offset_scaled)
    wavelengths = coef[0] + coef[1:] @ poly_features

    return wavelengths


def wavecal_model_order2_poly(refmodel, x, pwcpos):
    """
    Compute order 2 wavelengths.

    Parameters
    ----------
    refmodel : PastasossModel
        The reference model holding the wavecal models
        and scale extents
    x : float or numpy.ndarray
        The input pixel values for which the function
        will estimate wavelengths
    pwcpos : float
        The position of the pupil wheel; used to determine
        the difference between current and commanded position
        to rotate the model

    Returns
    -------
    wavelengths : numpy.ndarray
        The estimated wavelengths for the given pixel values.
    """
    x_scaler = partial(
        _min_max_scaler,
        **{
            "x_min": refmodel.wavecal_models[1].scale_extents[0][0],
            "x_max": refmodel.wavecal_models[1].scale_extents[1][0],
        },
    )

    pwcpos_offset_scaler = partial(
        _min_max_scaler,
        **{
            "x_min": refmodel.wavecal_models[1].scale_extents[0][1],
            "x_max": refmodel.wavecal_models[1].scale_extents[1][1],
        },
    )

    def get_poly_features(x, offset):
        """
        Polynomial features for the order 2 wavecal model.

        Parameters
        ----------
        x : float or numpy.ndarray
            The input pixel values for which the function will estimate wavelengths
        offset : float or numpy.ndarray
            The offset values for the pupil wheel position

        Returns
        -------
        numpy.ndarray
            The polynomial features for the order 2 wavecal model.
        """
        poly_features = np.array(
            [
                x,
                offset,
                x**2,
                x * offset,
                offset**2,
                x**3,
                x**2 * offset,
                x * offset**2,
                offset**3,
            ]
        )
        return poly_features

    # coef and intercept
    coef = refmodel.wavecal_models[1].coefficients

    # get pixel columns and then scaled
    x_scaled = x_scaler(x)

    offset = np.ones_like(x) * pwcpos
    offset_scaled = pwcpos_offset_scaler(offset)

    # polynomial features
    poly_features = get_poly_features(x_scaled, offset_scaled)
    wavelengths = coef[0] + coef[1:] @ poly_features

    return wavelengths


def _rotate(x, y, angle, origin=(0, 0)):
    """
    Apply a rotation transformation to a set of 2D points.

    Parameters
    ----------
    x : np.ndarray
        The x-coordinates of the points to be transformed.
    y : np.ndarray
        The y-coordinates of the points to be transformed.
    angle : float
        The angle (in degrees) by which to rotate the points.
    origin : Tuple[float, float], optional
        The point about which to rotate the points. Default is (0, 0).

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        The x and y coordinates of the rotated points.

    Examples
    --------
    >>> x = np.array([0, 1, 2, 3])
    >>> y = np.array([0, 1, 2, 3])
    >>> x_rot, y_rot = _rotate(x, y, 90)
    """
    # shift to rotate about center
    xy_center = np.atleast_2d(origin).T
    xy = np.vstack([x, y])

    # Rotation transform matrix
    radians = np.radians(angle)
    c, s = np.cos(radians), np.sin(radians)
    rotation_matrix = np.array([[c, -s], [s, c]])

    # apply transformation
    x_new, y_new = rotation_matrix @ (xy - xy_center) + xy_center

    # interpolate rotated positions onto x-pixel column values
    # interpolate new coordinates onto original x values and mask values
    # outside of the domain of the image 0<=x<=2047 and 0<=y<=255.
    y_new = interp1d(x_new, y_new, fill_value="extrapolate")(x)
    mask = np.where(y_new <= 255.0)
    x = x[mask]
    y_new = y_new[mask]
    return x, y_new


def _find_spectral_order_index(refmodel, order):
    """
    Return index of trace and wavecal dict corresponding to order.

    Parameters
    ----------
    refmodel : datamodel
        The reference file holding traces and wavelength calibration
        models, under `refmodel.traces` and `refmodel.wavecal_models`
    order : int
        The spectral order to find trace and wavecal model indices for.

    Returns
    -------
    int
        The index to provide the reference file lists of traces and wavecal
        models to retrieve the arrays for the desired spectral order
    """
    if order not in [1, 2]:
        error_message = f"Order {order} is not supported at this time."
        log.error(error_message)
        raise ValueError(error_message)

    for i, entry in enumerate(refmodel.traces):
        if entry.spectral_order == order:
            return i

    log.warning("Order not found in reference file trace list.")
    return -1


def get_soss_traces(refmodel, pwcpos, order, subarray):
    """
    Generate the traces given a pupil wheel position.

    This is the primary method for generating the gr700xd trace position given a
    pupil wheel position angle provided in the FITS header under keyword
    PWCPOS. The traces for a given spectral order are found by performing a
    rotation transformation using the reference trace positions at the commanded
    PWCPOS=245.76 degrees. This method yields sub-pixel performance and will be
    improved upon in later iterations as more NIRISS/SOSS observations become
    available.

    Parameters
    ----------
    refmodel : PastasossModel
        The reference file datamodel.
    pwcpos : float
        The pupil wheel positions angle provided in the FITS header under
        keyword PWCPOS.
    order : str or int
        The spectral order for which a trace is computed.
        Order 3 is currently unsupported.
    subarray : str
        Name of subarray in use, typically 'SUBSTRIP96' or 'SUBSTRIP256'.

    Returns
    -------
    order : str
        The spectral order for which a trace is computed.
    x_new, y_new : Tuple[np.ndarray, np.ndarray]]
        If `order` is '1', a tuple of the x and y coordinates of the rotated
        points for the first spectral order.
        If `order` is '2', a tuple of the x and y coordinates of the rotated
        points for the second spectral order.
    wavelengths : np.ndarray
        The wavelengths associated with the rotated points.

    Raises
    ------
    ValueError
        If `order` is not in ['1', '2'].
    """
    spectral_order_index = _find_spectral_order_index(refmodel, int(order))

    # reference trace data
    x, y = refmodel.traces[spectral_order_index].trace.T.copy()
    origin = (
        refmodel.traces[spectral_order_index].pivot_x,
        refmodel.traces[spectral_order_index].pivot_y,
    )

    # Offset for SUBSTRIP96
    if subarray == "SUBSTRIP96":
        y -= 10
    # rotated reference trace
    x_new, y_new = _rotate(x, y, pwcpos - refmodel.meta.pwcpos_cmd, origin)

    # wavelength associated to trace at given pwcpos value
    wavelengths = _get_wavelengths(refmodel, x_new, pwcpos, int(order))

    return order, x_new, y_new, wavelengths


def _extrapolate_to_wavegrid(w_grid, wavelength, quantity):
    """
    Extrapolate quantities on the right and the left of a given array of quantity.

    Parameters
    ----------
    w_grid : sequence
        The wavelength grid to interpolate onto
    wavelength : sequence
        The native wavelength values of the data
    quantity : sequence
        The data to interpolate

    Returns
    -------
    Array
        The interpolated quantities
    """
    sort_i = np.argsort(wavelength)
    q = quantity[sort_i]
    w = wavelength[sort_i]

    # Determine the slope on the right of the array
    slope_right = (q[-1] - q[-2]) / (w[-1] - w[-2])
    # extrapolate at wavelengths larger than the max on the right
    indright = np.where(w_grid > w[-1])[0]
    q_right = q[-1] + (w_grid[indright] - w[-1]) * slope_right
    # Determine the slope on the left of the array
    slope_left = (q[1] - q[0]) / (w[1] - w[0])
    # extrapolate at wavelengths smaller than the min on the left
    indleft = np.where(w_grid < w[0])[0]
    q_left = q[0] + (w_grid[indleft] - w[0]) * slope_left
    # Construct and extrapolated array of the quantity
    w = np.concatenate((w_grid[indleft], w, w_grid[indright]))
    q = np.concatenate((q_left, q, q_right))

    # resample at the w_grid everywhere
    return np.interp(w_grid, w, q)


def _calc_2d_wave_map(wave_grid, x_dms, y_dms, tilt, oversample=2, padding=0, maxiter=5, dtol=1e-2):
    """
    Compute the 2D wavelength map on the detector.

    Parameters
    ----------
    wave_grid : sequence
        The wavelength corresponding to the x_dms, y_dms, and tilt values.
    x_dms : sequence
        The trace x position on the detector in DMS coordinates.
    y_dms : sequence
        The trace y position on the detector in DMS coordinates.
    tilt : sequence
        The trace tilt angle in degrees.
    oversample : int
        The oversampling factor of the input coordinates.
    padding : int
        The native pixel padding around the edge of the detector.
    maxiter : int
        The maximum number of iterations used when solving for the wavelength at each pixel.
    dtol : float
        The tolerance of the iterative solution in pixels.

    Returns
    -------
    Array
        An array containing the wavelength at each pixel on the detector.
    """
    os = np.copy(oversample)
    xpad = np.copy(padding)
    ypad = np.copy(padding)

    # No need to compute wavelengths across the entire detector,
    # slightly larger than SUBSTRIP256 will do.
    dimx, dimy = SOSS_XDIM, SOSS_YDIM
    y_dms = y_dms + (dimy - SOSS_XDIM)  # Adjust y-coordinate to area of interest.

    # Generate the oversampled grid of pixel coordinates.
    x_vec = np.arange((dimx + 2 * xpad) * os) / os - (os - 1) / (2 * os) - xpad
    y_vec = np.arange((dimy + 2 * ypad) * os) / os - (os - 1) / (2 * os) - ypad
    x_grid, y_grid = np.meshgrid(x_vec, y_vec)

    # Iteratively compute the wavelength at each pixel.
    delta_x = 0.0  # A shift in x represents a shift in wavelength.
    for _niter in range(maxiter):
        # Assume all y have same wavelength.
        wave_iterated = np.interp(
            x_grid - delta_x, x_dms[::-1], wave_grid[::-1]
        )  # Invert arrays to get increasing x.

        # Compute the tilt angle at the wavelengths.
        tilt_tmp = np.interp(wave_iterated, wave_grid, tilt)

        # Compute the trace position at the wavelengths.
        x_estimate = np.interp(wave_iterated, wave_grid, x_dms)
        y_estimate = np.interp(wave_iterated, wave_grid, y_dms)

        # Project that back to pixel coordinates.
        x_iterated = x_estimate + (y_grid - y_estimate) * np.tan(np.deg2rad(tilt_tmp))

        # Measure error between requested and iterated position.
        delta_x = delta_x + (x_iterated - x_grid)

        # If the desired precision has been reached end iterations.
        if np.all(np.abs(x_iterated - x_grid) < dtol):
            break

    # Evaluate the final wavelength map, this time setting out-of-bounds values to NaN.
    wave_map_2d = np.interp(
        x_grid - delta_x, x_dms[::-1], wave_grid[::-1], left=np.nan, right=np.nan
    )

    # Extend to full detector size.
    tmp = np.full((os * (dimx + 2 * xpad), os * (dimx + 2 * xpad)), fill_value=np.nan)
    tmp[-os * (dimy + 2 * ypad) :] = wave_map_2d
    wave_map_2d = tmp

    return wave_map_2d


def get_soss_wavemaps(refmodel, pwcpos, subarray, padding=False, padsize=0, spectraces=False):
    """
    Generate order 1 and 2 2D wavemaps from the rotated SOSS trace positions.

    Parameters
    ----------
    pwcpos : float
        The pupil wheel position
    subarray : str
        The subarray name, ['FULL', 'SUBSTRIP256', 'SUBSTRIP96']
    padding : bool
        Include padding on map edges (only needed for reference files)
    padsize : int
        The size of the padding to include on each side
    spectraces : bool
        Return the interpolated spectraces as well

    Returns
    -------
    Array, Array
        The 2D wavemaps and corresponding 1D spectraces
    """
    _, order1_x, order1_y, order1_wl = get_soss_traces(
        refmodel, pwcpos, order="1", subarray=subarray
    )
    _, order2_x, order2_y, order2_wl = get_soss_traces(
        refmodel, pwcpos, order="2", subarray=subarray
    )

    # Make wavemap from trace center wavelengths, padding to shape (296, 2088)
    wavemin = WAVEMAP_WLMIN
    wavemax = WAVEMAP_WLMAX
    nwave = WAVEMAP_NWL
    wave_grid = np.linspace(wavemin, wavemax, nwave)

    # Extrapolate wavelengths for order 1 trace
    xtrace_order1 = _extrapolate_to_wavegrid(wave_grid, order1_wl, order1_x)
    ytrace_order1 = _extrapolate_to_wavegrid(wave_grid, order1_wl, order1_y)
    spectrace_1 = np.array([xtrace_order1, ytrace_order1, wave_grid])

    # Set cutoff for order 2 where it runs off the detector
    o2_cutoff = XTRACE_ORD2_LEN
    w_o2_tmp = order2_wl[:o2_cutoff]
    # Subtract 8 from FULL width to avoid reference pixels
    w_o2 = np.zeros(SOSS_XDIM - 8) * np.nan
    w_o2[:o2_cutoff] = w_o2_tmp
    y_o2_tmp = order2_y[:o2_cutoff]
    y_o2 = np.zeros(SOSS_XDIM - 8) * np.nan
    y_o2[:o2_cutoff] = y_o2_tmp
    x_o2 = np.copy(order1_x)

    # Fill for column > 1400 with linear extrapolation
    m = w_o2[o2_cutoff - 1] - w_o2[o2_cutoff - 2]
    dx = np.arange(SOSS_XDIM - 8 - o2_cutoff) + 1
    w_o2[o2_cutoff:] = w_o2[o2_cutoff - 1] + m * dx
    m = y_o2[o2_cutoff - 1] - y_o2[o2_cutoff - 2]
    dx = np.arange(SOSS_XDIM - 8 - o2_cutoff) + 1
    y_o2[o2_cutoff:] = y_o2[o2_cutoff - 1] + m * dx

    # Extrapolate wavelengths for order 2 trace
    xtrace_order2 = _extrapolate_to_wavegrid(wave_grid, w_o2, x_o2)
    ytrace_order2 = _extrapolate_to_wavegrid(wave_grid, w_o2, y_o2)
    spectrace_2 = np.array([xtrace_order2, ytrace_order2, wave_grid])

    # Make wavemap from wavelength solution for order 1
    wavemap_1 = _calc_2d_wave_map(
        wave_grid,
        xtrace_order1,
        ytrace_order1,
        np.zeros_like(xtrace_order1),
        oversample=1,
        padding=padsize,
    )

    # Make wavemap from wavelength solution for order 2
    wavemap_2 = _calc_2d_wave_map(
        wave_grid,
        xtrace_order2,
        ytrace_order2,
        np.zeros_like(xtrace_order2),
        oversample=1,
        padding=padsize,
    )

    # Extrapolate wavemap to FULL frame
    wavemap_1[: SUBARRAY_YMIN - padsize, :] = wavemap_1[SUBARRAY_YMIN - padsize]
    wavemap_2[: SUBARRAY_YMIN - padsize, :] = wavemap_2[SUBARRAY_YMIN - padsize]

    # Trim to subarray
    if subarray == "SUBSTRIP256":
        wavemap_1 = wavemap_1[SUBARRAY_YMIN - padsize : SOSS_XDIM + padsize, :]
        wavemap_2 = wavemap_2[SUBARRAY_YMIN - padsize : SOSS_XDIM + padsize, :]
    if subarray == "SUBSTRIP96":
        wavemap_1 = wavemap_1[SUBARRAY_YMIN - padsize : SUBARRAY_YMIN + 96 + padsize, :]
        wavemap_2 = wavemap_2[SUBARRAY_YMIN - padsize : SUBARRAY_YMIN + 96 + padsize, :]

    # Remove padding if necessary
    if not padding and padsize != 0:
        wavemap_1 = wavemap_1[padsize:-padsize, padsize:-padsize]
        wavemap_2 = wavemap_2[padsize:-padsize, padsize:-padsize]

    if spectraces:
        return np.array([wavemap_1, wavemap_2]), np.array([spectrace_1, spectrace_2])
    return np.array([wavemap_1, wavemap_2])
