import logging
from datetime import datetime
from functools import partial

import numpy as np
import stdatamodels.jwst.datamodels as dm
from astropy.modeling.polynomial import Polynomial2D
from scipy.interpolate import interp1d
from stpipe import crds_client

log = logging.getLogger(__name__)

SOSS_XDIM = 2048
SOSS_YDIM = 300
CUTOFFS = [SOSS_XDIM, 1783, 1134]  # max pixel x-index to consider for a given order
WAVEMAP_WLMIN = 0.5  # min wavelength for pastasoss 1-d wavelength grid
WAVEMAP_WLMAX = 5.5  # max wavelength for pastasoss 1-d wavelength grid
WAVEMAP_NWL = 5001  # number of wavelengths for pastasoss 1-d wavelength grid
SUBARRAY_YMIN = 2048 - 256  # pixel y-value defining the start of the subarray
PWCPOS_BOUNDS = (245.79 - 0.25, 245.79 + 0.25)  # reasonable PWC position limits
DEFAULT_CRDS_PARAMS = {
    "meta.instrument.name": "NIRISS",
    "meta.observation.date": datetime.today().strftime("%Y-%m-%d"),
    "meta.observation.time": datetime.today().strftime("%H:%M:%S.%f"),
    "meta.instrument.detector": "NIS",
    "meta.instrument.filter": "CLEAR",
    "meta.exposure.type": "NIS_SOSS",
}

__all__ = ["get_soss_traces", "get_soss_wavemaps", "retrieve_default_pastasoss_model"]


def _verify_requested_orders(orders_requested, refmodel_orders):
    """
    Verify that the requested orders are valid.

    Parameters
    ----------
    orders_requested : list
        A list of the spectral orders requested for extraction.
    refmodel_orders : list
        A list of the spectral orders available in the reference model.

    Returns
    -------
    list
        The validated list of requested orders.

    Raises
    ------
    ValueError
        If any of the requested orders are not available in the reference model.
    """
    orders_requested = np.array(orders_requested)
    refmodel_orders = np.array(refmodel_orders)
    not_in = np.isin(orders_requested, refmodel_orders, invert=True)
    if np.all(not_in):
        raise ValueError(
            "None of the requested orders are available in the PASTASOSS reference file: "
            f"{orders_requested}. Defined orders are {refmodel_orders}."
        )
    if np.any(not_in):
        orders_requested = orders_requested[~not_in]
        log.warning(
            "Some requested orders were not found in reference model. Skipping those orders "
            f"and proceeding with orders {orders_requested}."
        )
    return orders_requested.tolist()


def _convert_refmodel_poly_to_astropy(coefficients):
    """
    Reorder reference file 2-D polynomial coefficients to create Astropy polynomial.

    Ordering in reference files is expected to be
    C00 + C10 * x + C01 * y + C20 * x^2 + C11 * x * y + C02 * y^2 +
    C30 * x^3 + C21 * x^2 * y + C12 * x * y^2 + C03 * y^3 ...

    Astropy ordering is
    C00 + C10 * x + C20 * x^2 ...
    + C01 * y + C02 * y^2 + ...
    + C11 * x * y + C12 * x^2 * y + C13 * x^3 * y ...

    Parameters
    ----------
    coefficients : list
        A list of polynomial coefficients in the reference file format.

    Returns
    -------
    Polynomial2D
        The Astropy 2-D polynomial representation of the input coefficients.
    """
    # figure out the degree from the length by inverting triangle number formula
    degree = int(np.sqrt(2 * len(coefficients))) - 1
    poly = Polynomial2D(degree=degree)
    coeff_names = poly.param_names
    for name in coeff_names:
        # compute index in coefficients corresponding to that order
        xord, yord = (int(n) for n in name.strip("c").split("_"))
        ord_sum = xord + yord
        min_idx_ord_sum = ord_sum * (ord_sum + 1) // 2  # triangle number
        idx = min_idx_ord_sum + yord
        coeff = coefficients[idx]
        setattr(poly, name, coeff)

    return poly


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
    order_idx = _find_spectral_order_index(refmodel, order)
    x_scaler = partial(
        _min_max_scaler,
        **{
            "x_min": refmodel.wavecal_models[order_idx].scale_extents[0][0],
            "x_max": refmodel.wavecal_models[order_idx].scale_extents[1][0],
        },
    )

    pwcpos_offset_scaler = partial(
        _min_max_scaler,
        **{
            "x_min": refmodel.wavecal_models[order_idx].scale_extents[0][1],
            "x_max": refmodel.wavecal_models[order_idx].scale_extents[1][1],
        },
    )

    # extract model weights and intercept
    coef = refmodel.wavecal_models[order_idx].coefficients
    poly = _convert_refmodel_poly_to_astropy(coef)

    # scale x and pwcpos offset to be between 0 and 1
    x_scaled = x_scaler(x)
    if int(order) == 2:
        # Reference file has a bug for order 2 where the scale_extents are not
        # defined as offsets: they have not had the commanded position subtracted
        offset = np.ones_like(x) * pwcpos
    else:
        offset = np.ones_like(x) * (pwcpos - refmodel.meta.pwcpos_cmd)
    offset_scaled = pwcpos_offset_scaler(offset)

    wavelengths = poly(x_scaled, offset_scaled)

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


def _rotate(x, y, angle, origin=(0, 0)):
    """
    Apply a rotation transformation to a set of 2D points.

    Parameters
    ----------
    x : ndarray
        The x-coordinates of the points to be transformed.
    y : ndarray
        The y-coordinates of the points to be transformed.
    angle : float
        The angle (in degrees) by which to rotate the points.
    origin : Tuple[float, float], optional
        The point about which to rotate the points. Default is (0, 0).

    Returns
    -------
    tuple[ndarray, ndarray]
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
    if order not in [1, 2, 3]:
        error_message = f"Order {order} is not supported."
        log.error(error_message)
        raise ValueError(error_message)

    for i, entry in enumerate(refmodel.traces):
        if entry.spectral_order == order:
            return i

    log.warning("Order not found in reference file trace list.")
    return -1


def get_soss_traces(pwcpos, order, subarray="SUBSTRIP256", refmodel=None):
    """
    Get the SOSS traces for a given input model and spectral order.

    Parameters
    ----------
    pwcpos : float
        The pupil wheel position angle provided in the FITS header under keyword PWCPOS.
        Values are expected to be within +/- 0.25 degrees of the commanded position
        (245.76 degrees).
    order : int
        The spectral order for which to retrieve the traces.
    subarray : str
        Name of subarray in use, typically 'SUBSTRIP96' or 'SUBSTRIP256'.
    refmodel : PastasossModel, optional
        The reference model for the SOSS extraction. If not set, it will be fetched
        from CRDS.

    Returns
    -------
    order : str
        The spectral order for which a trace is computed.
    x : ndarray
        The x coordinates of the rotated points.
    y : ndarray
        The y coordinates of the rotated points.
    wavelengths : ndarray
        The wavelengths associated with the rotated points.
    """
    if refmodel is None:
        refmodel = retrieve_default_pastasoss_model()

    pwcpos_is_valid = _check_pwcpos_bounds(pwcpos)
    if not pwcpos_is_valid:
        raise ValueError(f"PWC position {pwcpos} is outside bounds ({PWCPOS_BOUNDS}).")

    return _get_soss_traces(refmodel, pwcpos, order, subarray)


def retrieve_default_pastasoss_model():
    """
    Retrieve the default PastasossModel reference file.

    This function fetches the default PastasossModel reference file from CRDS.
    It is used when no specific reference model is provided.

    Returns
    -------
    PastasossModel
        The default PastasossModel reference file.
    """
    ref_name = crds_client.get_reference_file(DEFAULT_CRDS_PARAMS, "pastasoss", "jwst")
    ref_file = crds_client.check_reference_open(ref_name)
    return dm.PastasossModel(ref_file)


def _get_soss_traces(refmodel, pwcpos, order, subarray):
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
    subarray : str
        Name of subarray in use, typically 'SUBSTRIP96' or 'SUBSTRIP256'.

    Returns
    -------
    order : str
        The spectral order for which a trace is computed.
    x : ndarray
        The x coordinates of the rotated points.
    y : ndarray
        The y coordinates of the rotated points.
    wavelengths : ndarray
        The wavelengths associated with the rotated points.
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


def get_soss_wavemaps(
    pwcpos,
    subarray="SUBSTRIP256",
    refmodel=None,
    padsize=None,
    spectraces=False,
    orders_requested=None,
):
    """
    Get the SOSS wavelength maps and (optionally) spectraces.

    Parameters
    ----------
    pwcpos : float
        The pupil wheel position angle, e.g. as provided in the FITS header under keyword PWCPOS.
        Values are expected to be within +/- 0.25 degrees of the commanded position
        (245.76 degrees).
    subarray : str, optional
        The subarray name, one of 'SUBSTRIP256', 'SUBSTRIP96', or 'FULL'.
    refmodel : PastasossModel, optional
        The reference model for the SOSS extraction. If not set, it will be fetched
        from CRDS.
    padsize : int, optional
        The padding to apply to the wavelength maps.
    spectraces : bool, optional
        If True, return the interpolated spectraces as well.
    orders_requested : list
        A list of the spectral orders requested for extraction.
        If None, all orders in the reference file will be used.

    Returns
    -------
    wavemaps : ndarray
        The 2D wavemaps. Will have shape (n_orders, array_x, array_y) with orders 1, 2, etc. being
        the elements of the first dimension. Wavemaps for all orders defined in the reference file
        will be returned.
    traces : ndarray, optional
        The corresponding 1D traces (if ``spectraces`` is True).
    """
    if refmodel is None:
        refmodel = retrieve_default_pastasoss_model()

    refmodel_orders = [int(trace.spectral_order) for trace in refmodel.traces]
    if orders_requested is None:
        orders_requested = refmodel_orders
    else:
        orders_requested = _verify_requested_orders(orders_requested, refmodel_orders)

    pwcpos_is_valid = _check_pwcpos_bounds(pwcpos)
    if not pwcpos_is_valid:
        raise ValueError(f"PWC position {pwcpos} is outside bounds ({PWCPOS_BOUNDS}).")

    if padsize is None:
        padsize = getattr(refmodel.traces[0], "padding", 0)
    if padsize > 0:
        do_padding = True
    else:
        do_padding = False

    # Make wavemap from trace center wavelengths, padding to shape (296, 2088)
    wavemin = WAVEMAP_WLMIN
    wavemax = WAVEMAP_WLMAX
    nwave = WAVEMAP_NWL
    wave_grid = np.linspace(wavemin, wavemax, nwave)

    wavemaps = []
    traces = []
    for order in orders_requested:
        _, x, y, wl = _get_soss_traces(refmodel, pwcpos, order=str(order), subarray=subarray)
        xtrace = _extrapolate_to_wavegrid(wave_grid, wl, x)
        ytrace = _extrapolate_to_wavegrid(wave_grid, wl, y)
        spectrace = np.array([xtrace, ytrace, wave_grid])

        # Make wavemap from wavelength solution
        wavemap = _calc_2d_wave_map(
            wave_grid,
            xtrace,
            ytrace,
            np.zeros_like(xtrace),
            oversample=1,
            padding=padsize,
        )
        # Extrapolate wavemap to FULL frame
        wavemap[: SUBARRAY_YMIN - padsize, :] = wavemap[SUBARRAY_YMIN - padsize]

        # Trim to subarray
        if subarray == "SUBSTRIP256":
            wavemap = wavemap[SUBARRAY_YMIN - padsize : SOSS_XDIM + padsize, :]
        if subarray == "SUBSTRIP96":
            wavemap = wavemap[SUBARRAY_YMIN - padsize : SUBARRAY_YMIN + 96 + padsize, :]

        # remove padding if necessary
        if not do_padding and padsize != 0:
            wavemap = wavemap[padsize:-padsize, padsize:-padsize]
        wavemaps.append(wavemap)
        traces.append(spectrace)

    # Combine wavemaps and spectraces into ndarray output
    if spectraces:
        return np.array(wavemaps), np.array(traces)
    return np.array(wavemaps)


def _check_pwcpos_bounds(pwcpos):
    """
    Check if the provided PWC position is within the bounds.

    Parameters
    ----------
    pwcpos : float
        The pupil wheel position angle.

    Returns
    -------
    bool
        True if the PWC position is within bounds, False otherwise.
    """
    return PWCPOS_BOUNDS[0] <= pwcpos <= PWCPOS_BOUNDS[1]
