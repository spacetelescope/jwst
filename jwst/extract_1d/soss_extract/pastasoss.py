from dataclasses import dataclass
from functools import partial
from typing import Any, Dict, Tuple
import logging
import numpy as np
from pkg_resources import resource_filename
from scipy.interpolate import interp1d
from ...datamodels import SossTraceModel

log = logging.getLogger(__name__)


def get_wavelengths(
    refmodel: SossTraceModel, x: np.ndarray, pwcpos: float, order: int
) -> np.ndarray:
    """Get the associated wavelength values for a given spectral order"""
    if order == 1:
        wavelengths = wavecal_model_order1_poly(x, pwcpos, wavecal_meta)
    elif order == 2:
        wavelengths = wavecal_model_order2_poly(x, pwcpos, wavecal_meta)
    else:
        error_message = f"Order {order} is not supported at this time."
        log.error(error_message)
        raise ValueError(error_message)

    return wavelengths


def min_max_scaler(x, x_min, x_max):
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
    # scaling the input x values
    x_scaled = (x - x_min) / (x_max - x_min)
    return x_scaled


def wavecal_model_order1_poly(refmodel, x, pwcpos):
    """compute order 1 wavelengths"""
    x_scaler = partial(
        min_max_scaler,
        **{
            "x_min": refmodel.wavecal_models.scale_extents[0][0],
            "x_max": refmodel.wavecal_models.scale_extents[1][0],
        },
    )

    pwcpos_offset_scaler = partial(
        min_max_scaler,
        **{
            "x_min": refmodel.wavecal_models.scale_extents[0][1],
            "x_max": refmodel.wavecal_models.scale_extents[1][1],
        },
    )

    def get_poly_features(x: np.array, offset: np.array) -> np.ndarray:
        """polynomial features for the order 1 wavecal model"""
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
    offset = np.ones_like(x) * (pwcpos - PWCPOS_CMD)
    offset_scaled = pwcpos_offset_scaler(offset)

    # polynomial features
    poly_features = get_poly_features(x_scaled, offset_scaled)
    wavelengths = coef[0] + coef[1:] @ poly_features

    return wavelengths


def wavecal_model_order2_poly(x, pwcpos, wavecal_meta: NIRISS_GR700XD_WAVECAL_META):
    """compute order 2 wavelengths"""
    x_scaler = partial(
        min_max_scaler,
        **{
            "x_min": wavecal_meta.scaler_data_min_[0],
            "x_max": wavecal_meta.scaler_data_max_[0],
        },
    )

    pwcpos_offset_scaler = partial(
        min_max_scaler,
        **{
            "x_min": wavecal_meta.scaler_data_min_[1],
            "x_max": wavecal_meta.scaler_data_max_[1],
        },
    )

    def get_poly_features(x: np.array, offset: np.array) -> np.ndarray:
        """Polynomial features for the order 2 wavecal model"""
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
    coef = wavecal_meta.coefficients

    # get pixel columns and then scaled
    x_scaled = x_scaler(x)

    # offset
    # offset = np.ones_like(x) * (pwcpos - PWCPOS_CMD)
    # this will need to get changed later...
    offset = np.ones_like(x) * pwcpos
    offset_scaled = pwcpos_offset_scaler(offset)

    # polynomial features
    poly_features = get_poly_features(x_scaled, offset_scaled)
    wavelengths = coef[0] + coef[1:] @ poly_features

    return wavelengths


def rotate(
    x: np.ndarray,
    y: np.ndarray,
    angle: float,
    origin: Tuple[float, float] = (0, 0),
    interp: bool = True,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Applies a rotation transformation to a set of 2D points.

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
    interp : bool, optional
        Whether to interpolate the rotated positions onto the original x-pixel
        column values. Default is True.

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        The x and y coordinates of the rotated points.

    Examples
    --------
    >>> x = np.array([0, 1, 2, 3])
    >>> y = np.array([0, 1, 2, 3])
    >>> x_rot, y_rot = rotate(x, y, 90)
    """

    # shift to rotate about center
    xy_center = np.atleast_2d(origin).T
    xy = np.vstack([x, y])

    # Rotation transform matrix
    radians = np.radians(angle)
    c, s = np.cos(radians), np.sin(radians)
    R = np.array([[c, -s], [s, c]])

    # apply transformation
    x_new, y_new = R @ (xy - xy_center) + xy_center

    # interpolate rotated positions onto x-pixel column values (default)
    if interp:
        # interpolate new coordinates onto original x values and mask values
        # outside of the domain of the image 0<=x<=2047 and 0<=y<=255.
        y_new = interp1d(x_new, y_new, fill_value="extrapolate")(x)
        mask = np.where(y_new <= 255.0)
        x = x[mask]
        y_new = y_new[mask]
        return x, y_new

    return x_new, y_new


def get_soss_traces(
    refmodel: SossTraceModel, pwcpos: float,
        order: str = "12", interp: bool = True
) -> Tuple[np.ndarray, np.ndarray]:
    """
    This is the primary method for generate the gr700xd trace position given a
    pupil wheel positions angle provided in the FITS header under keyword
    PWCPOS. The traces for a given spectral order are found by perform a
    rotation transformation using the refence trace positions at the commanded
    PWCPOS=245.76 degrees. This methods yield sub-pixel performance and will be
    improved upon in later interations as more NIRISS/SOSS observations become
    available.

    Parameters
    ----------
    pwcpos : float
        The pupil wheel positions angle provided in the FITS header under
        keyword PWCPOS.
    order : str, optional
        The spectral order to compute the new traces for. Default is '123'.
        Support for order 3 will be added at a later date.
    interp : bool, optional
        Whether to interpolate the rotated positions onto the original x-pixel
        column values. Default is True.

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]]
    If `order` is '1', a tuple of the x and y coordinates of the rotated
    points for the first spectral order.
    If `order` is '2', a tuple of the x and y coordinates of the rotated
    points for the second spectral order.
    If `order` is '3' or a combination of '1', '2', and '3', a list of
    tuples of the x and y coordinates of the rotated points for each
    spectral order.

    Raises
    ------
    ValueError
        If `order` is not '1', '2', '3', or a combination of '1', '2', and '3'.

    Examples
    --------
    >>> x_new, y_new = get_trace_from_reference_transform(2.3)
    """

    norders = len(order)

    if norders > 3:
        raise ValueError("order must be: 1,2, 3.")
    if norders > 1:
        # recursively compute the new traces for each order
        return [get_soss_traces(pwcpos, m) for m in order]

    # This might be an alternative way of writing this
    # if 'order'+order in REFERENCE_TRACE_FILES.keys():
    #     ref_file = REFERENCE_TRACE_FILES["order"+order]
    # This section can definite be refactored in a later version
    elif order == "1":
        ref_trace_file = REFERENCE_TRACE_FILES["order1"]
        wave_cal_model_meta = get_wavecal_meta_for_spectral_order("order1")

    elif order == "2":
        ref_trace_file = REFERENCE_TRACE_FILES["order2"]
        wave_cal_model_meta = get_wavecal_meta_for_spectral_order("order2")

    elif order == "3":
        print("The software currently does not support order 3 at this time.")
        return None

    # reference trace data
    x, y, origin = get_reference_trace(ref_trace_file)

    # rotated reference trace
    x_new, y_new = rotate(x, y, pwcpos - PWCPOS_CMD, origin, interp=interp)

    # wavelength associated to trace at given pwcpos value
    wavelengths = get_wavelengths(x_new, pwcpos, wave_cal_model_meta)

    # return x_new, y_new, wavelengths
    return TraceModel(order, x_new, y_new, wavelengths)