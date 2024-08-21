# Module to predict SOSS trace positions for a given spectral order(s)
# given the  a pupil wheel position angle taken from the "PWCPOS" fits
# header keyword.

from dataclasses import dataclass
from typing import Tuple

import numpy as np
from scipy.interpolate import interp1d
from pkg_resources import resource_filename

from pastasoss.wavecal import get_wavecal_meta_for_spectral_order
from pastasoss.wavecal import get_wavelengths

PWCPOS_CMD = 245.7600  # Commanded PWCPOS for the GR700XD


# TODO: order 3 currently unsupport ATM. Will be support in the future: TBD
REFERENCE_TRACE_FILES = {
    "order1": "jwst_niriss_gr700xd_order1_trace_refmodel.txt",
    "order2": "jwst_niriss_gr700xd_order2_trace_refmodel_002.txt",
    # order 3 currently unsupport ATM. Will be support in the future: TBD
}

REFERENCE_WAVECAL_MODELS = {
    "order1": resource_filename(
        __name__, "data/jwst_niriss_gr700xd_wavelength_model_order1.json"
    ),
    "order2": resource_filename(
        __name__, "data/jwst_niriss_gr700xd_wavelength_model_order2_002.json"
    ),
}


@dataclass
class TraceModel:
    order: str
    x: np.ndarray
    y: np.ndarray
    wavelength: np.ndarray


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


def get_reference_trace(
    file: str,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Load in the reference trace positions given a file associated with a given
    spectral order.

    Parameters
    ----------
    file : str
        The path to the file containing the reference trace positions.

    Returns
    -------
    Tuple[np.ndarray, np.ndarray, np.ndarray]
        A tuple containing the x, y positions and origin of the reference
        traces.

    Examples
    --------
    >>> x, y, origin = get_reference_traces_positions('ref_filename.txt')
    """
    filepath = resource_filename(__name__, f"data/{file}")
    traces = np.loadtxt(filepath)
    origin = traces[0]
    x = traces[1:, 0]
    y = traces[1:, 1]
    return x, y, origin


def get_soss_traces(
    pwcpos: float, order: str = "123", subarray: str = 'SUBSTRIP256', interp: bool = True
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
    subarray : str
        The subarray being used, ['SUBSTRIP96', 'SUBSTRIP256']
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
        return [get_soss_traces(pwcpos, m, subarray) for m in order]

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

    # Offset for SUBSTRIP96
    if subarray == 'SUBSTRIP96':
        y -= 10

    # rotated reference trace
    x_new, y_new = rotate(x, y, pwcpos - PWCPOS_CMD, origin, interp=interp)

    # wavelength associated to trace at given pwcpos value
    wavelengths = get_wavelengths(x_new, pwcpos, wave_cal_model_meta)

    # return x_new, y_new, wavelengths
    return TraceModel(order, x_new, y_new, wavelengths)
