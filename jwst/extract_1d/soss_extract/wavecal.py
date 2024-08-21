# Wavecal module for handling the wavelength calibration model. For now we
# will use the reference json files that has all the model information. Ideally
# we would like to use just the train model for production but for now this is
# this is fine.

import json
from dataclasses import dataclass
from functools import partial
from typing import Any, Dict

import numpy as np
import numpy.typing as npt
from pkg_resources import resource_filename

# explicitly defining the commanded position here as well (this is temp)
PWCPOS_CMD = 245.76

# # TODO: order 3 currently unsupport ATM. Will be support in the future: TBD
# REFERENCE_WAVECAL_MODELS = {
#     "order1": resource_filename(
#         __name__, "data/jwst_niriss_gr700xd_wavelength_model_order1.json"
#     ),
#     "order2": resource_filename(
#         __name__, "data/jwst_niriss_gr700xd_wavelength_model_order2_002.json"
#     ),
# }


# dataclass object to store the metadata from a wavecal model
@dataclass
class NIRISS_GR700XD_WAVECAL_META:
    """
    Datamodel object container for the wavecal meta data in the json reference
    file.
    """

    order: str
    coefficients: npt.NDArray[np.float64]
    intercept: npt.NDArray[np.float64]
    scaler_data_min_: float
    scaler_data_max_: float
    poly_degree: int


def load_wavecal_model(filename: str) -> Dict[str, Any]:
    """
    Load a wavecalibration model from a JSON file.

    Parameters
    ----------
    filename : str
        The path to the JSON file containing the wavecalibration model.

    Returns
    -------
    dict
        A dictionary representing the loaded wavecalibration model.

    Notes
    -----
    This function reads the specified JSON file and returns its contents as a
    dictionary, which typically contains information about a wavecalibration
    model used in data analysis.
    """
    with open(filename, "r") as file:
        wavecal_model = json.load(file)
    return wavecal_model


def get_wavecal_meta_for_spectral_order(
    order: str,
) -> NIRISS_GR700XD_WAVECAL_META:
    """
    Retrieve wavecalibration model metadata for a specific spectral order.

    Parameters
    ----------
    order : str
        The spectral order for which wavecalibration metadata is requested.
        Valid options are 'order 1', 'order 2', or 'order 3'.

    Returns
    -------
    NIRISS_GR700XD_WAVECAL_META
        An object containing wavecalibration metadata, including order, model
        coefficients, intercept, and scaler data bounds.

    Raises
    ------
    ValueError
        If the provided 'order' is not one of the valid options.

    Notes
    -----
    This function retrieves wavecalibration metadata for the specified spectral
    order. It first checks if the 'order' is valid and then loads the
    corresponding wavecalibration model. The model's coefficients, intercept,
    and scaler data bounds are extracted and returned as part of the metadata
    object.
    """
    # get the reference wavecal file name
    valid_orders = ['order 1', 'order 2']
    if order not in valid_orders:
        raise ValueError(f"valid orders are: {valid_orders}.")

    # get the appropiate reference file name given the order
    reference_filename = REFERENCE_WAVECAL_MODELS[order]

    # load in the model
    wavecal_model = load_wavecal_model(reference_filename)

    # model coefficients
    poly_degree = wavecal_model["model"]["poly_deg"]
    coefficients = wavecal_model["model"]["coef"]
    intercept = wavecal_model["model"]["intercept"]

    # info for scaling inputs
    scaler_data_min_ = wavecal_model["model"]["scaler"]["data_min_"]
    scaler_data_max_ = wavecal_model["model"]["scaler"]["data_max_"]

    return NIRISS_GR700XD_WAVECAL_META(
        order, coefficients, intercept, scaler_data_min_, scaler_data_max_, poly_degree
    )


def get_wavelengths(
    x: np.ndarray, pwcpos: float, wavecal_meta: NIRISS_GR700XD_WAVECAL_META
) -> np.ndarray:
    """Get the associated wavelength values for a given spectral order"""
    if wavecal_meta.order == "order1":
        wavelengths = wavecal_model_order1_poly(x, pwcpos, wavecal_meta)
    elif wavecal_meta.order == "order2":
        # raise NotImplementedError("Order 2 not implemented")
        wavelengths = wavecal_model_order2_poly(x, pwcpos, wavecal_meta)
    elif wavecal_meta.order == "order3":
        raise ValueError("Order 3 not supported at this time")
    else:
        raise ValueError("not a valid order")

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


def wavecal_model_order1_poly(x, pwcpos, wavecal_meta: NIRISS_GR700XD_WAVECAL_META):
    """compute order 1 wavelengths"""
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
    coef = wavecal_meta.coefficients
    intercept = wavecal_meta.intercept

    # get pixel columns and then scaled
    x_scaled = x_scaler(x)

    # offset
    offset = np.ones_like(x) * (pwcpos - PWCPOS_CMD)
    offset_scaled = pwcpos_offset_scaler(offset)

    # polynomial features
    poly_features = get_poly_features(x_scaled, offset_scaled)
    wavelengths = coef @ poly_features + intercept

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
    intercept = wavecal_meta.intercept

    # get pixel columns and then scaled
    x_scaled = x_scaler(x)

    # offset
    # offset = np.ones_like(x) * (pwcpos - PWCPOS_CMD)
    # this will need to get changed later...
    offset = np.ones_like(x) * pwcpos
    offset_scaled = pwcpos_offset_scaler(offset)

    # polynomial features
    poly_features = get_poly_features(x_scaled, offset_scaled)
    wavelengths = coef @ poly_features + intercept

    return wavelengths
