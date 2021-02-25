"""Test various utility functions"""
import pytest

import numpy as np

from jwst.datamodels import SlitModel, ImageModel, dqflags
from jwst.resample.resample_spec import find_dispersion_axis
from jwst.resample.resample_utils import build_mask, build_driz_weight


DO_NOT_USE = dqflags.pixel["DO_NOT_USE"]
GOOD = dqflags.pixel["GOOD"]


DQ = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8])
BITVALUES = 2**0 + 2**2
BITVALUES_STR = f'{2**0}, {2**2}'
BITVALUES_INV_STR = f'~{2**0}, {2**2}'
JWST_NAMES = 'DO_NOT_USE, JUMP_DET'
JWST_NAMES_INV = '~' + JWST_NAMES
@pytest.mark.parametrize(
    'dq, bitvalues, expected', [
        (DQ, 0,                 np.array([1, 0, 0, 0, 0, 0, 0, 0, 0])),
        (DQ, BITVALUES,         np.array([1, 1, 0, 0, 1, 1, 0, 0, 0])),
        (DQ, BITVALUES_STR,     np.array([1, 1, 0, 0, 1, 1, 0, 0, 0])),
        (DQ, BITVALUES_INV_STR, np.array([1, 0, 1, 0, 0, 0, 0, 0, 1])),
        (DQ, JWST_NAMES,        np.array([1, 1, 0, 0, 1, 1, 0, 0, 0])),
        (DQ, JWST_NAMES_INV,    np.array([1, 0, 1, 0, 0, 0, 0, 0, 1])),
        (DQ, None,              np.array([1, 1, 1, 1, 1, 1, 1, 1, 1])),
    ]
)
def test_build_mask(dq, bitvalues, expected):
    """Test logic of mask building

    Parameters
    ----------
    dq: numpy.array
        The input data quality array

    bitvalues: int or str
        The bitvalues to match against

    expected: numpy.array
        Expected mask array
    """
    result = build_mask(dq, bitvalues)
    np.testing.assert_array_equal(result, expected)


@pytest.mark.parametrize("weight_type, expected_value", [
    ("ivm", 0.1),
    ("exptime", 10),
    (None, 1.0),
    ]
)
def test_build_driz_weight(weight_type, expected_value):
    """Check that correct weight map is returned of different weight types"""
    model = ImageModel((10, 10))
    model.dq[0] = DO_NOT_USE
    model.meta.exposure.exposure_time = 10.0
    model.var_rnoise += 10.0

    weight_map = build_driz_weight(model, weight_type=weight_type, good_bits="GOOD")
    np.testing.assert_array_equal(weight_map[0], 0)
    np.testing.assert_array_equal(weight_map[1:], expected_value)
    assert weight_map.dtype == np.float32


@pytest.mark.parametrize("weight_type", ["ivm", "exptime"])
def test_build_driz_weight_zeros(weight_type):
    """Check that zero or not finite weight maps raise an error"""
    model = ImageModel((10, 10))
    model.var_rnoise += 0
    model.meta.exposure.exposure_time = 0.0

    with pytest.raises(ValueError):
        build_driz_weight(model, weight_type=weight_type)


def test_find_dispersion_axis():
    """
    Test the find_dispersion_axis() function
    """
    dm = SlitModel()

    dm.meta.wcsinfo.dispersion_direction = 1    # horizontal
    assert find_dispersion_axis(dm) == 0        # X axis for wcs functions

    dm.meta.wcsinfo.dispersion_direction = 2    # vertical
    assert find_dispersion_axis(dm) == 1        # Y axis for wcs functions
