"""Test various utility functions"""
import pytest

import numpy as np

from jwst.datamodels import SlitModel
from jwst.resample.resample_spec import find_dispersion_axis
from jwst.resample.resample_utils import build_mask


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
    assert np.array_equal(result, expected)


def test_find_dispersion_axis():
    """
    Test the find_dispersion_axis() function
    """
    dm = SlitModel()

    dm.meta.wcsinfo.dispersion_direction = 1    # horizontal
    assert find_dispersion_axis(dm) == 0        # X axis for wcs functions

    dm.meta.wcsinfo.dispersion_direction = 2    # vertical
    assert find_dispersion_axis(dm) == 1        # Y axis for wcs functions
