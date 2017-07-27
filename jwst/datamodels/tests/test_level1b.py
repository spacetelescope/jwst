"""Test Level1bModel"""

import numpy as np

from .. import Level1bModel


def test_no_zeroframe():
    """Test for default zeroframe"""
    nx = 10
    ny = 10
    ngroups = 5
    nints = 2

    data = np.zeros((nints, ngroups, ny, nx), np.int16)
    model = Level1bModel(data)
    assert model.data.shape == (nints, ngroups, ny, nx)
    assert model.zeroframe.shape == (nints, ny, nx)
