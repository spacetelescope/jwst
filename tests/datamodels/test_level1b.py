"""Test Level1bModel"""

import pytest
import numpy as np

from jwst.datamodels import Level1bModel


@pytest.mark.xfail
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
