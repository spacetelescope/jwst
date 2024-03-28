
import numpy as np
import pytest

import jwst.ami.utils as utils


@pytest.mark.parametrize("shape, center", [
    ((10, 10), (4.5, 4.5)),
    ((11, 11), (5, 5)),
])
def test_centerpoint(shape, center):
    assert utils.centerpoint(shape) == center


def test_find_centroid():
    arr = np.zeros((30, 30), dtype='f4')
    arr[15, 15] = 1
    thresh = 0.02
    assert np.allclose(utils.find_centroid(arr, thresh), (0.5, 0.5))


@pytest.mark.parametrize("mas, rad", [
    (206264.8062471, 0.001),
    (103132403.12355, 0.5),
])
def test_mas2rad(mas, rad):
    assert np.isclose(utils.mas2rad(mas), rad)
