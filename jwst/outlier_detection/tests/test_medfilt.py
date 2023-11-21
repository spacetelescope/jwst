import warnings

import pytest
import numpy as np
import scipy.signal

from jwst.outlier_detection.outlier_detection_ifu import medfilt



@pytest.mark.parametrize("shape,kern_size", [
    ([7, 7], [3, 3]),
    ([7, 7], [3, 1]),
    ([7, 7], [1, 3]),
    ([7, 5], [3, 3]),
    ([5, 7], [3, 3]),
    ([42, 42], [7, 7]),
    ([42, 42], [7, 5]),
    ([42, 42], [5, 7]),
    ([42, 7, 5], [3, 3, 3]),
    ([5, 7, 42], [5, 5, 5]),
])
def test_medfilt_against_scipy(shape, kern_size):
    arr = np.arange(np.prod(shape), dtype='uint32').reshape(shape)
    result = medfilt(arr, kern_size)
    expected = scipy.signal.medfilt(arr, kern_size)
    np.testing.assert_allclose(result, expected)


@pytest.mark.parametrize("arr,kern_size,expected", [
    ([2, np.nan, 0], [3], [1, 1, 0]),
    ([np.nan, np.nan, np.nan], [3], [0, np.nan, 0]),
])
def test_medfilt_nan(arr, kern_size, expected):
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="All-NaN slice",
            category=RuntimeWarning
        )
        result = medfilt(arr, kern_size)
    np.testing.assert_allclose(result, expected)
