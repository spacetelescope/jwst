import numpy as np

from jwst.outlier_detection.tso import moving_median_over_zeroth_axis


def test_rolling_median():
    time_axis = np.arange(10)
    expected_time_axis = np.array([1, 1, 2, 3, 4, 5, 6, 7, 8, 8])
    spatial_axis = np.ones((5, 5))
    arr = time_axis[:, np.newaxis, np.newaxis] * spatial_axis[np.newaxis, :, :]

    w = 3
    result = moving_median_over_zeroth_axis(arr, w)
    expected = expected_time_axis[:, np.newaxis, np.newaxis] * spatial_axis[np.newaxis, :, :]
    assert np.allclose(result, expected)
