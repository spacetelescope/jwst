import numpy as np
import pytest

from jwst.outlier_detection.tso import moving_median_over_zeroth_axis


@pytest.mark.parametrize(
    "w, expected_time_axis",
    [
        (3, [1, 1, 2, 3, 4, 5, 6, 7, 8, 8]),
        (4, [1.5, 1.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 7.5]),
    ],
)
def test_rolling_median(w, expected_time_axis):
    time_axis = np.arange(10)
    spatial_axis = np.ones((5, 5))
    arr = time_axis[:, np.newaxis, np.newaxis] * spatial_axis[np.newaxis, :, :]

    result = moving_median_over_zeroth_axis(arr, w)
    expected = (
        np.array(expected_time_axis)[:, np.newaxis, np.newaxis] * spatial_axis[np.newaxis, :, :]
    )
    assert np.allclose(result, expected)
