import numpy as np
import pytest

from jwst.resample.resample_spec import find_dispersion_axis


def test_find_dispersion_axis():
    """
    Test the find_dispersion_axis() function
    """
    wavelengths = np.arange(100) * 0.1 + np.exp(0.1) * 13.0
    # [14.36722193, 14.46722193, 14.56722193, ... 24.16722193, 24.26722193]

    wavelengths_horizontal = np.tile(wavelengths, 15).reshape(15, 100)
    assert find_dispersion_axis(wavelengths_horizontal) == 0

    wavelengths_vertical = np.repeat(wavelengths, 15).reshape(100, 15)
    assert find_dispersion_axis(wavelengths_vertical) == 1

    # Make sure it works for decreasing wavelengths
    assert find_dispersion_axis(np.fliplr(wavelengths_horizontal)) == 0
    assert find_dispersion_axis(np.flipud(wavelengths_vertical)) == 1

    # Make sure it works if there are NaNs
    wavelengths_horizontal[:,0] = np.nan
    assert find_dispersion_axis(wavelengths_horizontal) == 0

    wavelengths_vertical[:,0] = np.nan
    assert find_dispersion_axis(wavelengths_vertical) == 1

    # Make sure if wavelengths don't change it produces an error
    wavelengths_zeros = np.zeros((15, 100))
    with pytest.raises(RuntimeError):
        find_dispersion_axis(wavelengths_zeros)
