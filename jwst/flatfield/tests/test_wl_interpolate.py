"""
Test for flat_field.wl_interpolate
"""
import math

import numpy as np

from jwst.flatfield.flat_field import wl_interpolate


def test_wl_interpolate():

    tab_flat = np.array([2.99794324, 3.61644058, 4.12867112, 4.57353729,
                         4.95631171, 5.27168452, 5.51551932, 5.68937161,
                         5.80078189, 5.86148773], dtype=np.float64)

    tab_neg_wl = np.array([-1., -0.418, 0.128, 0.638, 1.112, 1.55, 1.952,
                           2.318, 2.648, 2.942], dtype=np.float64)
    # Check that this is not a valid wavelength, even though
    # wavelength > tab_neg_wl[0].
    wavelength = 0.
    value = wl_interpolate(wavelength, tab_neg_wl, tab_flat)
    assert value is None

    tab_wl = np.array([5., 5.582, 6.128, 6.638, 7.112, 7.55, 7.952,
                       8.318, 8.648, 8.942], dtype=np.float64)
    wavelength = 4.                     # out of bounds
    value = wl_interpolate(wavelength, tab_wl, tab_flat)
    assert value is None

    wavelength = 9.                     # out of bounds
    value = wl_interpolate(wavelength, tab_wl, tab_flat)
    assert value is None

    p = 0.25
    q = 1. - p
    wavelength = q * tab_wl[2] + p * tab_wl[3]
    expected_value = q * tab_flat[2] + p * tab_flat[3]
    value = wl_interpolate(wavelength, tab_wl, tab_flat)
    assert math.isclose(value, expected_value, rel_tol=1.e-12, abs_tol=1.e-12)
