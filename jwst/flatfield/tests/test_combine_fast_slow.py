"""
Test for flat_field.combine_fast_slow
"""
import numpy as np
from scipy.integrate import quad
from astropy.modeling import polynomial
from stdatamodels.jwst.datamodels import dqflags

from jwst.flatfield.flat_field import combine_fast_slow


def test_combine_fast_slow():
    """Assign an array of values to tab_wl; these are strictly increasing,
    but they're not uniformly spaced.
    The abscissas and weights in combine_fast_slow are for three-point
    Gaussian integration, which will (in principle) integrate a fifth-order
    polynomial exactly.  So for a test, we use a set of six polynomial
    coefficients to create the tab_flat array by evaluating the polynomial
    with those coefficients over the tab_wl array.
    """

    # Generate an array (tab_wl) of wavelengths, not uniformly spaced.
    wl_coeff = {'c0': 5., 'c1': 1., 'c2': -0.05}
    wl_poly = polynomial.Polynomial1D(degree=2, **wl_coeff)
    tab_index = np.arange(6000, dtype=np.float64) / 1000.
    tab_wl = wl_poly(tab_index)
    del tab_index

    coeff = {'c0': -41.9,
             'c1': 30.7,
             'c2': -8.3,
             'c3': 1.15,
             'c4': -7.8e-2,
             'c5': 2.0e-3}
    poly = polynomial.Polynomial1D(degree=5, **coeff)
    tab_flat = poly(tab_wl)

    wl0 = 7.37
    dwl0 = 2.99
    lower = wl0 - dwl0 / 2.
    upper = wl0 + dwl0 / 2.
    # Compute the integral of the polynomial over wavelength.
    # This function returns a tuple, the integral and an error estimate.
    integral = quad(poly, lower, upper)[0]
    # We want the average over the interval, not the integral itself.
    correct_value = integral / (upper - lower)

    # Make a simple flat with the expected bin width in wavelength
    wl_1d = np.arange(-5, 5) * dwl0 + wl0
    wl = np.tile(wl_1d, (10, 1))
    flat_2d = np.ones((10, 10), dtype=float)
    flat_dq = np.zeros((10, 10), dtype=np.uint32)
    dispaxis = 1

    value, new_dq = combine_fast_slow(wl, flat_2d, flat_dq, tab_wl,
                                      tab_flat, dispaxis)

    # Column 5 is the expected value
    assert np.allclose(value[:, 5], correct_value, rtol=1.e-8, atol=1.e-8)
    assert np.all(new_dq[:, 5] == 0)

    # Columns 0-2 are bad (negative wavelengths, not marked in DQ)
    assert np.all(value[:, :3] == 1)
    assert np.all(new_dq[:, :3] == 0)

    # Columns 3, 4, 6-9 are not covered by the tabular data
    # (missing values, marked in DQ)
    bad_value = dqflags.pixel['NO_FLAT_FIELD'] | dqflags.pixel['DO_NOT_USE']
    assert np.all(value[:, (3, 4, 6, 7, 8, 9)] == 1)
    assert np.all(new_dq[:, (3, 4, 6, 7, 8, 9)] == bad_value)
