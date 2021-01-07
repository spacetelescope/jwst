"""
Test for flat_field.g_average
"""
import math

import numpy as np
from scipy.integrate import quad
from astropy.modeling import polynomial

from jwst.flatfield.flat_field import g_average


def test_g_average():

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

    # These three lines were copied from combine_fast_slow in flat_field.py
    d = math.sqrt(0.6) / 2.
    dx = np.array([-d, 0., d])
    wgt = np.array([5., 8., 5.]) / 18.

    value = g_average(wl0, dwl0, tab_wl, tab_flat, dx, wgt)
    assert math.isclose(value, correct_value, rel_tol=1.e-8, abs_tol=1.e-8)
