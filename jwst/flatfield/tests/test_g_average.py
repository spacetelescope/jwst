"""
Test for flat_field.g_average
"""
import math

import numpy as np

from jwst.flatfield.flat_field import g_average

def evaluate_polynomial(x, coeff):
    """Evaluate a polynomial."""
    n = len(coeff)
    sum = np.zeros_like(x) + coeff[0]
    for i in range(1, n):
        sum = sum * x + coeff[i]
    return sum


def integrate_polynomial(lower, upper, coeff):
    """Integrate a polynomial."""

    # Compute the indefinite integral by modifying the coefficients,
    # setting the constant term to 0.
    n = len(coeff)
    icoeff = np.zeros(n + 1, dtype=np.float64)
    p = float(n - 1)
    # Don't explicitly assign a value to icoeff[-1], leave it at 0.
    for i in range(n):
        icoeff[i] = coeff[i] / (p + 1.)
        p -= 1.

    # Compute the definite integral.
    value = (evaluate_polynomial(upper, icoeff) -
             evaluate_polynomial(lower, icoeff))

    return value


def test_g_average():

    """Assign an array of values to tab_wl; these are strictly increasing,
    but they're not uniformly spaced.
    The abscissas and weights in combine_fast_slow are for three-point
    Gaussian integration, which will (in principle) integrate a fifth-order
    polynomial exactly.  So for a test, we use a set of six polynomial
    coefficients to create the tab_flat array by evaluating the polynomial
    with those coefficients over the tab_wl array.
    """

    # Fifth-order polynomial, created using x = 5 to 11 inclusive.
    coeff = np.array([2.0e-03, -7.8e-02, 1.15, -8.3, 30.7, -41.9],
                     dtype=np.float64)

    tab_index = np.arange(6000, dtype=np.float64) / 1000.
    wl_coeff = np.array([-0.05, 1., 5.], dtype=np.float64)
    tab_wl = evaluate_polynomial(tab_index, wl_coeff)
    tab_flat = evaluate_polynomial(tab_wl, coeff)

    wl0 = 7.37
    dwl0 = 2.99
    lower = wl0 - dwl0 / 2.
    upper = wl0 + dwl0 / 2.
    # This is the actual integral of the polynomial over wavelength, not
    # over the non-uniformly spaced tab_wl.
    integral = integrate_polynomial(lower, upper, coeff)
    # We want the average over the interval, not the integral itself.
    correct_value = integral / (upper - lower)

    # These three lines were copied from combine_fast_slow in flat_field.py
    d = math.sqrt(0.6) / 2.
    dx = np.array([-d, 0., d])
    wgt = np.array([5., 8., 5.]) / 18.

    value = g_average(wl0, dwl0, tab_wl, tab_flat, dx, wgt)
    assert math.isclose(value, correct_value, rel_tol=1.e-8, abs_tol=1.e-8)
