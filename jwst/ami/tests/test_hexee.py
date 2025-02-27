"""Unit tests for AMI hexee module."""

import pytest
import numpy as np
from numpy.testing import assert_allclose

from jwst.ami import hexee


@pytest.fixture(scope="module")
def hexee_params():
    """
    Initialize values for parameters needed for the hexee tests.

    Returns
    -------
    xi : 2D float array
        hexagon's coordinate center at center of symmetry, along flat edge
    eta : 2D float array
        hexagon's coordinate center at center of symmetry, normal to xi;
        (not currently used)
    kwargs : dict
        c : tuple(float, float)
            coordinates of center
        pixel : float
            pixel scale
        d : float
            flat-to-flat distance across hexagon
        lam : float
            wavelength
        minus : bool
            if set, use flipped sign of xi in calculation
    """
    xdim, ydim = 3, 3
    xi = np.zeros(ydim * xdim).reshape((ydim, xdim))
    eta = np.zeros(ydim * xdim).reshape((ydim, xdim))

    for ii in range(ydim):
        xi[ii, :] = ii
        eta[:, ii] = ii

    kwargs = {
        "d": 0.8,
        "c": (28.0, 28.0),
        "lam": 2.3965000082171173e-06,
        "pixel": 1.0375012775744072e-07,
        "minus": False,
    }

    return xi, eta, kwargs


def test_hexee_g_eeag(hexee_params):
    """
    Test of g_eeag() in the hexee module.

    Calculate the Fourier transform of one half of a hexagon that is
    bisected from one corner to its diametrically opposite corner.
    """
    xi, eta, kwargs = hexee_params
    g = hexee.g_eeag(xi, eta, **kwargs)

    true_g = np.array(
        [
            [-0.04454286 + 0.05015766j, -0.04164985 + 0.06041733j, -0.03830953 + 0.07099764j],
            [-0.04072437 + 0.05375103j, -0.03729262 + 0.06415232j, -0.03340318 + 0.07486623j],
            [-0.03657856 + 0.05703437j, -0.03258885 + 0.06754246j, -0.02813134 + 0.07835476j],
        ]
    )

    assert_allclose(g, true_g)


def test_hexee_glimit(hexee_params):
    """
    Test of glimit() in the hexee module.

    Calculate the analytic limit of the Fourier transform of one half of the
    hexagon along eta=0.
    """
    xi, _, kwargs = hexee_params
    g = hexee.glimit(xi, **kwargs)

    true_g = np.array(
        [
            [0.07105571 + 0.28088478j, 0.07105571 + 0.28088478j, 0.07105571 + 0.28088478j],
            [0.08609692 + 0.28598645j, 0.08609692 + 0.28598645j, 0.08609692 + 0.28598645j],
            [0.10178022 + 0.29008864j, 0.10178022 + 0.29008864j, 0.10178022 + 0.29008864j],
        ]
    )

    assert_allclose(g, true_g)
