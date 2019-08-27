"""
Test for extract_1d._fit_background_model
"""
import math

import numpy as np

from jwst.extract_1d import extract1d

def test_fit_background_model():

    shape = (9, 5)
    image = np.arange(shape[0] * shape[1], dtype=np.float32).reshape(shape)
    x = 2
    j = 2

    b_lower = np.zeros(shape[1], dtype=np.float64) + 3.5    # 4, inclusive
    b_upper = np.zeros(shape[1], dtype=np.float64) + 4.5    # 4, inclusive
    bkglim = [[b_lower, b_upper]]
    bkg_order = 0

    (bkg_model, npts) = extract1d._fit_background_model(
                        image, x, j, bkglim, bkg_order)

    assert math.isclose(bkg_model(0.), 22.0, rel_tol=1.e-8, abs_tol=1.e-8)
    assert math.isclose(bkg_model(8.), 22.0, rel_tol=1.e-8, abs_tol=1.e-8)

    assert npts == 2

    image[:, 2] = np.nan

    (bkg_model, npts) = extract1d._fit_background_model(
                        image, x, j, bkglim, bkg_order)

    assert math.isclose(bkg_model(0.), 0.0, rel_tol=1.e-8, abs_tol=1.e-8)
    assert math.isclose(bkg_model(8.), 0.0, rel_tol=1.e-8, abs_tol=1.e-8)

    assert npts == 0
