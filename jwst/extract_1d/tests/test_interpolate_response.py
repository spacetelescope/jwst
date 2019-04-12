"""
Test for extract_1d.interpolate_response
"""
import math

import numpy as np

from jwst.datamodels import SlitModel
from jwst.extract_1d import extract

def test_interpolate_response():

    wavelength = np.array([0.5, 5., 7., 18., 26.], dtype=np.float32)
    wl1 = np.arange(1, 21, dtype=np.float64)
    wl2 = np.arange(1, 21, dtype=np.float64)[::-1]      # decreasing order
    response = np.arange(20, dtype=np.float64) * 2. + 5.

    slit = SlitModel()          # this is only used for the dtype of relsens
    tbl1 = np.array(list(zip(wl1, response)), dtype=slit.relsens.dtype)
    tbl2 = np.array(list(zip(wl2, response)), dtype=slit.relsens.dtype)

    reciprocal_r1 = extract.interpolate_response(wavelength, tbl1,
                                                 verbose=False)
    assert len(reciprocal_r1) == len(wavelength)
    assert math.isclose(1. / reciprocal_r1[1], 13.0, rel_tol=1.e-6)
    assert reciprocal_r1[0] == 0.       # because wavelength is out of range
    assert reciprocal_r1[4] == 0.

    reciprocal_r2 = extract.interpolate_response(wavelength, tbl2,
                                                 verbose=False)
    assert len(reciprocal_r2) == len(wavelength)
    assert math.isclose(1. / reciprocal_r2[1], 35.0, rel_tol=1.e-6)
    assert reciprocal_r2[0] == 0.       # because wavelength is out of range
    assert reciprocal_r2[4] == 0.
