"""
Test for extract_1d.interpolate_response
"""
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

    r1 = extract.interpolate_response(wavelength, tbl1, verbose=False)
    assert len(r1) == len(wavelength)
    assert r1[1] == 13.0

    r2 = extract.interpolate_response(wavelength, tbl2, verbose=False)
    assert len(r2) == len(wavelength)
    assert r2[1] == 35.0
