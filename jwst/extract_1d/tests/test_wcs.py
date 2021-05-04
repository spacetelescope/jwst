import numpy as np
from numpy.testing import assert_allclose

from jwst.extract_1d import spec_wcs

import pytest


def test_spec_wcs_monotonic():
    ra, dec = 5.6, -72.3
    wave = np.array([1.2, 3.4, 4.5, 5.6, 6.9])
    wcsobj = spec_wcs.create_spectral_wcs(ra, dec, wave)

    inp_pix = [1.2, 3.5, 0, 1]
    res = (np.array([ra] * 4), np.array([dec] * 4), np.array([3.62, 6.25, 1.2, 3.4]))

    assert_allclose(wcsobj(inp_pix), res, atol=10**-14)
    assert_allclose(wcsobj.invert(*res), inp_pix, atol=10**-14)


def test_spec_wcs_nonmonotonic():
    wave = np.array([1.2, 3.5, 0.1, 5.6, 3.4])
    wcsobj = spec_wcs.create_spectral_wcs(1.2, 23, wave)
    with pytest.raises((NotImplementedError, ValueError)):
        wcsobj.invert(1.8)
