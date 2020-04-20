# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Test jwst.transforms
"""
import pytest
from numpy.testing import assert_allclose

from jwst.transforms import models


#_RANDOM_SEED = 0x1337

"""
def test_logical():

    with NumpyRNGContext(0x1337):
        compareto = np.random.randn(10)
    with NumpyRNGContext(0x1338):
        val = np.random.randn(10)
    with NumpyRNGContext(0x1339):
        x = np.random.randn(10)
    l = models.Logical('GT', .5, 10)

    res = l(x)
    y = x.copy()
    y[np.greater(x, .5)] = 10
    assert_allclose(res, npres)
    l = models.Logical('lt', compareto, val)
    cond = np.less(x, compareto)
    y = x.copy()
    y[cond] = val[cond]
    assert_allclose(res, npres)
"""

def test_ideal_to_v23_roundtrip():
    """
    Test roundtripping of the transforms.
    """
    v2i = models.V2V3ToIdeal(.4, 450, 730, 1)
    x, y = 450, 730
    assert_allclose(v2i.inverse(*v2i(x, y)), (x, y))


@pytest.mark.parametrize(('wavelength', 'n'),
                         [(1e-6, 1.43079543),
                          (2e-6,  1.42575377),
                          (5e-6, 1.40061966)
                          ])
def test_refraction_index(wavelength, n):
    """
    Tests the computation of the refraction index.
    True values are from the ESA pipeline.
    Reference values are from the PRISM reference file from CV3.
    """
    temp_sys = 37.06107795068881 # in K
    tref = 35 # in K
    pref = 0 # in atm
    pressure_sys = 0 # in atm
    kcoef =  [0.58339748, 0.46085267, 3.8915394]
    lcoef = [0.00252643, 0.010078333, 1200.556]
    tcoef = [-2.66e-05, 0.0, 0.0, 0.0, 0.0, 0.0]
    n_pipeline = models.Snell.compute_refraction_index(wavelength, temp_sys,
                                                       tref, pref, pressure_sys,
                                                       kcoef, lcoef, tcoef)
    assert_allclose(n_pipeline, n)
