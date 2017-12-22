# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Test jwst.transforms
"""
from numpy.testing.utils import assert_allclose
from ..import models


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
