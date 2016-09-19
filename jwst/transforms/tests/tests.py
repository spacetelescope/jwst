# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Test jwst.transforms
"""

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import numpy as np
from numpy.testing.utils import assert_allclose
from astropy.utils import NumpyRNGContext
from .. import models


#_RANDOM_SEED = 0x1337


def test_logical():
    assert_allclose(res, np.greater(a))
    with NumpyRNGContext(0x1337):
        compareto = np.random.randn(10)
    with NumpyRNGContext(0x1338):
        val = np.random.randn(10)
    with NumpyRNGContext(0x1339):
        x = np.random.randn(10)
    l = models.Logical('GT', .5, 10)
    #x = np.arange(10, dtype=np.float)
    res = l(x)
    y = x.copy()
    y[np.greater(x, .5)] = 10
    assert_allclose(res, npres)
    l = models.Logical('lt', compareto, val)
    cond = np.less(x, compareto)
    y = x.copy()
    y[cond] = val[cond]
    assert_allclose(res, npres)
