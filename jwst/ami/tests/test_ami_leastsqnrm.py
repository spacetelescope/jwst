import numpy as np

from jwst.ami import leastsqnrm


def test_weighted_operations():
    img = np.arange(4, dtype='f4').reshape((2, 2))
    model = np.arange(8, dtype='f4').reshape((2, 2, 2))
    dq = np.zeros((2, 2), dtype='int32')
    dq[1, 1] = 1
    x, res, cond, singvals = leastsqnrm.weighted_operations(img, model, dq)

    np.testing.assert_almost_equal(x, [0.19512195, 0.24390244])
    assert cond is None
    np.testing.assert_almost_equal(singvals, [4.52769249, 0.])
