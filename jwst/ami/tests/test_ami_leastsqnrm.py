import numpy as np

from jwst.ami import leastsqnrm


def test_weighted_operations():
    img = np.arange(1, 5, dtype='f4').reshape((2, 2))
    model = np.arange(8, dtype='f4').reshape((2, 2, 2))
    dq = np.zeros((2, 2), dtype='int32')
    dq[1, 1] = 1
    x, res, cond, singvals = leastsqnrm.weighted_operations(img, model, dq)

    np.testing.assert_almost_equal(x, [-0.5, 1])
    assert cond is None
    np.testing.assert_almost_equal(singvals, [4.4870429, 0.1819676])
