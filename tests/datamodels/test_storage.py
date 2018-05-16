import numpy as np
from .. import util


def test_gentle_asarray():
    x = np.array([('abc', 1.0)], dtype=[
        (str('FOO'), str('S3')),
        (str('BAR'), str('>f8'))])

    new_dtype = [(str('foo'), str('|S3')), (str('bar'), str('<f8'))]

    y = util.gentle_asarray(x, new_dtype)

    assert y['BAR'][0] == 1.0
