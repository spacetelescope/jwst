import numpy as np
from .. import util


def test_gentle_asarray():
    x = np.array([('abc', 1.0)], dtype=[
        ('FOO', 'S3'),
        ('BAR', '>f8')])

    new_dtype = [('foo', '|S3'), ('bar', '<f8')]

    y = util.gentle_asarray(x, new_dtype)

    assert y['BAR'][0] == 1.0
