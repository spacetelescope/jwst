"""Unit tests for AMI hextransformee module."""

import numpy as np

from jwst.ami import hextransformee
from jwst.ami.bp_fix import filtwl_d
from jwst.ami.utils import Affine2d
from .conftest import PXSC_RAD


def test_hextransform(example_model):
    """Test the hextransform function."""

    shape = example_model.data.shape[1:]
    # make the center a bit off of data array center for testing
    # need to copy ctr before passing into hextransform because it gets modified inside
    ctr = [shape[0]//2 + 4, shape[1]//2 -3]
    c = ctr.copy()

    lam = filtwl_d[example_model.meta.instrument.filter]
    affine2d = Affine2d(rotradccw=0.4)

    hex_complex = hextransformee.hextransform(s=shape, c=c, d=0.8, lam=lam, pitch=PXSC_RAD, affine2d=affine2d)
    
    # test shape and dtype
    assert hex_complex.dtype == np.complex128
    assert hex_complex.shape == shape

    # test proper normalization and centering
    hex_real = np.abs(hex_complex)
    assert np.isclose(np.max(hex_real), np.sqrt(3) / 2.0, atol=1.0e-8)
    assert np.argmax(hex_real) == ctr[0] * shape[1] + ctr[1]
    assert (hex_real <= np.sqrt(3) / 2.0).all()
