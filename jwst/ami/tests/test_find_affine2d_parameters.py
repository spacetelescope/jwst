"""Unit tests for AMI find_affine2d_parameters module."""

import numpy as np

from jwst.ami import find_affine2d_parameters, utils
from .conftest import PXSC_DEG


def test_create_afflist_rot():
    rotdegs = np.arange(1, 16)
    alist = find_affine2d_parameters.create_afflist_rot(rotdegs)
    assert len(alist) == len(rotdegs)
    for aff, rot in zip(alist, rotdegs):
        assert isinstance(aff, utils.Affine2d)
        assert np.isclose(aff.rotradccw, np.pi * rot / 180.0)

    # check that singularity has been avoided at 15 degrees
    assert alist[-1].rotradccw != np.pi * rotdegs[-1] / 180.0


def test_find_rotation(example_model, nrm_model, bandpass):
    imagedata = example_model.data[0]
    psf_offset = np.array([0.0, 0.0])
    rotdegs = np.arange(-3, 4)
    npix = imagedata.shape[0]
    over = 3
    holeshape = "hex"

    new_affine2d = find_affine2d_parameters.find_rotation(
        imagedata, nrm_model, psf_offset, rotdegs, PXSC_DEG, npix, bandpass, over, holeshape
    )

    assert isinstance(new_affine2d, utils.Affine2d)
    retrieved_rot_deg = new_affine2d.rotradccw * 180 / np.pi
    assert rotdegs[0] < retrieved_rot_deg
    assert retrieved_rot_deg < rotdegs[-1]
