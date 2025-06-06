"""Unit tests for AMI analyticnrm2 module."""

import pytest
import numpy as np
from numpy.testing import assert_allclose

from jwst.ami import analyticnrm2, utils


@pytest.fixture
def setup_sf():
    """
    Initialize values for these parameters needed for the analyticnrm2 tests.

    Returns
    -------
    pixel (optional, via **kwargs) : float
        pixel scale
    fov : integer
        number of detector pixels on a side
    oversample : integer
        oversampling factor
    ctrs : float, float
        coordinates of hole centers
    d : float
        hole diameter
    lam : float
        wavelength
    phi : float
        distance of fringe from hole center in units of waves
    centering : string
        if set to 'PIXELCENTERED' or unspecified, the offsets will be set to
        (0.5,0.5); if set to 'PIXELCORNER', the offsets will be set to
        (0.0,0.0).
    aff : Affine2d object
        Affine2d object
    """
    pixel = 3.1125038327232215e-07
    fov = 2
    oversample = 3
    ctrs = np.array(
        [
            [0.06864653, -2.63910736],
            [-2.28553695, -0.05944972],
            [2.31986022, -1.26010406],
            [-2.31986022, 1.26010406],
            [-1.19424838, 1.94960579],
            [2.25121368, 1.3790035],
            [1.09127858, 2.00905525],
        ]
    )
    d = 0.8
    lam = 2.3965000082171173e-06
    phi = np.zeros(7, dtype=np.float32)
    centering = (0.5, 0.5)
    aff_obj = utils.Affine2d(rotradccw=0.4)

    return pixel, fov, oversample, ctrs, d, lam, phi, centering, aff_obj


@pytest.mark.parametrize("holeshape", ["circ", "circonly", "hex", "hexonly", "fringeonly"])
def test_analyticnrm2_psf(setup_sf, holeshape):
    """Test of psf() in the analyticnrm2 module"""
    pixel, fov, oversample, ctrs, d, lam, phi, psf_offset, aff_obj = setup_sf

    computed_psf = analyticnrm2.psf(
        pixel, fov, oversample, ctrs, d, lam, phi, psf_offset, aff_obj, shape=holeshape
    )
    # basic checks of type and shape
    assert isinstance(computed_psf, np.ndarray)
    assert computed_psf.shape == (fov * oversample, fov * oversample)

    # ensure PSF offset is applied
    assert np.argmax(computed_psf) == 28


def test_analyticnrm2_asf_hex(setup_sf):
    """Test of asf_hex() in the analyticnrm2 module FOR HEX"""
    pixel, fov, oversample, _ctrs, d, lam, _phi, psf_offset, aff_obj = setup_sf

    asf = analyticnrm2.asf_hex(pixel, fov, oversample, d, lam, psf_offset, aff_obj)

    true_asf = np.array(
        [
            [
                0.82125698 + 7.84095011e-16j,
                0.83091456 + 2.48343013e-14j,
                0.83785899 - 2.49800181e-16j,
                0.84204421 - 1.80411242e-16j,
                0.8434424 - 2.91433544e-16j,
                0.84204424 - 1.24900090e-16j,
            ],
            [
                0.83091447 + 4.09394740e-16j,
                0.84064761 + 1.38777878e-15j,
                0.8476463 + 1.29063427e-15j,
                0.85186417 - 6.17561557e-16j,
                0.85327325 + 2.98372438e-16j,
                0.85186418 + 1.90125693e-15j,
            ],
            [
                0.83785894 + 1.68268177e-16j,
                0.84764629 + 1.07552856e-16j,
                0.8546839 + 6.38378239e-16j,
                0.8589252 - 1.65145675e-15j,
                0.8603421 - 9.29811783e-16j,
                0.8589252 + 1.15185639e-15j,
            ],
            [
                0.84204421 - 6.59194921e-17j,
                0.85186417 - 6.70470623e-16j,
                0.8589252 + 8.91214186e-16j,
                0.86318061 - 3.46944695e-16j,
                0.86460222 + 2.08166817e-17j,
                0.86318061 - 5.34294831e-16j,
            ],
            [
                0.84344243 + 2.28983499e-16j,
                0.85327326 + 2.98719383e-15j,
                0.8603421 + 5.02722863e-15j,
                0.86460222 + 5.48866508e-15j,
                0.8660254 + 0.00000000e00j,
                0.86460222 + 5.29611077e-15j,
            ],
            [
                0.84204425 - 1.48492330e-15j,
                0.85186418 + 6.03683770e-16j,
                0.8589252 + 5.68989300e-16j,
                0.86318061 + 2.77555756e-16j,
                0.86460222 - 1.72431514e-15j,
                0.86318061 - 5.54070678e-15j,
            ],
        ]
    )

    assert_allclose(asf, true_asf, atol=1e-7)


def test_analyticnrm2_interf(setup_sf):
    """Test of interf() in the analyticnrm2 module"""
    ASIZE = 4
    kx = np.arange(ASIZE * ASIZE).reshape((ASIZE, ASIZE))
    ky = np.arange(ASIZE * ASIZE).reshape((ASIZE, ASIZE))
    vv = np.arange(ASIZE)

    for ii in np.arange(ASIZE):
        kx[:, ii] = vv
        ky[ii, :] = vv

    pixel, _fov, oversample, ctrs, _d, lam, phi, _centering, aff_obj = setup_sf

    pitch = pixel / float(oversample)

    c = (ASIZE / 2.0, ASIZE / 2)

    interference = analyticnrm2.interf(
        kx, ky, ctrs=ctrs, phi=phi, lam=lam, pitch=pitch, c=c, affine2d=aff_obj
    )

    true_interference = np.array(
        [
            [
                2.6870043 + 1.24219632j,
                4.01721904 + 0.66189711j,
                4.2132531 + 0.21372447j,
                3.18675131 - 0.03818252j,
            ],
            [
                3.8517604 + 1.53442862j,
                5.71582424 + 0.84829672j,
                6.24380079 + 0.2201634j,
                5.25470657 - 0.31113349j,
            ],
            [
                4.02194801 + 1.32112798j,
                6.1888738 + 0.66733046j,
                7.0 + 0.0j,
                6.1888738 - 0.66733046j,
            ],
            [
                3.07194559 + 0.75829976j,
                5.25470657 + 0.31113349j,
                6.24380079 - 0.2201634j,
                5.71582424 - 0.84829672j,
            ],
        ]
    )

    assert_allclose(interference, true_interference, atol=1e-7)


def test_analyticnrm2_phasor():
    """Test of phasor() in the analyticnrm2 module"""
    ASIZE = 4
    kx = np.arange(ASIZE * ASIZE).reshape((ASIZE, ASIZE))

    for ii in np.arange(ASIZE):
        kx[:, ii] = ii

    ky = kx.transpose()

    hx = 0.06864653345335156
    hy = -2.6391073592116028

    lam = 2.3965000082171173e-06
    phi = 0.0
    pitch = 1.0375012775744072e-07

    aff_obj = utils.Affine2d(rotradccw=0.4)

    result = analyticnrm2.phasor(kx, ky, hx, hy, lam, phi, pitch, aff_obj)

    true_result = np.array(
        [
            [
                1.0 + 0.0j,
                0.96578202 + 0.25935515j,
                0.86546981 + 0.50096108j,
                0.70592834 + 0.70828326j,
            ],
            [
                0.78476644 + 0.61979161j,
                0.59716716 + 0.80211681j,
                0.36870018 + 0.92954837j,
                0.11500085 + 0.99336539j,
            ],
            [
                0.23171672 + 0.97278331j,
                -0.02850852 + 0.99959355j,
                -0.28678275 + 0.95799564j,
                -0.52543073 + 0.85083638j,
            ],
            [
                -0.42107943 + 0.90702377j,
                -0.64191223 + 0.76677812j,
                -0.81881514 + 0.57405728j,
                -0.93968165 + 0.34205027j,
            ],
        ]
    )

    assert_allclose(result, true_result, atol=1e-7)
