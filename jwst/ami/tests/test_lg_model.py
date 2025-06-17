"""Unit tests for AMI lg_model module."""

import pytest
import numpy as np
from numpy.testing import assert_allclose

from jwst.ami import lg_model
from .conftest import PXSC_RAD


PSF_FOV = 21


@pytest.fixture(scope="module", params=[1, 3])  # parametrize oversample
def lgmodel(request, nrm_model, bandpass):
    """Create an LgModel object."""
    model = lg_model.LgModel(
        nrm_model,
        PXSC_RAD,
        bandpass,
        mask="jwst_ami",
        holeshape="hex",
        over=request.param,
        phi=None,
        chooseholes=None,
        affine2d=None,
    )
    return model


@pytest.mark.parametrize("psf_offset", [(0, 0), (1, -3)])
def test_simulate(request, lgmodel, bandpass, psf_offset):
    over = lgmodel.over
    psf = lgmodel.simulate(PSF_FOV, psf_offset=psf_offset)
    psf_over = lgmodel.psf_over

    # check attribute setting by this method
    assert np.all(psf == lgmodel.psf)
    assert np.all(bandpass == lgmodel.bandpass)

    # Check shapes and basic values
    assert psf.shape == (PSF_FOV, PSF_FOV)
    assert psf_over.shape == (PSF_FOV * over, PSF_FOV * over)
    assert np.all(psf >= 0)
    assert np.sum(np.isnan(psf)) == 0
    assert np.sum(np.isinf(psf)) == 0

    # Check that the PSF is centered, respecting offset
    expected_center = (int(PSF_FOV / 2 + psf_offset[1]), int(PSF_FOV / 2 + psf_offset[0]))
    max_idx = np.unravel_index(np.argmax(psf), psf.shape)
    assert np.all(max_idx == expected_center)


@pytest.mark.parametrize("psf_offset", [(0, 0), (1, -3)])
def test_make_model(lgmodel, psf_offset):
    fringemodel = lgmodel.make_model(
        PSF_FOV,
        psf_offset=psf_offset,
    )

    # Check attribute setting by this method
    assert np.all(fringemodel == lgmodel.model)
    assert lgmodel.fov == PSF_FOV
    assert hasattr(lgmodel, "model_beam")
    assert hasattr(lgmodel, "fringes")

    # Check shapes and basic values
    n_coeffs = lgmodel.N * (lgmodel.N - 1) + 2
    sz = PSF_FOV * lgmodel.over
    assert fringemodel.shape == (PSF_FOV, PSF_FOV, n_coeffs)
    assert lgmodel.model_beam.shape == (sz, sz)
    assert lgmodel.fringes.shape == (n_coeffs - 1, sz, sz)
    assert np.sum(np.isnan(fringemodel)) == 0
    assert np.sum(np.isinf(fringemodel)) == 0

    # Check that the beam is centered, respecting offset
    expected_center = (
        int(sz / 2 + psf_offset[1] * lgmodel.over),
        int(sz / 2 + psf_offset[0] * lgmodel.over),
    )
    max_idx = np.unravel_index(np.argmax(lgmodel.model_beam), lgmodel.model_beam.shape)
    assert np.all(max_idx == expected_center)


@pytest.mark.parametrize("weighted", [True, False])
def test_fit_image(example_model, lgmodel, weighted):
    """
    Test fit_image and create_modelpsf methods.

    These two methods are intended to be called in sequence.
    """
    image = example_model.data[0]
    fov = image.shape[0]

    fringemodel = lgmodel.make_model(
        fov,
    )
    lgmodel.fit_image(
        image,
        model_in=fringemodel,
        weighted=weighted,
    )
    lgmodel.create_modelpsf()

    # Check attribute setting by this method
    for attr in [
        "rawDC",
        "flux",
        "soln",
        "fringeamp",
        "fringephase",
        "fringepistons",
        "redundant_cps",
        "t3_amplitudes",
        "redundant_cas",
        "q4_phases",
        "modelpsf",
    ]:
        assert hasattr(lgmodel, attr)

    # Check shapes and basic values
    assert lgmodel.soln.size == 44
    assert lgmodel.fringeamp.size == 21
    assert lgmodel.fringephase.size == 21
    assert lgmodel.fringepistons.size == 7
    assert lgmodel.redundant_cps.size == 35
    assert lgmodel.t3_amplitudes.size == 35
    assert lgmodel.redundant_cas.size == 35
    assert lgmodel.q4_phases.size == 35
    assert lgmodel.modelpsf.shape == image.shape

    # check soln normalized
    assert lgmodel.soln[0] == 1.0
