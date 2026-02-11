"""Fixtures for AMI tests."""

import numpy as np
import pytest
import stdatamodels.jwst.datamodels as dm

from jwst.ami.bp_fix import calc_pupil_support, filthp_d
from jwst.ami.tests import helpers
from jwst.stpipe import Step


@pytest.fixture(scope="package")
def example_model():
    """
    Create a simple CubeModel simulating input to the ami3 pipeline.

    Returns
    -------
    CubeModel
        A simple CubeModel simulating NIRISS AMI data.
    """
    return helpers.example_model()


@pytest.fixture(scope="package")
def circular_pupil():
    """
    Make a simple circular pupil mask.

    Returns
    -------
    ndarray
        A 1024x1024 array with a circular pupil of radius 0.2x the array shape.
    """
    return helpers.circular_pupil()


@pytest.fixture(scope="package")
def nrm_model(example_model):
    """
    Retrieve a real NRMModel reference file from CRDS.

    Returns
    -------
    NRMModel
        The NRMModel reference file for the filter defined by example_model
    """
    nrm_reffile = Step().get_reference_file(example_model, "nrm")
    return dm.NRMModel(nrm_reffile)


@pytest.fixture(scope="package")
def nrm_model_circular(circular_pupil, nrm_model):
    """
    Make a simple NRMModel with a circular pupil.

    Returns
    -------
    NRMModel
        A simple NRMModel with a circular pupil.
    """
    model = nrm_model.copy()
    model.nrm = circular_pupil
    return model


@pytest.fixture(scope="package")
def nrm_psf(example_model, nrm_model):
    """
    Compute the PSF from the nrm model.

    Returns
    -------
    psf : ndarray
        The PSF computed from the nrm model.
    """
    filt = example_model.meta.instrument.filter
    fov_npix = example_model.data.shape[1]
    return calc_pupil_support(filt, fov_npix, helpers.PXSC_RAD, nrm_model.nrm)


@pytest.fixture(scope="package")
def bandpass(example_model):
    """
    Simulate the bandpass of the example_model with a top-hat function.

    Returns
    -------
    ndarray
        An Nx2 array with the first column being the throughput and the second
        column being the wavelength in meters.
    """
    sz = 99
    fractional_pad = 0.25
    filt = example_model.meta.instrument.filter
    wlhp_l, wlhp_h = filthp_d[filt]

    # wlhp_l and wlhp_h are the cut-on and cut-off points.
    # outside of this range, throughput should be low. pad with 0.001 out to 0.5x
    # the range of the filter
    wl_low = wlhp_l - fractional_pad * (wlhp_h - wlhp_l)
    wl_high = wlhp_h + fractional_pad * (wlhp_h - wlhp_l)
    wls = np.linspace(wl_low, wl_high, sz).astype(np.float32)
    weights = np.ones_like(wls)
    weights[wls < wlhp_l] = 0.001
    weights[wls > wlhp_h] = 0.001
    return np.array(list(zip(weights, wls, strict=True)))


@pytest.fixture(scope="package")
def throughput_model(example_model, bandpass):
    """
    Mock a datamodels.ThroughputModel.

    Returns
    -------
    ThroughputModel
        A ThroughputModel with the same filter as the example_model and a top-hat
        bandpass defined by the bandpass parameter.
    """
    model = dm.ThroughputModel()
    model.meta.instrument.filter = example_model.meta.instrument.filter

    # build filter table from bandpass we already have
    # wavelengths are assumed to be angstrom by synphot, so convert them here
    wl_angstrom = bandpass[:, 1] * 1e10
    filter_table = np.rec.array(
        list(zip(wl_angstrom, bandpass[:, 0], strict=True)),
        names="wavelength,throughput",
        formats="f4,f4",
    )
    model.filter_table = filter_table
    return model
