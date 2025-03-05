"""Fixtures for AMI tests."""

import pytest
import numpy as np
from jwst.stpipe import Step
import stdatamodels.jwst.datamodels as dm
from jwst.ami.bp_fix import filthp_d


PXSC_DEG = 65.6 / (60.0 * 60.0 * 1000)
PXSC_RAD = PXSC_DEG * np.pi / (180)
PXSC_MAS = PXSC_DEG * 3600 * 1000


@pytest.fixture(scope="package")
def example_model():
    """
    Create a simple CubeModel simulating input to the ami3 pipeline.

    Returns
    -------
    CubeModel
        A simple CubeModel simulating NIRISS AMI data.
    """
    model = dm.CubeModel((2, 81, 81))
    # some non-zero data is required as this step will center
    # the image and find the centroid (both fail with all zeros)
    model.data[:, 24, 24] = 1
    model.data[:, 28, 28] = 1
    model.meta.instrument.name = "NIRISS"
    model.meta.instrument.filter = "F277W"
    model.meta.subarray.name = "SUB80"
    model.meta.observation.date = "2021-12-26"
    model.meta.observation.time = "00:00:00"
    model.meta.target.proposer_name = ""
    model.meta.program.pi_name = "someone"
    model.meta.target.catalog_name = ""
    model.meta.visit.start_time = "2022-06-05 12:15:41.5020000"
    model.meta.wcsinfo.roll_ref = 171.8779402866089
    model.meta.wcsinfo.v3yangle = 0.56126717
    model.meta.filename = "test_calints.fits"
    model.meta.instrument.pupil = "NRM"
    model.meta.exposure.type = "NIS_AMI"
    return model


@pytest.fixture(scope="package")
def circular_pupil():
    """
    Make a simple circular pupil mask.

    Returns
    -------
    np.ndarray
        A 1024x1024 array with a circular pupil of radius 0.2x the array shape.
    """
    shape = (1024, 1024)
    r = 0.2
    x = np.linspace(-1, 1, shape[0])
    y = np.linspace(-1, 1, shape[1])
    xx, yy = np.meshgrid(x, y)
    rr = np.sqrt(xx**2 + yy**2)
    pupil = np.zeros(shape)
    pupil[rr < r] = 1
    return pupil


@pytest.fixture(scope="package")
def nrm_model_circular(circular_pupil):
    """
    Make a simple NRMModel with a circular pupil.

    Returns
    -------
    NRMModel
        A simple NRMModel with a circular pupil.
    """
    return dm.NRMModel(nrm=circular_pupil)


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
def bandpass(example_model):
    """
    Simulate the bandpass of the example_model with a top-hat function.

    Returns
    -------
    np.ndarray
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
