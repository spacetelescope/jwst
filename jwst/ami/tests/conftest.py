"""Fixtures for AMI tests."""

import pytest
import numpy as np
from jwst.stpipe import Step
import stdatamodels.jwst.datamodels as dm


@pytest.fixture()
def example_model():
    """Create a simple CubeModel simulating input to the ami3 pipeline."""
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


@pytest.fixture
def circular_pupil():
    """Make a simple circular pupil mask."""
    shape = (1024, 1024)
    r = 0.2
    x = np.linspace(-1, 1, shape[0])
    y = np.linspace(-1, 1, shape[1])
    xx, yy = np.meshgrid(x, y)
    rr = np.sqrt(xx**2 + yy**2)
    pupil = np.zeros(shape)
    pupil[rr < r] = 1
    return pupil


@pytest.fixture
def nrm_model_circular(circular_pupil):
    """Make a simple NRMModel with a circular pupil."""
    return dm.NRMModel(nrm=circular_pupil)


PXSC_DEG = 65.6 / (60.0 * 60.0 * 1000)
PXSC_RAD = PXSC_DEG * np.pi / (180)
PXSC_MAS = PXSC_DEG * 3600 * 1000


@pytest.fixture
def nrm_model(example_model):
    """Retrieve a real NRMModel reference file from CRDS."""
    nrm_reffile = Step().get_reference_file(example_model, "nrm")
    nrm_model = dm.NRMModel(nrm_reffile)
    return nrm_model
