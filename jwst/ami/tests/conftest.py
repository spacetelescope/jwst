"""Fixtures for AMI tests."""

from stdatamodels.jwst import datamodels
import pytest


@pytest.fixture()
def example_model():
    """Create a simple CubeModel simulating input to the ami3 pipeline."""
    model = datamodels.CubeModel((2, 81, 81))
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
