"""
Unit tests for gain_scale correction
"""

import numpy as np
import pytest
from stdatamodels.jwst.datamodels import CubeModel

from jwst.gain_scale import GainScaleStep
from jwst.gain_scale.gain_scale import do_correction


def test_correction(make_cubemodel):
    """Make sure correct gain factor is applied"""
    datmod = make_cubemodel(2, 50, 50)
    gf = datmod.meta.exposure.gain_factor
    # run on a copy (the detection to make the copy or not would have happened at _step.py)
    output = do_correction(datmod.copy(), gain_factor=gf)

    assert output.meta.cal_step.gain_scale == "COMPLETE"
    assert np.all(output.data == datmod.data * gf)
    assert np.all(output.err == datmod.err * gf)
    assert np.all(output.var_poisson == datmod.var_poisson * gf * gf)
    assert np.all(output.var_rnoise == datmod.var_rnoise * gf * gf)


def test_step(make_cubemodel):
    """Make sure correct gain factor is applied at step level"""
    datmod = make_cubemodel(2, 50, 50)
    gf = datmod.meta.exposure.gain_factor
    output = GainScaleStep.call(datmod)

    assert output.meta.cal_step.gain_scale == "COMPLETE"
    assert np.all(output.data == datmod.data * gf)
    assert np.all(output.err == datmod.err * gf)
    assert np.all(output.var_poisson == datmod.var_poisson * gf * gf)
    assert np.all(output.var_rnoise == datmod.var_rnoise * gf * gf)

    # Input is not modified
    assert output is not datmod
    assert datmod.meta.cal_step.gain_scale is None


def test_step_from_reference(make_cubemodel):
    """Test application from reference model."""
    datmod = make_cubemodel(2, 50, 50)
    datmod.meta.exposure.gain_factor = None

    output = GainScaleStep.call(datmod)

    assert output.meta.cal_step.gain_scale == "COMPLETE"
    assert not np.all(output.data == datmod.data)

    # Input is not modified
    assert output is not datmod
    assert datmod.meta.cal_step.gain_scale is None


def test_step_reference_missing(make_cubemodel):
    """Test skip for missing reference model."""
    datmod = make_cubemodel(2, 50, 50)
    datmod.meta.exposure.gain_factor = None

    output = GainScaleStep.call(datmod, override_gain="N/A")

    assert output.meta.cal_step.gain_scale == "SKIPPED"
    assert np.all(output.data == datmod.data)

    # Input is not modified
    assert output is not datmod
    assert datmod.meta.cal_step.gain_scale is None


@pytest.fixture(scope="function")
def make_cubemodel():
    """Cube model for testing"""

    def _dm(nints, ysize, xsize):
        # create the data arrays
        csize = (nints, ysize, xsize)
        data = np.random.randint(low=1, high=50, size=csize)
        err = np.random.randint(low=0.1, high=5.0, size=csize)
        var_p = err * err
        var_r = np.ones(csize) * 12

        # create a JWST datamodel
        dm = CubeModel(data=data, err=err, var_poisson=var_p, var_rnoise=var_r)

        dm.meta.instrument.name = "NIRSPEC"
        dm.meta.date = "2018-01-01"
        dm.meta.instrument.detector = "NRS1"
        dm.meta.observation.date = "2018-01-01"
        dm.meta.observation.time = "00:00:00"

        dm.meta.exposure.gain_factor = 2

        return dm

    return _dm
