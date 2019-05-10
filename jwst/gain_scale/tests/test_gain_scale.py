"""
Unit tests for gain_scale correction
"""

from jwst.datamodels import RampModel
from jwst.gain_scale.gain_scale import do_correction
from jwst.gain_scale import GainScaleStep
import numpy as np
import pytest

def test_correction(make_rampmodel):
    """Make sure correct gain factor is applied
    """
    datmod = make_rampmodel(2, 50, 50)
    output = do_correction(datmod, gain_factor=datmod.meta.exposure.gain_factor)

    assert(output.meta.cal_step.gain_scale == 'COMPLETE')
    assert np.all(output.err == datmod.err * datmod.meta.exposure.gain_factor)
    assert np.all(output.data == datmod.data * datmod.meta.exposure.gain_factor)


def test_step(make_rampmodel):
    """Make sure correct gain factor is applied at step level
    """
    datmod = make_rampmodel(2, 50, 50)
    output = GainScaleStep.call(datmod)

    assert(output.meta.cal_step.gain_scale == 'COMPLETE')
    assert np.all(output.err == datmod.err * datmod.meta.exposure.gain_factor)
    assert np.all(output.data == datmod.data * datmod.meta.exposure.gain_factor)


@pytest.fixture(scope='function')
def make_rampmodel():
    '''Ramp model for testing'''

    # NRS1 and NRS2 are size  2048x2048 pixels
    def _dm(ngroups, ysize, xsize):
        # create the data and groupdq arrays
        nints = 2
        csize = (nints, ngroups, ysize, xsize)
        data = np.random.randint(low=1, high=50, size=csize)
        err = np.random.randint(low=0.1, high=5.0, size=csize)

        # create a JWST datamodel for NIRSPEC data
        dm = RampModel(data=data, err=err)

        dm.meta.instrument.name = 'NIRSPEC'
        dm.meta.date = '2018-01-01'
        dm.meta.instrument.detector= 'NRS1'
        dm.meta.observation.date = '2018-01-01'

        dm.meta.exposure.gain_factor = 2

        return dm

    return _dm
