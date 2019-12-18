"""
Unit tests for straylight step configuration
"""

from jwst.datamodels import IFUImageModel
from jwst.straylight import StraylightStep
import numpy as np
import pytest


@pytest.fixture(scope='module')
def miri_mrs_long():
    """Set up MIRI MRS Long data"""

    image = IFUImageModel((30, 30))
    image.data = np.random.random((30, 30))
    image.meta.instrument.name = 'MIRI'
    image.meta.instrument.detector = 'MIRIFULONG'
    image.meta.exposure.type = 'MIR_MRS'
    image.meta.instrument.channel = '34'
    image.meta.instrument.band = 'LONG'
    image.meta.filename = 'test_miri_long.fits'
    image.meta.observation.date = '2019-01-01'
    image.meta.observation.time = '10:10:10'
    return image


def test_call_straylight_mrslong(_jail, miri_mrs_long):
    """Test step is skipped for MRS IFULONG data"""
    result = StraylightStep.call(miri_mrs_long)
    assert result.meta.cal_step.straylight == 'SKIPPED'
