"""
Unit tests for straylight step configuration
"""

from jwst.datamodels import CubeModel
from jwst.straylight import StraylightStep
import numpy as np
import pytest


@pytest.fixture(scope='module')
def miri_mrs_short_tso():
    """Set up MIRI MRS SHORT TSO data"""

    image = CubeModel((5, 1024, 1032))
    image.data = np.random.random((5, 1024, 1032))
    image.meta.instrument.name = 'MIRI'
    image.meta.instrument.detector = 'MIRIFUSHORT'
    image.meta.exposure.type = 'MIR_MRS'
    image.meta.instrument.channel = '12'
    image.meta.instrument.band = 'SHORT'
    image.meta.filename = 'test_miri_short_tso.fits'
    image.meta.observation.date = '2019-01-01'
    image.meta.observation.time = '10:10:10'
    return image


def test_call_straylight_mrsshort_tso(_jail, miri_mrs_short_tso):
    """Test step is skipped for MRS IFUSHORT TSO data"""
    result = StraylightStep.call(miri_mrs_short_tso)
    assert result.meta.cal_step.straylight == 'SKIPPED'
