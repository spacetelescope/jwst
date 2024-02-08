import stdatamodels.jwst.datamodels as dm
from jwst.extract_1d import Extract1dStep
from jwst.photom import PhotomStep
import pytest
import numpy as np


@pytest.fixture(scope='module')
def mock_niriss_full():
    model = dm.CubeModel((3, 3, 3))
    model.meta.instrument.name = 'NIRISS'
    model.meta.instrument.detector = 'NIS'
    model.meta.observation.date = '2023-07-22'
    model.meta.observation.time = '06:24:45.569'
    model.meta.instrument.name = 'NIRISS'
    model.meta.instrument.detector = 'NIS'
    model.meta.instrument.filter = 'CLEAR'
    model.meta.exposure.type = 'NIS_SOSS'
    model.meta.subarray.name = 'FULL'
    model.data = np.arange(27).reshape((3, 3, 3))
    return model


def test_expected_skip_niriss_soss_full(mock_niriss_full):

    with mock_niriss_full as model:
        result = Extract1dStep().process(model)
        result2 = PhotomStep().process(result)
        assert result2.meta.cal_step.photom == 'SKIPPED'
        assert result2.meta.cal_step.extract_1d == 'SKIPPED'
        assert np.all(result2.data == model.data)
