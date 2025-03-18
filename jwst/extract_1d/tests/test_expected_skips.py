import stdatamodels.jwst.datamodels as dm
from jwst.combine_1d import Combine1dStep
from jwst.extract_1d import Extract1dStep
from jwst.photom import PhotomStep
import pytest
import numpy as np


@pytest.fixture(scope="module")
def mock_niriss_full():
    model = dm.CubeModel((3, 3, 3))
    model.meta.instrument.name = "NIRISS"
    model.meta.instrument.detector = "NIS"
    model.meta.observation.date = "2023-07-22"
    model.meta.observation.time = "06:24:45.569"
    model.meta.instrument.name = "NIRISS"
    model.meta.instrument.detector = "NIS"
    model.meta.instrument.filter = "CLEAR"
    model.meta.exposure.type = "NIS_SOSS"
    model.meta.subarray.name = "FULL"
    model.data = np.arange(27).reshape((3, 3, 3))
    return model


@pytest.fixture(scope="module")
def mock_niriss_f277w():
    model = dm.CubeModel((3, 3, 3))
    model.meta.instrument.name = "NIRISS"
    model.meta.instrument.detector = "NIS"
    model.meta.observation.date = "2023-07-22"
    model.meta.observation.time = "06:24:45.569"
    model.meta.instrument.name = "NIRISS"
    model.meta.instrument.detector = "NIS"
    model.meta.instrument.filter = "F277W"
    model.meta.exposure.type = "NIS_SOSS"
    model.meta.subarray.name = "FULL"
    model.data = np.arange(27).reshape((3, 3, 3))
    return model


def test_expected_skip_niriss_soss_full(mock_niriss_full):
    with mock_niriss_full as model:
        result = Extract1dStep().process(model)
        result2 = PhotomStep().process(result)
        result3 = Combine1dStep().process(result)
        assert result2.meta.cal_step.photom == "SKIPPED"
        assert result2.meta.cal_step.extract_1d == "SKIPPED"
        assert result3.meta.cal_step.combine_1d == "SKIPPED"
        assert np.all(result2.data == model.data)


def test_expected_skip_niriss_soss_f277w(mock_niriss_f277w):
    with mock_niriss_f277w as model:
        result = Extract1dStep().process(model)
        result2 = PhotomStep().process(result)
        result3 = Combine1dStep().process(result)
        assert result2.meta.cal_step.photom == "SKIPPED"
        assert result2.meta.cal_step.extract_1d == "SKIPPED"
        assert result3.meta.cal_step.combine_1d == "SKIPPED"
        assert np.all(result2.data == model.data)


def test_expected_skip_multi_int_multi_slit():
    model = dm.MultiSlitModel()
    model.slits.append(dm.SlitModel(np.zeros((10, 10, 10))))
    result = Extract1dStep().process(model)
    assert result.meta.cal_step.extract_1d == "SKIPPED"
    assert np.all(result.slits[0].data == model.slits[0].data)
    model.close()
    result.close()


def test_expected_skip_unexpected_model():
    model = dm.MultiExposureModel()
    result = Extract1dStep().process(model)
    assert result.meta.cal_step.extract_1d == "SKIPPED"
    model.close()
    result.close()
