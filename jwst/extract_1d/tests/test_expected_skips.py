import numpy as np
import stdatamodels.jwst.datamodels as dm

from jwst.combine_1d import Combine1dStep
from jwst.extract_1d import Extract1dStep
from jwst.extract_1d.tests import helpers
from jwst.photom import PhotomStep


def test_expected_skip_niriss_soss_full():
    with helpers.mock_niriss_soss_full_func() as model:
        result = Extract1dStep().process(model)
        result2 = PhotomStep().process(result)
        result3 = Combine1dStep().process(result)
        assert result2.meta.cal_step.photom == "SKIPPED"
        assert result2.meta.cal_step.extract_1d == "SKIPPED"
        assert result3.meta.cal_step.combine_1d == "SKIPPED"
        assert np.all(result2.data == model.data)

        # make sure input is not modified
        assert result is not model
        assert result2 is not model
        assert result3 is not model
        assert model.meta.cal_step.extract_1d is None
        assert model.meta.cal_step.photom is None
        assert model.meta.cal_step.combine_1d is None


def test_expected_skip_niriss_soss_f277w():
    with helpers.mock_niriss_soss_f277w_func() as model:
        result = Extract1dStep().process(model)
        result2 = PhotomStep().process(result)
        result3 = Combine1dStep().process(result)
        assert result2.meta.cal_step.photom == "SKIPPED"
        assert result2.meta.cal_step.extract_1d == "SKIPPED"
        assert result3.meta.cal_step.combine_1d == "SKIPPED"
        assert np.all(result2.data == model.data)

        # make sure input is not modified
        assert result is not model
        assert result2 is not model
        assert result3 is not model
        assert model.meta.cal_step.extract_1d is None
        assert model.meta.cal_step.photom is None
        assert model.meta.cal_step.combine_1d is None


def test_expected_skip_multi_int_multi_slit():
    model = dm.MultiSlitModel()
    model.slits.append(dm.SlitModel(np.zeros((10, 10, 10))))
    result = Extract1dStep().process(model)
    assert result.meta.cal_step.extract_1d == "SKIPPED"
    assert np.all(result.slits[0].data == model.slits[0].data)

    # make sure input is not modified
    assert result is not model
    assert model.meta.cal_step.extract_1d is None

    model.close()
    result.close()


def test_expected_skip_unexpected_model():
    model = dm.MultiExposureModel()
    result = Extract1dStep().process(model)
    assert result.meta.cal_step.extract_1d == "SKIPPED"

    # make sure input is not modified
    assert result is not model
    assert model.meta.cal_step.extract_1d is None

    model.close()
    result.close()
