"""Unit tests for firstframe  module and step."""

import pytest
from stdatamodels.jwst.datamodels import ImageModel

from jwst.firstframe.firstframe_step import FirstFrameStep
from jwst.tests.helpers import _help_pytest_warns


def test_firstframe_deprecated(create_miri_model):
    model = create_miri_model()
    image_list = []
    for i in range(model.data.shape[0]):
        image_list.append(ImageModel(model.data[i, :, :]))

    with _help_pytest_warns(), pytest.deprecated_call():
        result = FirstFrameStep.call(image_list)
        assert result is not image_list
        assert result.meta.cal_step.firstframe == "COMPLETE"
        for model in image_list:
            assert result is not model
            assert model.meta.cal_step.firstframe is None
