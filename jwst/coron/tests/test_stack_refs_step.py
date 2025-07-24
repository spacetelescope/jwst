import pytest
from stdatamodels.jwst import datamodels

from jwst.coron.stack_refs_step import StackRefsStep
from jwst.datamodels import ModelContainer


@pytest.mark.parametrize("use_container", [True, False])
def test_stack_refs_step(target_model, psf_model, use_container):
    input_data = [target_model, psf_model]
    if use_container:
        input_data = ModelContainer(input_data)

    result = StackRefsStep.call(input_data)

    expected_shape = (
        target_model.shape[0] + psf_model.shape[0],
        target_model.shape[1],
        target_model.shape[2],
    )
    assert isinstance(result, datamodels.CubeModel)
    assert result.data.shape == expected_shape

    assert result.meta.cal_step.stack_psfs == "COMPLETE"

    # Make sure input was not modified
    assert result is not input_data
    assert result is not target_model
    assert result is not psf_model
    assert target_model.meta.cal_step.stack_psfs is None
    assert psf_model.meta.cal_step.stack_psfs is None


def test_stack_refs_shape_mismatch():
    model1 = datamodels.CubeModel((3, 10, 10))
    model2 = datamodels.CubeModel((2, 5, 5))
    with pytest.raises(ValueError, match="must have the same x/y dimensions"):
        StackRefsStep.call([model1, model2])
