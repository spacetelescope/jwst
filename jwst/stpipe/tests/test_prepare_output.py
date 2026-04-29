import numpy as np
import pytest
from astropy.utils.data import get_pkg_data_filename
from stdatamodels.jwst import datamodels

from jwst.datamodels import ModelContainer, ModelLibrary
from jwst.stpipe.tests.steps import (
    PrepareOutputAsRampStep,
    PrepareOutputAsTypeStep,
    PrepareOutputForceCopyNoOpenStep,
    PrepareOutputForceCopyPipeline,
    PrepareOutputForceCopyStep,
    PrepareOutputNoCopyStep,
    PrepareOutputNoOpenStep,
    PrepareOutputPipeline,
    PrepareOutputStep,
)


@pytest.fixture()
def mock_copy(monkeypatch):
    """Monkeypatch the ImageModel copy to raise an error if called."""

    def _mock(*args, **kwargs):
        raise RuntimeError("Copy was called")

    monkeypatch.setattr(datamodels.image.ImageModel, "copy", _mock)


@pytest.mark.parametrize("step", [PrepareOutputStep, PrepareOutputForceCopyStep])
def test_prepare_output_from_str(step):
    input_data = get_pkg_data_filename("data/science.fits", package="jwst.stpipe.tests")
    result = step.call(input_data)

    # Input is str, output is model
    assert result is not input_data
    assert isinstance(input_data, str)
    assert isinstance(result, datamodels.ImageModel)
    assert result.meta.cal_step.prepare_output == "COMPLETE"


def test_prepare_output_from_str_no_open():
    input_data = get_pkg_data_filename("data/science.fits", package="jwst.stpipe.tests")

    # Step does not call datamodels.open for the input string
    result = PrepareOutputNoOpenStep.call(input_data)

    assert result is input_data
    assert isinstance(input_data, str)
    assert isinstance(result, str)


def test_prepare_output_from_str_force_copy_no_open(monkeypatch):
    input_data = get_pkg_data_filename("data/science.fits", package="jwst.stpipe.tests")

    # Raises a type error: input is not opened so it can't be copied
    with pytest.raises(TypeError, match="Copy is not possible"):
        PrepareOutputForceCopyNoOpenStep.call(input_data)


@pytest.mark.parametrize("step", [PrepareOutputStep, PrepareOutputForceCopyStep])
def test_prepare_output_from_model(step):
    input_data = datamodels.ImageModel()
    result = step.call(input_data)

    # Output is a copy of the input
    assert result is not input_data
    assert isinstance(result, datamodels.ImageModel)
    assert result.meta.cal_step.prepare_output == "COMPLETE"
    assert not input_data.meta.cal_step.hasattr("prepare_output")


def test_prepare_output_from_model_no_copy():
    input_data = datamodels.ImageModel()
    result = PrepareOutputNoCopyStep.call(input_data)

    # Output is the same as the input
    assert result is input_data
    assert isinstance(result, datamodels.ImageModel)

    assert result.meta.cal_step.prepare_output == "COMPLETE"
    assert input_data.meta.cal_step.prepare_output == "COMPLETE"


@pytest.mark.parametrize("step", [PrepareOutputStep, PrepareOutputForceCopyStep])
def test_prepare_output_from_container(step):
    container = ModelContainer([datamodels.ImageModel()])
    input_model = container[0]
    result = step.call(container)

    # Output is a copy of the input
    assert result is not container
    assert result[0] is not input_model
    assert isinstance(result, ModelContainer)
    assert isinstance(result[0], datamodels.ImageModel)
    assert result[0].meta.cal_step.prepare_output == "COMPLETE"
    assert not input_model.meta.cal_step.hasattr("prepare_output")


def test_prepare_output_from_container_no_copy():
    container = ModelContainer([datamodels.ImageModel()])
    input_model = container[0]
    result = PrepareOutputNoCopyStep.call(container)

    # Output container is the same as the input
    assert result is container
    assert result[0] is input_model
    assert result[0].meta.cal_step.prepare_output == "COMPLETE"
    assert input_model.meta.cal_step.prepare_output == "COMPLETE"


@pytest.mark.parametrize("step", [PrepareOutputStep, PrepareOutputForceCopyStep])
def test_prepare_output_from_list_of_model(step):
    input_model = datamodels.ImageModel()
    input_data = [input_model]
    result = step.call(input_data)

    # Output is a ModelContainer with a copy of the input
    assert isinstance(result, ModelContainer)
    assert result[0] is not input_model
    assert isinstance(result[0], datamodels.ImageModel)
    assert result[0].meta.cal_step.prepare_output == "COMPLETE"
    assert not input_model.meta.cal_step.hasattr("prepare_output")


def test_prepare_output_from_list_of_model_no_copy():
    input_model = datamodels.ImageModel()
    input_data = [input_model]
    result = PrepareOutputNoCopyStep.call(input_data)

    # Output is a ModelContainer with a shallow copy of the input
    assert isinstance(result, ModelContainer)
    assert result[0] is not input_model
    assert result[0].meta.cal_step.prepare_output == "COMPLETE"
    assert input_model.meta.cal_step.prepare_output == "COMPLETE"


def test_prepare_output_from_list_of_str():
    input_data = [get_pkg_data_filename("data/science.fits", package="jwst.stpipe.tests")]
    result = PrepareOutputStep.call(input_data)

    # Input is list of str, output is ModelContainer
    assert result is not input_data
    assert isinstance(input_data, list)
    assert isinstance(input_data[0], str)
    assert isinstance(result, ModelContainer)
    assert isinstance(result[0], datamodels.ImageModel)
    assert result[0].meta.cal_step.prepare_output == "COMPLETE"


def test_prepare_output_pipeline_from_str(mock_copy):
    input_data = get_pkg_data_filename("data/science.fits", package="jwst.stpipe.tests")

    # Should run with no errors: input data is not copied after being opened
    result = PrepareOutputPipeline.call(input_data)
    assert isinstance(result, datamodels.ImageModel)
    assert result.meta.cal_step.prepare_output == "COMPLETE"


def test_prepare_output_pipeline_force_copy(mock_copy):
    input_data = get_pkg_data_filename("data/science.fits", package="jwst.stpipe.tests")

    # Should generate an error from the mocked copy function: copy is
    # called by one of the steps.
    with pytest.raises(RuntimeError, match="Copy was called"):
        PrepareOutputForceCopyPipeline.call(input_data)


def test_prepare_output_pipeline_from_model():
    input_data = datamodels.ImageModel()

    result = PrepareOutputPipeline.call(input_data)

    # Output is a copy of the input
    assert result is not input_data
    assert isinstance(result, datamodels.ImageModel)
    assert result.meta.cal_step.prepare_output == "COMPLETE"
    assert not input_data.meta.cal_step.hasattr("prepare_output")


def test_prepare_output_pipeline_from_container():
    container = ModelContainer([datamodels.ImageModel()])
    input_model = container[0]
    result = PrepareOutputPipeline.call(container)

    # Output is a copy of the input
    assert result is not container
    assert result[0] is not input_model
    assert isinstance(result, ModelContainer)
    assert isinstance(result[0], datamodels.ImageModel)
    assert result[0].meta.cal_step.prepare_output == "COMPLETE"
    assert not input_model.meta.cal_step.hasattr("prepare_output")


def test_prepare_output_pipeline_from_list_of_model():
    input_model = datamodels.ImageModel()
    input_data = [input_model]
    result = PrepareOutputPipeline.call(input_data)

    # Output is a ModelContainer with a copy of the input
    assert isinstance(result, ModelContainer)
    assert result[0] is not input_model
    assert isinstance(result[0], datamodels.ImageModel)
    assert result[0].meta.cal_step.prepare_output == "COMPLETE"
    assert not input_model.meta.cal_step.hasattr("prepare_output")


def test_prepare_output_pipeline_from_list_of_str(mock_copy):
    input_data = [get_pkg_data_filename("data/science.fits", package="jwst.stpipe.tests")]

    # Copy is never called (mock_copy produces no errors)
    result = PrepareOutputPipeline.call(input_data)

    # Input is still list of str, output is ModelContainer
    assert result is not input_data
    assert isinstance(input_data, list)
    assert isinstance(input_data[0], str)
    assert isinstance(result, ModelContainer)
    assert isinstance(result[0], datamodels.ImageModel)
    assert result[0].meta.cal_step.prepare_output == "COMPLETE"


@pytest.mark.parametrize("step", [PrepareOutputStep, PrepareOutputPipeline])
@pytest.mark.parametrize(
    "input_data",
    [
        datamodels.ImageModel(),
        get_pkg_data_filename("data/science.fits", package="jwst.stpipe.tests"),
    ],
)
def test_prepare_output_from_library(mock_copy, step, input_data):
    library = ModelLibrary([input_data])

    # Copies are never made for libraries
    result = step.call(library)

    assert result is library
    with library:
        model = library.borrow(0)
        if isinstance(input_data, datamodels.ImageModel):
            assert model is input_data
        assert model.meta.cal_step.prepare_output == "COMPLETE"
        library.shelve(model, modify=False)


def test_prepare_output_open_as_type_from_model():
    input_data = datamodels.ImageModel((10, 10))
    result = PrepareOutputAsTypeStep.call(input_data)

    # Output is a deep copy of the input
    assert result is not input_data
    assert result.data is not input_data.data

    # Output has the specified type
    assert type(input_data) is datamodels.ImageModel
    assert type(result) is datamodels.IFUImageModel

    assert result.meta.cal_step.prepare_output == "COMPLETE"
    assert not input_data.meta.cal_step.hasattr("prepare_output")


def test_prepare_output_open_as_type_from_str():
    input_data = get_pkg_data_filename("data/science.fits", package="jwst.stpipe.tests")
    result = PrepareOutputAsTypeStep.call(input_data)

    # Input is str, output is expected model type
    assert result is not input_data
    assert isinstance(input_data, str)
    assert isinstance(result, datamodels.IFUImageModel)
    assert result.meta.cal_step.prepare_output == "COMPLETE"


@pytest.mark.parametrize("from_str", [True, False])
@pytest.mark.parametrize(
    "expected_type",
    [datamodels.RampModel, datamodels.SuperstripeRampModel, datamodels.GuiderRawModel],
)
def test_prepare_output_open_as_ramp_from_model(tmp_path, expected_type, from_str):
    input_data = datamodels.RampModel()
    input_data.data = np.zeros((6, 3, 4, 5), dtype=np.float32)
    if expected_type is datamodels.SuperstripeRampModel:
        input_data.meta.subarray.num_superstripe = 3
    elif expected_type is datamodels.GuiderRawModel:
        input_data.meta.exposure.type = "FGS_ID-IMAGE"

    if from_str:
        filename = str(tmp_path / "test_ramp.fits")
        input_data.save(filename)
        input_data.close()
        input_data = filename

    result = PrepareOutputAsRampStep.call(input_data)

    if from_str:
        # Input is str, output is expected model type
        assert result is not input_data
        assert isinstance(input_data, str)
    else:
        # Output is a deep copy of the input
        assert result is not input_data
        assert result.data is not input_data.data
        assert type(input_data) is datamodels.RampModel
        # The cal_step status is not updated in the input
        assert not input_data.meta.cal_step.hasattr("prepare_output")

    # Output has the specified type and cal_step status is updated
    assert type(result) is expected_type
    assert result.meta.cal_step.prepare_output == "COMPLETE"


def test_prepare_output_open_as_ramp_incompatible_pixeldq():
    # Make a ramp model from a shape, so it has a 2D pixeldq
    input_data = datamodels.RampModel((6, 3, 4, 5))

    # Add superstripe configuration
    input_data.meta.subarray.num_superstripe = 3

    # Input data can't be converted. The 2D pixeldq needs to be removed.
    with pytest.raises(TypeError, match="not compatible"):
        PrepareOutputAsRampStep.call(input_data)

    # Now it should succeed
    input_data.pixeldq = None
    result = PrepareOutputAsRampStep.call(input_data)
    assert type(result) is datamodels.SuperstripeRampModel
    assert result.data.shape == (6, 3, 4, 5)
    assert result.pixeldq.shape == (3, 4, 5)


def test_prepare_output_open_as_ramp_level1b():
    # Make a level1b model from a shape, so it has a uint16 data
    input_data = datamodels.Level1bModel((6, 3, 4, 5))

    # Input data can't be converted because of the data type
    with pytest.raises(TypeError, match="not compatible"):
        PrepareOutputAsRampStep.call(input_data)
