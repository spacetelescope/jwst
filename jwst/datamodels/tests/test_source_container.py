import numpy as np
import pytest
from stdatamodels.jwst import datamodels

from jwst.datamodels import SourceModelContainer


@pytest.fixture()
def multi_exp():
    exposures = [datamodels.SlitModel(), datamodels.SlitModel(), datamodels.SlitModel()]
    model = datamodels.MultiExposureModel()
    for i, exposure in enumerate(exposures):
        exposure.meta.filename = f"test_{i}"
        exposure.data = np.full((10, 10), i)
        model.exposures.append(exposure)
    return model


def test_source_container_invalid_init():
    with pytest.raises(TypeError, match="None cannot initialize a SourceModelContainer"):
        SourceModelContainer(None)


def test_source_container_from_empty_multiexp():
    model = datamodels.MultiExposureModel()
    container = SourceModelContainer(model)
    assert len(container._models) == 0
    assert container.multiexposure is model


def test_source_container_from_multiexp(multi_exp):
    container = SourceModelContainer(multi_exp)
    assert len(container._models) == len(multi_exp.exposures)
    assert container.multiexposure is multi_exp
    for in_exp, out_exp in zip(multi_exp.exposures, container._models):
        assert in_exp is not out_exp
        assert in_exp == out_exp


def test_source_container_from_source_container(multi_exp):
    container1 = SourceModelContainer(multi_exp)

    # Make a shallow copy from another container
    container2 = SourceModelContainer(container1)

    assert container1 is not container2
    for in_exp, out_exp in zip(container1._models, container2._models):
        assert in_exp is not out_exp
        assert in_exp == out_exp
    assert container2._multiexposure is container1._multiexposure


def test_save(tmp_path, multi_exp):
    container = SourceModelContainer(multi_exp)
    out_path = tmp_path / "test_multiexp.fits"
    container.save(str(out_path))
    assert out_path.exists()

    with datamodels.open(out_path) as model:
        assert isinstance(model, datamodels.MultiExposureModel)
        assert len(model.exposures) == len(container._models)
        for in_exp, out_exp in zip(container._models, model.exposures):
            assert in_exp is not out_exp
            assert in_exp.meta.filename == out_exp.meta.filename
            assert np.all(in_exp.data == out_exp.data)

            # Equality comparison fails when underlying instance is not a shallow copy
            with pytest.raises(ValueError, match="truth value of an array"):
                assert in_exp == out_exp


def test_save_func(tmp_path, multi_exp):
    def save_func(*args, **kwargs):
        out_path = tmp_path / "test_save_func.txt"
        with open(out_path, "w") as fh:
            fh.write("test\n")

    container = SourceModelContainer(multi_exp)
    out_path = tmp_path / "test_save_func.txt"
    container.save(save_model_func=save_func)
    assert out_path.exists()

    with open(out_path) as fh:
        lines = fh.readlines()
    assert lines == ["test\n"]


def test_copy(multi_exp):
    # make a deep copy of a container
    container = SourceModelContainer(multi_exp)
    container_copy = container.copy()

    assert len(container._models) == len(container_copy._models)
    assert container._multiexposure is not container_copy._multiexposure
    assert len(container._multiexposure.exposures) == len(container_copy._multiexposure.exposures)
    for in_exp, out_exp in zip(container._models, container_copy._models, strict=True):
        assert in_exp is not out_exp
        assert in_exp.meta.filename == out_exp.meta.filename
        assert np.all(in_exp.data == out_exp.data)

        # Equal comparison fails when underlying instance is not a shallow copy
        with pytest.raises(ValueError, match="truth value of an array"):
            assert in_exp == out_exp
