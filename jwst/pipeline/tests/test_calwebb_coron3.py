import pytest

from stdatamodels.jwst.datamodels import CubeModel

from jwst.pipeline.calwebb_coron3 import to_container
import numpy as np


# Generate data
def make_container():
    """Create the container to test"""
    size = 5
    cube = CubeModel((size, size, size))
    cube.meta.target.proposer_name = "JWST regression test"
    container = to_container(cube)
    return cube, container


cube, container = make_container()
models = [(cube, model) for model in container]


def test_to_container():
    """Cover bug where IndexError would be raised when area extension of CubeModel
    has shape (x,y) instead of (nints,x,y). In this case area extension should
    be copied to each ImageModel in the ModelContainer.
    """
    shp = (10, 5, 5)
    cube = CubeModel(shp)
    extensions_3d = [
        "data",
        "dq",
        "err",
    ]
    for extension in extensions_3d:
        setattr(cube, extension, np.random.rand(*shp))
    cube.area = np.random.rand(5, 5)

    container = to_container(cube)
    for i, model in enumerate(container):
        for extension in extensions_3d:
            assert np.all(getattr(model, extension) == getattr(cube, extension)[i])
        assert hasattr(model, "area")
        assert np.all(model.area == cube.area)


@pytest.mark.parametrize("cube, container", [(cube, container)])
def test_container_shape(cube, container):
    """Test container shape"""
    assert len(container) == cube.shape[0]


@pytest.mark.parametrize("cube, model", models)
def test_meta(cube, model):
    """Test meta equivalency"""
    assert model.meta.target.proposer_name == cube.meta.target.proposer_name


@pytest.mark.parametrize("array", ["data", "dq", "err"])
@pytest.mark.parametrize("cube, model", models)
def test_shape(cube, model, array):
    """Test array shapes"""
    assert model[array].shape == cube[array][0].shape


@pytest.mark.parametrize("array", ["zeroframe", "area", "con", "wht"])
@pytest.mark.parametrize("cube, model", models)
def test_nonexistent_arrays(cube, model, array):
    """Test for non-existent arrays"""
    with pytest.raises(AttributeError):
        model.getarray_noinit(array)
