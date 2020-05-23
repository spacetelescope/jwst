"""Test aspects of CubeModel"""
import pytest

from jwst.datamodels import CubeModel


# Generate data
def make_container():
    """Create the container to test"""
    size = 5
    cube = CubeModel((size, size, size))
    cube.meta.target.proposer_name = 'JWST regression test'
    container = cube.to_container()
    return cube, container


cube, container = make_container()
models = [
    (cube, model)
    for model in container
]


@pytest.mark.parametrize('cube, container', [(cube, container)])
def test_container_shape(cube, container):
    """Test container shape"""
    assert len(container) == cube.shape[0]


@pytest.mark.parametrize('cube, model', models)
def test_meta(cube, model):
    """Test meta equivalency"""
    assert model.meta.target.proposer_name == cube.meta.target.proposer_name


@pytest.mark.parametrize('array', ['data', 'dq', 'err'])
@pytest.mark.parametrize('cube, model', models)
def test_shape(cube, model, array):
    """Test array shapes"""
    assert model[array].shape == cube[array][0].shape


@pytest.mark.parametrize('array', ['zeroframe', 'area', 'con', 'wht'])
@pytest.mark.parametrize('cube, model', models)
def test_nonexistant_arrays(cube, model, array):
    """Test for non-existant arrays"""
    with pytest.raises(AttributeError):
        model.getarray_noinit(array)
