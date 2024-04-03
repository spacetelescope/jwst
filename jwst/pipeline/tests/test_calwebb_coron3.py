import pytest

from stdatamodels.jwst.datamodels import CubeModel

from jwst.pipeline.calwebb_coron3 import to_models


CUBE_SIZE = 5

@pytest.fixture
def cube():
    return CubeModel((CUBE_SIZE, CUBE_SIZE, CUBE_SIZE))


@pytest.fixture
def cube_models(cube):
    """Create the container to test"""
    cube.meta.target.proposer_name = 'JWST regression test'
    return to_models(cube)


@pytest.fixture(params=range(CUBE_SIZE))
def single_cube_model(request, cube_models):
    return cube_models[request.param]


def test_container_shape(cube, cube_models):
    """Test container shape"""
    assert len(cube_models) == cube.shape[0]


def test_meta(cube, single_cube_model):
    """Test meta equivalency"""
    assert single_cube_model.meta.target.proposer_name == cube.meta.target.proposer_name


@pytest.mark.parametrize('array', ['data', 'dq', 'err'])
def test_shape(cube, single_cube_model, array):
    """Test array shapes"""
    assert single_cube_model[array].shape == cube[array][0].shape


@pytest.mark.parametrize('array', ['zeroframe', 'area', 'con', 'wht'])
def test_nonexistent_arrays(cube, single_cube_model, array):
    """Test for non-existant arrays"""
    with pytest.raises(AttributeError):
        single_cube_model.getarray_noinit(array)
