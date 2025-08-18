import pytest
from stdatamodels.jwst import datamodels

from jwst.coron.tests.helpers import psf_miri, psf_nircam, target_miri, target_nircam


@pytest.fixture()
def target_model():
    """
    Make a NIRCam model with a Gaussian target, offset from center.

    Yields
    ------
    CubeModel
        The target model.
    """
    model = target_nircam()
    yield model
    model.close()


@pytest.fixture()
def psf_model():
    """
    Make a NIRCam model with a Gaussian target at center.

    Yields
    ------
    CubeModel
        The PSF model.
    """
    model = psf_nircam()
    yield model
    model.close()


@pytest.fixture()
def target_model_miri():
    """
    Make a MIRI model with a Gaussian target, offset from center.

    Yields
    ------
    CubeModel
        The target model.
    """
    model = target_miri()
    yield model
    model.close()


@pytest.fixture()
def psf_model_miri():
    """
    Make a MIRI model with a Gaussian target at center.

    Yields
    ------
    CubeModel
        The PSF model.
    """
    model = psf_miri()
    yield model
    model.close()


@pytest.fixture()
def target_image():
    """
    Make an image model from the first slice of the target cube.

    Yields
    ------
    ImageModel
        The target model.
    """
    cube = target_nircam()
    image = datamodels.ImageModel()
    image.data = cube.data[0]
    image.err = cube.err[0]
    image.dq = cube.dq[0]
    image.update(cube)
    yield image
    image.close()


@pytest.fixture()
def psf_image():
    """
    Make an image model from the first slice of the psf cube.

    Yields
    ------
    ImageModel
        The psf model.
    """
    cube = psf_nircam()
    image = datamodels.ImageModel()
    image.data = cube.data[0]
    image.err = cube.err[0]
    image.dq = cube.dq[0]
    image.update(cube)
    yield image
    image.close()
