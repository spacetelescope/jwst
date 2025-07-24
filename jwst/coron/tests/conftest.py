import numpy as np
import pytest
from astropy.modeling.functional_models import Gaussian2D
from stdatamodels.jwst import datamodels


def nircam_model():
    """
    Create a NRC_CORON cube model.

    Returns
    -------
    CubeModel
        The NIRCam coron model.
    """
    model = datamodels.CubeModel()
    model.meta.exposure.type = "NRC_CORON"
    model.meta.instrument.name = "NIRCAM"
    model.meta.instrument.filter = "F182M"
    model.meta.instrument.coronagraph = "MASKA335R"
    model.meta.observation.date = "2024-01-01"
    model.meta.observation.time = "00:00:00"
    model.meta.subarray.name = "SUB320A335R"

    shape = (3, 320, 320)
    model.data = np.zeros(shape, dtype=np.float32)
    model.err = np.full(shape, 0.1, dtype=np.float32)
    model.dq = np.zeros(shape, dtype=np.uint32)
    return model


def miri_model():
    """
    Create a MIR_LYOT cube model.

    Returns
    -------
    CubeModel
        The MIRI coron model.
    """
    model = datamodels.CubeModel()
    model.meta.exposure.type = "MIR_LYOT"
    model.meta.instrument.name = "MIRI"
    model.meta.instrument.filter = "F2300C"
    model.meta.instrument.coronagraph = "LYOT_2300"
    model.meta.observation.date = "2024-01-01"
    model.meta.observation.time = "00:00:00"
    model.meta.subarray.name = "MASKLYOT"

    shape = (3, 304, 320)
    model.data = np.zeros(shape, dtype=np.float32)
    model.err = np.full(shape, 0.1, dtype=np.float32)
    model.dq = np.zeros(shape, dtype=np.uint32)
    return model


def gaussian_source(shape, amplitude=10.0, x_mean=None, y_mean=None, x_stddev=30.0, y_stddev=30.0):
    """
    Model a Gaussian source.

    Returns
    -------
    ndarray
        An array of the expected shape, containing the modeled Gaussian.
    """
    y, x = np.mgrid[: shape[0], : shape[1]]
    if y_mean is None:
        y_mean = shape[0] / 2
    if x_mean is None:
        x_mean = shape[1] / 2
    source = Gaussian2D(amplitude, x_mean, y_mean, x_stddev, y_stddev)
    return source(y, x)


@pytest.fixture()
def target_model():
    """
    Make a NIRCam model with a Gaussian target, offset from center.

    Yields
    ------
    CubeModel
        The target model.
    """
    model = nircam_model()
    shape = model.shape[1:]
    y_mean = shape[0] / 2 + 1.0
    x_mean = shape[1] / 2 + 1.0
    amplitude = 5.0
    model.data += gaussian_source(shape, amplitude=amplitude, x_mean=x_mean, y_mean=y_mean)

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
    model = nircam_model()
    model.data += gaussian_source(model.shape[1:])
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
    model = miri_model()
    shape = model.shape[1:]
    y_mean = shape[0] / 2 + 1.0
    x_mean = shape[1] / 2 + 1.0
    amplitude = 5.0
    model.data += gaussian_source(shape, amplitude=amplitude, x_mean=x_mean, y_mean=y_mean)

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
    model = miri_model()
    model.data += gaussian_source(model.shape[1:])
    yield model
    model.close()


@pytest.fixture()
def target_image(target_model):
    """
    Make an image model from the first slice of the target cube.

    Yields
    ------
    ImageModel
        The target model.
    """
    target_image = datamodels.ImageModel()
    target_image.data = target_model.data[0]
    target_image.err = target_model.err[0]
    target_image.dq = target_model.dq[0]
    target_image.update(target_model)
    yield target_image
    target_image.close()


@pytest.fixture()
def psf_image(psf_model):
    """
    Make an image model from the first slice of the psf cube.

    Yields
    ------
    ImageModel
        The psf model.
    """
    psf_image = datamodels.ImageModel()
    psf_image.data = psf_model.data[0]
    psf_image.err = psf_model.err[0]
    psf_image.dq = psf_model.dq[0]
    psf_image.update(psf_model)
    yield psf_image
    psf_image.close()
