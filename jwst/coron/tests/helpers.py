import numpy as np
from astropy.modeling.functional_models import Gaussian2D
from stdatamodels.jwst import datamodels

__all__ = [
    "nircam_model",
    "miri_model",
    "gaussian_source",
    "target_nircam",
    "psf_nircam",
    "target_miri",
    "psf_miri",
]


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


def target_nircam():
    """
    Make a NIRCam model with a Gaussian target, offset from center.

    Returns
    -------
    CubeModel
        The target model.
    """
    model = nircam_model()
    shape = model.shape[1:]
    y_mean = shape[0] / 2 + 1.0
    x_mean = shape[1] / 2 + 1.0
    amplitude = 5.0
    model.data += gaussian_source(shape, amplitude=amplitude, x_mean=x_mean, y_mean=y_mean)

    return model


def psf_nircam():
    """
    Make a NIRCam model with a Gaussian target at center.

    Returns
    -------
    CubeModel
        The PSF model.
    """
    model = nircam_model()
    model.data += gaussian_source(model.shape[1:])
    return model


def target_miri():
    """
    Make a MIRI model with a Gaussian target, offset from center.

    Returns
    -------
    CubeModel
        The target model.
    """
    model = miri_model()
    shape = model.shape[1:]
    y_mean = shape[0] / 2 + 1.0
    x_mean = shape[1] / 2 + 1.0
    amplitude = 5.0
    model.data += gaussian_source(shape, amplitude=amplitude, x_mean=x_mean, y_mean=y_mean)

    return model


def psf_miri():
    """
    Make a MIRI model with a Gaussian target at center.

    Returns
    -------
    CubeModel
        The PSF model.
    """
    model = miri_model()
    model.data += gaussian_source(model.shape[1:])
    return model
