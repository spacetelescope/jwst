import numpy as np

from jwst import datamodels

__all__ = ["make_guider_image"]


def make_guider_image():
    """
    Generate science image.

    Returns
    -------
    image : `stdatamodels.jwst.datamodels.GuiderRawModel`
        The guider image.
    """
    image = datamodels.GuiderRawModel()
    image.meta.filename = "test_guider_rate.fits"

    image.meta.instrument.name = "FGS"
    image.meta.instrument.detector = "GUIDER1"
    image.meta.observation.date = "2016-04-07"
    image.meta.observation.time = "14:44:57"
    image.meta.exposure.frame_time = 234.3423235
    image.meta.exposure.ngroups = 4
    image.meta.exposure.group_time = 465.643643
    image.meta.exposure.type = "FGS_FINEGUIDE"

    rng = np.random.default_rng(seed=42)
    image.data = rng.random((4, 10, 10, 10))
    image.meta.subarray.xstart = 1226
    image.meta.subarray.ystart = 209
    image.meta.subarray.xsize = 10
    image.meta.subarray.ysize = 10

    return image
