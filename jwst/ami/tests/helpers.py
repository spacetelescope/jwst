"""Useful constants and functions for AMI tests."""

import numpy as np
import stdatamodels.jwst.datamodels as dm

PXSC_DEG = 65.6 / (60.0 * 60.0 * 1000)
PXSC_RAD = PXSC_DEG * np.pi / (180)
PXSC_MAS = PXSC_DEG * 3600 * 1000

__all__ = ["example_model", "circular_pupil"]


def example_model():
    """
    Create a simple CubeModel simulating input to the ami3 pipeline.

    Returns
    -------
    CubeModel
        A simple CubeModel simulating NIRISS AMI data.
    """
    model = dm.CubeModel((5, 81, 81))

    # make a simple data array, pure noise but with one bad pixel
    rng = np.random.default_rng(0)
    data = rng.normal(size=model.data.shape).astype(np.float32)
    data[0, 20, 20] = 100

    # add a "real" source that is not time varying
    # leave it a bit off-center to test centering in instrument_data.NIRISS.read_data_model
    data[:, 35, 35] = 10.0

    model.data = data
    model.meta.instrument.name = "NIRISS"
    model.meta.instrument.filter = "F277W"
    model.meta.subarray.name = "SUB80"
    model.meta.observation.date = "2021-12-26"
    model.meta.observation.time = "00:00:00"
    model.meta.target.proposer_name = ""
    model.meta.program.pi_name = "someone"
    model.meta.target.catalog_name = ""
    model.meta.visit.start_time = "2022-06-05 12:15:41.5020000"
    model.meta.wcsinfo.roll_ref = 171.8779402866089
    model.meta.wcsinfo.v3yangle = 0.56126717
    model.meta.wcsinfo.vparity = 1
    model.meta.filename = "test_calints.fits"
    model.meta.instrument.pupil = "NRM"
    model.meta.exposure.type = "NIS_AMI"

    # Pointing information copied from regression test data
    model.meta.target.ra = 82.18735041666667
    model.meta.target.dec = -65.44798333333335
    model.meta.target.proper_motion_ra = 0.02914965139708291
    model.meta.target.proper_motion_dec = 0.1644210871070447
    model.meta.target.ra_uncertainty = 0.000125653517585233
    model.meta.target.dec_uncertainty = 0.000153823202400172

    return model


def circular_pupil():
    """
    Make a simple circular pupil mask.

    Returns
    -------
    ndarray
        A 1024x1024 array with a circular pupil of radius 0.2x the array shape.
    """
    shape = (1024, 1024)
    r = 0.2
    x = np.linspace(-1, 1, shape[0])
    y = np.linspace(-1, 1, shape[1])
    xx, yy = np.meshgrid(x, y)
    rr = np.sqrt(xx**2 + yy**2)
    pupil = np.zeros(shape)
    pupil[rr < r] = 1
    return pupil
