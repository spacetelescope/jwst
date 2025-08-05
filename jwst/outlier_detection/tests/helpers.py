import numpy as np
from stdatamodels.jwst import datamodels

from jwst.assign_wcs import AssignWcsStep
from jwst.resample.resample import ResampleImage

EXPTYPE_TO_INSTRUMENT = {
    "MIR": "MIRI",
    "NRC": "NIRCAM",
    "NIS": "NIRISS",
    "NRS": "NIRSPEC",
    "FGS": "FGS",
}
SHAPE = (21, 20)
BACKGROUND = 1.5
SIGMA = 0.02
SIGNAL = 7.0
SIGNAL_LOC = (7, 7)

OUTLIER_DO_NOT_USE = np.bitwise_or(
    datamodels.dqflags.pixel["DO_NOT_USE"], datamodels.dqflags.pixel["OUTLIER"]
)


def mock_data():
    """
    Make some mock data with a "real" source at 7,7.

    Returns
    -------
    data, err : ndarray
        Data and error arrays.
    """
    rng = np.random.default_rng(99)
    j, k = SIGNAL_LOC
    data = rng.normal(loc=BACKGROUND, size=SHAPE, scale=SIGMA)
    err = np.zeros(SHAPE) + SIGMA
    data[j, k] += SIGNAL
    # update the noise for this source to include the photon/measurement noise
    err[j, k] = np.sqrt(SIGMA**2 + SIGNAL)
    return data, err


def assign_wcs_to_models(models, exptype, tsovisit, detector="ANY"):
    """
    Assign the same WCS to all models.

    Parameters
    ----------
    models : list or ModelContainer
        Models to update.
    exptype : str
        Exposure type to set in metadata.
    tsovisit : bool
        TSO visit status to set.
    detector : str, optional
        Detector to set.

    Returns
    -------
    models : list or ModelContainer
        Updated models.
    """
    for m in models:
        m.meta.exposure.type = exptype
        m.meta.instrument.name = EXPTYPE_TO_INSTRUMENT[exptype.split("_")[0]]
        m.meta.instrument.detector = detector
        m.meta.visit.tsovisit = tsovisit

    # Only need to call AssignWCS once because all are identical
    model = AssignWcsStep.call(models[0])
    wcs = model.meta.wcs
    wcsinfo = model.meta.wcsinfo

    for m in models:
        m.meta.wcs = wcs
        m.meta.wcsinfo = wcsinfo
    return models


def container_to_cube(container):
    """
    Convert a test container to a cube.

    Parameters
    ----------
    container : ModelContainer
        Models to stack.

    Returns
    -------
    CubeModel
        Stacked image models.
    """
    cube_data = np.array([i.data for i in container])
    cube_err = np.array([i.err for i in container])
    cube_dq = np.array([i.dq for i in container])
    cube_var_noise = np.array([i.var_rnoise for i in container])
    cube = datamodels.CubeModel(data=cube_data, err=cube_err, dq=cube_dq, var_noise=cube_var_noise)

    # update metadata of cube to match the first image
    cube.meta = container[0].meta
    return cube


def make_resamp(input_models):
    """
    Make a resampler from a set of models.

    Parameters
    ----------
    input_models : ModelLibrary
        Image models to resample.

    Returns
    -------
    ResampleImage
        The resampler with default settings for outlier_detection.
    """
    resamp = ResampleImage(
        input_models,
        output="",
        blendheaders=False,
        weight_type="ivm",
        pixfrac=1.0,
        kernel="square",
        fillval="INDEF",
        good_bits="~DO_NOT_USE",
        asn_id="test",
        enable_var=False,
        enable_ctx=False,
        compute_err="driz_err",
    )
    return resamp
