import numpy as np
from stdatamodels.jwst.datamodels import RampModel, ReadnoiseModel

__all__ = ["niriss_ramp_model", "niriss_ramp_model_with_signal", "mock_readnoise_model"]


def niriss_ramp_model(shape=(2, 10, 5, 5)):
    """
    Create a mock NIRISS ramp model.

    Arrays have the specified dimensions for the data, pixel DQ, group DQ, and
    ERR extensions. All groups in all arrays have initial values set to 0.

    Parameters
    ----------
    shape : tuple of int, optional
        Shape for the data arrays.

    Returns
    -------
    ramp_model : `~stdatamodels.jwst.datamodels.RampModel`
        The mock ramp model.
    """
    nints, ngroups, nrows, ncols = shape
    err = np.zeros(shape, dtype=np.float32)
    data = np.zeros(shape, dtype=np.float32)
    pixdq = np.zeros((nrows, ncols), dtype=np.uint32)
    gdq = np.zeros(shape, dtype=np.uint8)

    # Create and populate ramp model
    ramp_model = RampModel(data=data, err=err, pixeldq=pixdq, groupdq=gdq)
    ramp_model.meta.instrument.name = "NIRISS"
    ramp_model.meta.instrument.detector = "NIS"

    ramp_model.meta.observation.date = "2025-10-13"
    ramp_model.meta.observation.time = "00:00:00"

    ramp_model.meta.subarray.name = "FULL"
    ramp_model.meta.subarray.xstart = 1
    ramp_model.meta.subarray.ystart = 1
    ramp_model.meta.subarray.xsize = ncols
    ramp_model.meta.subarray.ysize = nrows

    ramp_model.meta.exposure.nints = nints
    ramp_model.meta.exposure.ngroups = ngroups
    ramp_model.meta.exposure.nframes = 1
    ramp_model.meta.exposure.groupgap = 0
    ramp_model.meta.exposure.drop_frames1 = 0
    ramp_model.meta.exposure.type = "NIS_IMAGE"

    return ramp_model


def niriss_ramp_model_with_signal(chg_thresh, pixel, shape=(1, 10, 5, 5), start_group=3):
    """
    Create a mock ramp model with flaggable signal in one pixel.

    Set up ramp array that is by default 5, 5 with 10 groups.
    The flagging starts in group 3 (zero based) in the pixel tested.

    Parameters
    ----------
    chg_thresh : float
        Signal threshold for flagging.
    pixel : tuple of int
        Pixel to be flagged.
    shape : tuple of int, optional
        Shape for the data arrays.
    start_group : int, optional
        Starting group for adding flaggable signal to.

    Returns
    -------
    ramp_model : `~stdatamodels.jwst.datamodels.RampModel`
        The mock ramp model.
    """
    nints, ngroups, nrows, ncols = shape
    ramp_model = niriss_ramp_model(shape)

    # Set up mock data
    base = chg_thresh * 0.05
    base_arr = [float(k + 1) * base for k in range(ngroups)]
    for row in range(nrows):
        for col in range(ncols):
            ramp_model.data[0, :, row, col] = np.array(base_arr, dtype=np.float32)

    # Make CHARGELOSS threshold starting at group 3
    in_row, in_col = pixel
    ramp_model.data[0, start_group:, in_row, in_col] += chg_thresh

    return ramp_model


def mock_readnoise_model(shape, value=1.0):
    """
    Mock a readnoise model for NIRISS.

    Parameters
    ----------
    shape : tuple of int
        Shape for the 2D array.
    value : float, optional
        Readnoise value to fill the array with.

    Returns
    -------
    readnoise_model : `~stdatamodels.jwst.datamodels.ReadnoiseModel`
        The mock readnoise model.
    """
    rnoise = np.full(shape, value, dtype=np.float64)
    readnoise_model = ReadnoiseModel(data=rnoise)
    readnoise_model.meta.instrument.name = "NIRISS"
    readnoise_model.meta.subarray.xstart = 1
    readnoise_model.meta.subarray.ystart = 1
    readnoise_model.meta.subarray.xsize = shape[1]
    readnoise_model.meta.subarray.ysize = shape[0]
    return readnoise_model
