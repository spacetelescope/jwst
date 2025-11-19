import numpy as np
from stdatamodels.jwst import datamodels

from jwst.pipeline import Detector1Pipeline


def make_miri_rampmodel(nints=1, ngroups=5, ysize=1024, xsize=1032):
    """
    Make a MIRI image ramp model with minimal metadata.

    Parameters
    ----------
    nints : int, optional
        Number of integrations.
    ngroups : int, optional
        Number of groups.
    ysize : int, optional
        Y size.
    xsize : int, optional
        X size.

    Returns
    -------
    ramp : `stdatamodels.jwst.datamodels.RampModel`
        A ramp model with only the data array and minimal metadata set.
    """
    csize = (nints, ngroups, ysize, xsize)
    data = np.full(csize, 1.0)

    ramp = datamodels.RampModel(data=data)
    ramp.meta.filename = "test_miri.fits"
    ramp.meta.instrument.name = "MIRI"
    ramp.meta.instrument.detector = "MIRIMAGE"
    ramp.meta.exposure.frame_time = 1.0
    ramp.meta.exposure.groupgap = 0
    ramp.meta.exposure.group_time = 1.0
    ramp.meta.exposure.nframes = 1
    ramp.meta.exposure.nints = nints
    ramp.meta.exposure.ngroups = ngroups
    ramp.meta.exposure.readpatt = "FASTR1"
    ramp.meta.observation.date = "2024-01-01"
    ramp.meta.observation.time = "00:00:00"
    ramp.meta.subarray.name = "FULL"
    ramp.meta.subarray.xstart = 1
    ramp.meta.subarray.xsize = xsize
    ramp.meta.subarray.ystart = 1
    ramp.meta.subarray.ysize = ysize

    return ramp


def test_detector1_miri(tmp_path):
    input_model = make_miri_rampmodel()
    input_model_copy = input_model.copy()
    Detector1Pipeline.call(input_model, save_results=True, output_dir=str(tmp_path))

    # Check for expected output
    assert (tmp_path / "test_miri_rate.fits").exists()
    assert (tmp_path / "test_miri_rateints.fits").exists()

    # Input is not modified
    np.testing.assert_allclose(input_model.data, input_model_copy.data)
    assert input_model.meta.cal_step._instance == input_model_copy.meta.cal_step._instance
