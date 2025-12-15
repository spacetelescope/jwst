import numpy as np
import pytest
from stdatamodels.jwst.datamodels import RampModel


@pytest.fixture()
def create_miri_model():
    """
    Make MIRI Ramp model for testing.

    Returns
    -------
    callable
        A function with optional parameters nints=2, ngroups=10,
        ysize=128, xsize=136. When called, it creates a RampModel
        with MIRI metadata.
    """

    def _create_model(nints=2, ngroups=10, ysize=128, xsize=136):
        # create the data and groupdq arrays
        csize = (nints, ngroups, ysize, xsize)
        data = np.full(csize, 1.0, dtype=np.float32)
        groupdq = np.zeros(csize, dtype=np.uint8)

        # create a JWST datamodel for MIRI data
        dm_ramp = RampModel(data=data, groupdq=groupdq)

        dm_ramp.meta.instrument.name = "MIRI"
        dm_ramp.meta.instrument.detector = "MIRIMAGE"
        dm_ramp.meta.observation.date = "2024-01-01"
        dm_ramp.meta.observation.time = "00:00:00"
        dm_ramp.meta.subarray.name = "SUB128"
        dm_ramp.meta.subarray.xstart = 1
        dm_ramp.meta.subarray.xsize = xsize
        dm_ramp.meta.subarray.ystart = 1
        dm_ramp.meta.subarray.ysize = ysize
        dm_ramp.meta.exposure.readpatt = "FASTR1"
        dm_ramp.meta.exposure.ngroups = ngroups
        dm_ramp.meta.exposure.nints = nints

        return dm_ramp

    return _create_model
