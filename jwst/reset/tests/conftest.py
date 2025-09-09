import numpy as np
import pytest
from stdatamodels.jwst.datamodels import RampModel, ResetModel


@pytest.fixture(scope="function")
def make_rampmodel():
    """
    Make MIRI Ramp model for testing.

    Returns
    -------
    callable
        A function with optional parameters nints=1, ngroups=15,
        ysize=128, xsize=136. When called, it creates a RampModel
        with MIRI metadata.
    """

    def _ramp(nints=1, ngroups=15, ysize=128, xsize=136):
        # create the data and groupdq arrays
        csize = (nints, ngroups, ysize, xsize)
        data = np.full(csize, 1.0)  # default = 1.0

        # create a JWST datamodel for MIRI data
        dm_ramp = RampModel(data=data)

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
        dm_ramp.meta.exposure.nints = nints
        dm_ramp.meta.exposure.ngroups = ngroups
        dm_ramp.meta.description = "Fake data."

        return dm_ramp

    return _ramp


@pytest.fixture(scope="function")
def make_resetmodel():
    """
    Make MIRI Reset model for testing.

    Returns
    -------
    callable
        A function with optional parameters ngroups=15, ysize=128, xsize=136.
        When called, it creates a ResetModel with MIRI metadata.
    """

    def _reset(ngroups=15, ysize=128, xsize=136):
        # create the data and groupdq arrays
        nints = 2
        csize = (nints, ngroups, ysize, xsize)
        data = np.full(csize, 1.0)  # default = 1.0

        # create a JWST datamodel for MIRI data
        reset = ResetModel(data=data)
        reset.meta.exposure.nints = nints
        reset.meta.exposure.ngroups = ngroups
        reset.meta.instrument.name = "MIRI"
        reset.meta.date = "2018-01-01"
        reset.meta.time = "00:00:00"
        reset.meta.description = "Fake data."
        reset.meta.reftype = "ResetModel"
        reset.meta.author = "Jane Morrison"
        reset.meta.pedigree = "Dummy"
        reset.meta.useafter = "2015-10-01T00:00:00"
        reset.meta.instrument.detector = "MIRIMAGE"
        reset.meta.subarray.name = "SUB128"
        return reset

    return _reset
