import numpy as np
import pytest
from astropy.io import fits
from stdatamodels.jwst import datamodels


@pytest.fixture(scope="module")
def create_sci_model():
    """
    Science datamodel for testing.

    Returns
    -------
    sci_datamodel : RampModel DataModel
        Sample science datamodel for testing
    """

    def _dm(nints=3, ngroups=12, nrows=512, ncols=512, xstart=257, ystart=769):
        scidata = np.ones((nints, ngroups, nrows, ncols), dtype=np.float32)
        dm = datamodels.RampModel(data=scidata)
        dm.meta.subarray.xstart = xstart
        dm.meta.subarray.ystart = ystart
        dm.meta.exposure.frame_time = 10.73677
        dm.meta.exposure.group_time = 21.47354
        dm.meta.exposure.ngroups = ngroups
        dm.meta.exposure.nframes = 1
        dm.meta.exposure.groupgap = 1
        dm.meta.exposure.nresets_at_start = 0
        dm.meta.exposure.nresets_between_ints = 1
        dm.meta.instrument.detector = "NRCA1"
        dm.meta.instrument.name = "NIRCAM"
        dm.meta.exposure.start_time = 59672.47338658912
        dm.meta.observation.date = "2022-04-03"
        dm.meta.observation.time = "11:21:40"
        dm.meta.filename = "test.fits"

        return dm

    return _dm


@pytest.fixture(scope="module")
def create_traps_filled_model():
    """
    Trapsfilled model for testing.

    Returns
    -------
    trapsfilled_model : TrapsFilledModel DataModel
        Sample TrapsFilledModel for testing
    """

    def _dm(nrows=2048, ncols=2048):
        tfdata = np.ones((3, nrows, ncols), dtype=np.float32)
        tf_dm = datamodels.TrapsFilledModel(data=tfdata)
        tf_dm.meta.subarray.xstart = 1
        tf_dm.meta.subarray.ystart = 1
        tf_dm.meta.exposure.end_time = 59672.47462927084
        return tf_dm

    return _dm


@pytest.fixture(scope="module")
def create_trap_density_model():
    """
    Trap density model for testing.

    Returns
    -------
    trapdensity_model : TrapDensityModel DataModel
        Sample TrapDensityModel for testing
    """

    def _dm(nrows=2048, ncols=2048):
        tddata = np.ones((nrows, ncols), dtype=np.float32)
        td_dm = datamodels.TrapDensityModel(data=tddata)
        td_dm.meta.subarray.xstart = 1
        td_dm.meta.subarray.ystart = 1
        return td_dm

    return _dm


@pytest.fixture()
def create_trappars_model():
    """
    Create TrapPars datamodel for testing.

    Returns
    -------
    trappars_model : TrapParsModel DataModel
        Sample TrapParsModel for testing
    """
    capture0 = fits.Column(name="capture0", array=np.array([180.0, 270.0, 80.0]), format="D")
    capture1 = fits.Column(name="capture1", array=np.array([-0.0004, -0.004, -0.0009]), format="D")
    capture2 = fits.Column(name="capture2", array=np.array([290.0, 140.0, 320.0]), format="D")
    decay = fits.Column(name="decay_param", array=np.array([-0.01, -0.001, -0.0002]), format="D")
    column_definitions = fits.ColDefs([capture0, capture1, capture2, decay])
    trappars_table = fits.BinTableHDU.from_columns(column_definitions, nrows=3)
    trappars_model = datamodels.TrapParsModel()
    trappars_model.trappars_table = trappars_table.data
    return trappars_model
