import pytest
from stdatamodels.jwst.datamodels import RampModel, SuperBiasModel


@pytest.fixture(scope="function")
def setup_full_cube():
    """
    Set up mock NIRCam FULL data to test.

    Returns
    -------
    callable
        A function with arguments nroups, nrows, ncols.
        When called, it creates a RampModel with NIRCam metadata.
    """

    def _cube(ngroups, nrows, ncols):
        nints = 1

        # create a JWST datamodel for NIRCam FULL data
        data_model = RampModel((nints, ngroups, nrows, ncols))
        data_model.meta.subarray.xstart = 1
        data_model.meta.subarray.ystart = 1
        data_model.meta.subarray.xsize = ncols
        data_model.meta.subarray.ysize = nrows
        data_model.meta.exposure.ngroups = ngroups
        data_model.meta.instrument.name = "NIRCAM"
        data_model.meta.instrument.detector = "NRCA1"
        data_model.meta.observation.date = "2017-10-01"
        data_model.meta.observation.time = "00:00:00"

        # create a superbias model for the superbias step
        bias_model = SuperBiasModel((2048, 2048))
        bias_model.meta.subarray.xstart = 1
        bias_model.meta.subarray.ystart = 1
        bias_model.meta.subarray.xsize = 2048
        bias_model.meta.subarray.ysize = 2048
        bias_model.meta.instrument.name = "NIRCAM"
        bias_model.meta.description = "Fake data."
        bias_model.meta.telescope = "JWST"
        bias_model.meta.reftype = "SuperBiasModel"
        bias_model.meta.author = "Alicia"
        bias_model.meta.pedigree = "Dummy"
        bias_model.meta.useafter = "2015-10-01T00:00:00"

        return data_model, bias_model

    return _cube


@pytest.fixture(scope="function")
def setup_subarray_cube():
    """
    Set up mock NIRCam subarray data to test.

    Returns
    -------
    callable
        A function with arguments nroups, nrows, ncols.
        When called, it creates a RampModel with NIRCam metadata.
    """

    def _cube(xstart, ystart, ngroups, nrows, ncols):
        nints = 1

        # create a JWST datamodel for NIRCam SUB320A335R data
        data_model = RampModel((nints, ngroups, nrows, ncols))
        data_model.meta.subarray.name = "SUB320A335R"
        data_model.meta.subarray.xstart = xstart
        data_model.meta.subarray.ystart = ystart
        data_model.meta.subarray.xsize = ncols
        data_model.meta.subarray.ysize = nrows
        data_model.meta.exposure.ngroups = ngroups
        data_model.meta.instrument.name = "NIRCAM"
        data_model.meta.instrument.detector = "NRCALONG"
        data_model.meta.observation.date = "2019-10-14"
        data_model.meta.observation.time = "16:44:12.000"

        # create a superbias model for the superbias step
        bias_model = SuperBiasModel((2048, 2048))
        bias_model.meta.subarray.xstart = 1
        bias_model.meta.subarray.ystart = 1
        bias_model.meta.subarray.xsize = 2048
        bias_model.meta.subarray.ysize = 2048
        bias_model.meta.instrument.name = "NIRCAM"
        bias_model.meta.description = "Fake data."
        bias_model.meta.telescope = "JWST"
        bias_model.meta.reftype = "SuperBiasModel"
        bias_model.meta.author = "Alicia"
        bias_model.meta.pedigree = "Dummy"
        bias_model.meta.useafter = "2015-10-01T00:00:00"

        return data_model, bias_model

    return _cube
