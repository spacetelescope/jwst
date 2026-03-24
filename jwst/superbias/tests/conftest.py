import numpy as np
import pytest
from stdatamodels.jwst.datamodels import RampModel, SuperBiasModel

from jwst.dq_init.tests.helpers import make_superstripe_model


def add_test_refmodel_metadata(refmodel):
    """
    Add some basic required metadata to a reference file model.

    TODO: consolidate this function with similar ones elsewhere.
    """
    refmodel.meta.telescope = "JWST"
    refmodel.meta.description = "filler"
    refmodel.meta.reftype = "filler"
    refmodel.meta.author = "Py Test"
    refmodel.meta.pedigree = "Pytest"
    refmodel.meta.useafter = "2015-01-01T01:00:00"


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
        add_test_refmodel_metadata(bias_model)

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
        add_test_refmodel_metadata(bias_model)

        return data_model, bias_model

    return _cube


@pytest.fixture(scope="function")
def setup_nis_superstripe_cube():
    """
    Set up mock NIRISS superstripe data to test.

    TODO: this fixture is copied from saturation tests. The two copies should be consolidated.

    Returns
    -------
    callable
        A function that creates a RampModel and a SuperBiasModel when called.
    """

    def _cube():
        data_model = make_superstripe_model()
        num_stripes = data_model.meta.subarray.num_superstripe
        data_model.pixeldq = np.zeros((num_stripes, *data_model.data.shape[-2:]), dtype=np.uint32)

        # create a bias model
        bias_model = SuperBiasModel((2028, 2048))
        bias_model.data = np.ones((2048, 2048)) * 15000  # bias for every pixel is 15000
        bias_model.meta.instrument.name = "NIRISS"
        bias_model.meta.instrument.detector = "NIS"
        bias_model.meta.subarray.xstart = 1
        bias_model.meta.subarray.xsize = 2048
        bias_model.meta.subarray.ystart = 1793
        bias_model.meta.subarray.ysize = 256
        bias_model.meta.subarray.fastaxis = -2
        bias_model.meta.subarray.slowaxis = -1
        bias_model.meta.exposure.readpatt = "ANY"
        add_test_refmodel_metadata(bias_model)

        return data_model, bias_model

    return _cube
