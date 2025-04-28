"""Test the routines in persistence.py."""

import numpy as np
from astropy.io import fits
import pytest

from jwst import datamodels

from jwst.persistence import persistence


def test_no_NaN():
    input_model = datamodels.RampModel()
    input_model.data = 2.0 * np.ones((5, 10, 10, 10))
    # Datamodel with no nans or zeros should be unchanged
    result = persistence.no_NaN(input_model, 1.0, zap_nan=True, zap_zero=True)
    np.testing.assert_equal(result.data, input_model.data)

    input_model.data[0] = np.nan
    input_model.data[1] = 0.0
    # Test replacing both NaNs and zeros
    result = persistence.no_NaN(input_model, 12.0, zap_nan=True, zap_zero=True)
    np.testing.assert_equal(result.data[0], 12.0)
    np.testing.assert_equal(result.data[1], 12.0)
    np.testing.assert_equal(result.data[2:], 2.0)
    # Test replacing only NaNs
    result = persistence.no_NaN(input_model, 11.0, zap_nan=True)
    np.testing.assert_equal(result.data[0], 11.0)
    np.testing.assert_equal(result.data[1], 0.0)
    np.testing.assert_equal(result.data[2:], 2.0)


@pytest.fixture(scope="module")
def create_sci_model():
    """
    Science datamodel for testing.

    Returns
    -------
    sci_datamodel : RampModel DataModel
        Sample science datamodel for testing
    """

    def _dm(nints, ngroups, nrows, ncols, xstart, ystart):
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
        return dm

    return _dm


@pytest.fixture(scope="module")
def create_traps_filled_model():
    """
    Trapsfilled model for testing.

    Returns
    -------
    trapsfilled_model : TrapsFilledModel DataModel
        Sample TrapsFilledModel datamodel for testing
    """

    def _dm(nrows, ncols):
        tfdata = np.ones((3, nrows, ncols), dtype=np.float32)
        tf_dm = datamodels.TrapsFilledModel(data=tfdata)
        tf_dm.meta.subarray.xstart = 1
        tf_dm.meta.subarray.ystart = 1
        return tf_dm

    return _dm


@pytest.fixture(scope="module")
def create_trap_density_model():
    """
    Trap density model for testing.

    Returns
    -------
    trapdensity_model : TrapDensityModel DataModel
        Sample TrapDensityModel datamodel for testing
    """

    def _dm(nrows, ncols):
        tddata = np.ones((nrows, ncols), dtype=np.float32)
        td_dm = datamodels.TrapDensityModel(data=tddata)
        td_dm.meta.subarray.xstart = 1
        td_dm.meta.subarray.ystart = 1
        return td_dm

    return _dm


def create_trappars_model():
    """
    Create TrapPars datamodel for testing.

    Returns
    -------
    trappars_model : TrapParsModel DataModel
        Sample TrapParsModel datamodel for testing
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


def test_get_slice(create_sci_model, create_traps_filled_model):
    """Test the get_slice method."""
    # Create empty DataSet
    ds = persistence.DataSet(None, None, 40.0, False, None, None, None)
    sci = create_sci_model(3, 12, 512, 512, 257, 769)
    tf = create_traps_filled_model(2048, 2048)
    returned_slice = ds.get_slice(tf, sci)
    assert returned_slice == (slice(768, 1280, None), slice(256, 768, None))


def test_ref_matches_sci(create_traps_filled_model):
    # Create empty DataSet
    ds = persistence.DataSet(None, None, 40.0, False, None, None, None)
    tf = create_traps_filled_model(2048, 2048)
    test_slices = [
        (slice(0, 2048, None), slice(0, 2048, None)),
        (slice(256, 768, None), slice(1024, 1536, None)),
    ]
    expected_results = [True, False]
    for this_slice, expected in zip(test_slices, expected_results):
        assert ds.ref_matches_sci(tf, this_slice) == expected


def test_get_subarray(create_trap_density_model):
    """
    Test the get_subarray method.

    This method is only run with 2-d reference file datamodels
    (trapdensity and persestencesat)
    """
    # Create empty DataSet
    ds = persistence.DataSet(None, None, 40.0, False, None, None, None)
    td = create_trap_density_model(2048, 2048)
    td.data[256:768, 1024:1536] = 12.0
    subarray = ds.get_subarray(td, (slice(256, 768, None), slice(1024, 1536, None)))
    np.testing.assert_equal(subarray.data, 12.0)


@pytest.fixture(scope="class")
def generate_trap_pars():
    par0 = np.array([180.0, 270.0, 80.0])
    par1 = np.array([-0.0004, -0.004, -0.0009])
    par2 = np.array([290.0, 140.0, 320.0])
    par3 = np.array([-0.01, -0.001, -0.0002])
    pars = (par0, par1, par2, par3)
    return pars


class TrapParsTester:
    """Class to handle testing trap pars."""

    def test_get_parameters(self):
        # Create empty DataSet
        ds = persistence.DataSet(None, None, 40.0, False, None, None, None)
        ds.trappars_model = create_trappars_model()
        pars = ds.get_parameters()
        np.testing.assert_equal(pars[0], np.array([180.0, 270.0, 80.0]))
        np.testing.assert_equal(pars[1], np.array([-0.0004, -0.004, -0.0009]))
        np.testing.assert_equal(pars[2], np.array([290.0, 140.0, 320.0]))
        np.testing.assert_equal(pars[3], np.array([-0.01, -0.001, -0.0002]))

    def test_get_capture_param(self, generate_trap_pars):
        """
        Test the get_capture_param method.

        Input data is a tuple of 4 float arrays of length 3
        But still need to create an empty dataset
        """
        ds = persistence.DataSet(None, None, 40.0, False, None, None, None)
        pars = generate_trap_pars()
        index = 1
        returned_pars = ds.get_capture_param(pars, index)
        assert returned_pars[0] == 270.0
        assert returned_pars[1] == -0.004
        assert returned_pars[2] == 140.0

    def test_get_decay_param(self, generate_trap_pars):
        ds = persistence.DataSet(None, None, 40.0, False, None, None, None)
        pars = generate_trap_pars()
        index = 1
        capture_param = ds.get_capture_param(pars, index)
        assert capture_param == -0.001


def test_get_group_info(create_sci_model):
    output_model = create_sci_model(3, 12, 512, 512, 257, 769)
    ds = persistence.DataSet(output_model, None, 40.0, False, None, None, None)
    integration = 1
    ds.get_group_info(integration)
    assert ds.tframe == 10.73677
    assert ds.tgroup == 21.47354
    assert ds.ngroups == 12
    assert ds.nframes == 1
    assert ds.groupgap == 1
    assert ds.nresets == 1
