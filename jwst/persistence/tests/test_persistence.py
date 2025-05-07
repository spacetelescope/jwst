"""Test the routines in persistence.py."""

import numpy as np
from astropy.io import fits
import pytest

from jwst import datamodels

from jwst.persistence import persistence


def test_no_nan():
    """Test the no_nan function."""
    input_model = datamodels.RampModel()
    input_model.data = 2.0 * np.ones((5, 10, 10, 10))
    # Datamodel with no nans or zeros should be unchanged
    result = persistence.no_nan(input_model, 1.0, zap_nan=True, zap_zero=True)
    np.testing.assert_equal(result.data, input_model.data)

    input_model.data[0] = np.nan
    input_model.data[1] = 0.0
    # Test replacing both NaNs and zeros
    result = persistence.no_nan(input_model, 12.0, zap_nan=True, zap_zero=True)
    np.testing.assert_equal(result.data[0], 12.0)
    np.testing.assert_equal(result.data[1], 12.0)
    np.testing.assert_equal(result.data[2:], 2.0)
    # Test replacing only NaNs
    result = persistence.no_nan(input_model, 11.0, zap_nan=True)
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
        dm.meta.exposure.start_time = 59672.47338658912
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

    def _dm(nrows, ncols):
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


def test_get_slice(create_sci_model, create_traps_filled_model):
    """Test the get_slice method."""
    # Create empty DataSet
    ds = persistence.DataSet(None, None, 40.0, False, None, None, None)
    sci = create_sci_model(3, 12, 512, 512, 257, 769)
    tf = create_traps_filled_model(2048, 2048)
    returned_slice = ds.get_slice(tf, sci)
    assert returned_slice == (slice(768, 1280, None), slice(256, 768, None))


def test_ref_matches_sci(create_traps_filled_model):
    """
    Test the ref_matches_sci method.

    Parameters
    ----------
    create_traps_filled_model : pytest fixture
        Fixture that returns a TrapsFilledModel for testing
    """
    # Create empty DataSet
    ds = persistence.DataSet(None, None, 40.0, False, None, None, None)
    tf = create_traps_filled_model(2048, 2048)
    test_slices = [
        (slice(0, 2048, None), slice(0, 2048, None)),
        (slice(256, 768, None), slice(1024, 1536, None)),
    ]
    expected_results = [True, False]
    for this_slice, expected in zip(test_slices, expected_results, strict=True):
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
    """
    Fixture to create arrays of trappar parameters.

    Returns
    -------
    pars : tuple of 4 numpy.ndarrays
        The tuple of parameter arrays
    """
    par0 = np.array([180.0, 270.0, 80.0])
    par1 = np.array([-0.0004, -0.004, -0.0009])
    par2 = np.array([290.0, 140.0, 320.0])
    par3 = np.array([-0.01, -0.001, -0.0002])
    pars = (par0, par1, par2, par3)
    return pars


class TrapParsTester:
    """Class to handle testing trap pars."""

    def test_get_parameters(self):
        """Test the get_parameters method."""
        # Create empty DataSet
        ds = persistence.DataSet(None, None, 40.0, False, None, None, None)
        ds.trappars_model = create_trappars_model()
        pars = ds.get_parameters()
        assert np.allclose(pars[0], np.array([180.0, 270.0, 80.0]), rtol=1.0e-7, atol=1.0e-7)
        assert np.allclose(pars[1], np.array([-0.0004, -0.004, -0.0009]), rtol=1.0e-7, atol=1.0e-7)
        assert np.allclose(pars[2], np.array([290.0, 140.0, 320.0]), rtol=1.0e-7, atol=1.0e-7)
        assert np.allclose(pars[3], np.array([-0.01, -0.001, -0.0002]), rtol=1.0e-7, atol=1.0e-7)

    def test_get_capture_param(self, generate_trap_pars):
        """
        Test the get_capture_param method.

        Parameters
        ----------
        generate_trap_pars : pytest fixture
            Fixture that returns trappars parameter arrays for testing
        """
        ds = persistence.DataSet(None, None, 40.0, False, None, None, None)
        pars = generate_trap_pars()
        index = 1
        returned_pars = ds.get_capture_param(pars, index)
        assert np.allclose(returned_pars[0], 270.0, rtol=1.0e-7, atol=1.0e-7)
        assert np.allclose(returned_pars[1], -0.004, rtol=1.0e-7, atol=1.0e-7)
        assert np.allclose(returned_pars[2], 140.0, rtol=1.0e-7, atol=1.0e-7)

    def test_get_decay_param(self, generate_trap_pars):
        """
        Test the get_decay_param method.

        Parameters
        ----------
        generate_trap_pars : pytest fixture
            Fixture that returns TrapPars parameter arrays for testing
        """
        ds = persistence.DataSet(None, None, 40.0, False, None, None, None)
        pars = generate_trap_pars()
        index = 1
        capture_param = ds.get_capture_param(pars, index)
        assert np.allclose(capture_param, -0.001, rtol=1.0e-7, atol=1.0e-7)


def test_get_group_info(create_sci_model):
    """
    Test the get_group_info method.

    Parameters
    ----------
    create_sci_model : pytest fixture
        Fixture that returns a science DataModel for testing
    """
    output_model = create_sci_model(3, 12, 512, 512, 257, 769)
    ds = persistence.DataSet(output_model, None, 40.0, False, None, None, None)
    integration = 1
    ds.get_group_info(integration)
    assert np.allclose(ds.tframe, 10.73677, rtol=1.0e-7, atol=1.0e-7)
    assert np.allclose(ds.tgroup, 21.47354, rtol=1.0e-7, atol=1.0e-7)
    assert ds.ngroups == 12
    assert ds.nframes == 1
    assert ds.groupgap == 1
    assert ds.nresets == 1


def create_persistencesat_model():
    """
    Create a PersistencesatModel for testing.

    Returns
    -------
    persat_model : DataModel
        The minimal PersistencesatModel for testing
    """
    scidata = np.ones((2048, 2048), dtype=np.float32)
    scidata = scidata[:] * 40.0
    persat_model = datamodels.PersistenceSatModel(data=scidata)
    persat_model.meta.subarray.xstart = 1
    persat_model.meta.subarray.ystart = 1
    return persat_model


def test_do_all(create_sci_model, create_traps_filled_model, create_trap_density_model):
    """
    Test the do_all method.

    Parameters
    ----------
    create_sci_model : pytest fixture
        Fixture that returns a science DataModel for testing
    create_traps_filled_model : pytest fixture
        Fixture that returns a TrapsFilledModel for testing
    create_trap_density_model : pytest fixture
        Fixture that returns a TrapDensityModel for testing
    """
    input_model = create_sci_model(3, 12, 512, 512, 257, 769)
    traps_filled_model = create_traps_filled_model(2048, 2048)
    trappars_model = create_trappars_model()
    trap_density_model = create_trap_density_model(2048, 2048)
    persistencesat_model = create_persistencesat_model()
    ds = persistence.DataSet(
        input_model,
        traps_filled_model,
        40.0,
        False,
        trap_density_model,
        trappars_model,
        persistencesat_model,
    )
    # Check that the do_all method runs to completion
    output_model, output_trapsfilled_model, output_persistence, skipped = ds.do_all()
    # Check a couple of the pixel values returned by the method
    assert np.allclose(output_model.data[0, -1, 100, 100], -2.0079181, rtol=1.0e-7, atol=1.0e-7)
    assert np.allclose(
        output_trapsfilled_model.data[0, 100, 100], 0.0010368865, rtol=1.0e-7, atol=1.0e-7
    )
