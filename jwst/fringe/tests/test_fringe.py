"""Unit tests for fringe correction."""

import pytest
import numpy as np
from numpy.testing import assert_array_equal

from stdatamodels.jwst.datamodels import IFUImageModel, FringeModel

from jwst.fringe import fringe


def test_data_correction(setup_inputs):
    """Test both good and NaN pixels."""

    FRINGE_CONSTANT = 2.0  # correction will be input data divided by this factor
    FRINGE_CONSTANT_SQ = 4.0

    input_model, fringe_model = setup_inputs(
        shape=(4, 5), bad_corner=True, fringe_constant=FRINGE_CONSTANT
    )

    # Do the correction.
    output_model = fringe.do_correction(input_model, fringe_model)

    # Check that correction was done on pixels with valid values for all
    # SCI, ERR, and VAR arrays.
    good_pix = np.where(np.isfinite(input_model.data))

    assert_array_equal(output_model.data[good_pix], (input_model.data * FRINGE_CONSTANT)[good_pix])
    assert_array_equal(output_model.err[good_pix], (input_model.err * FRINGE_CONSTANT)[good_pix])
    assert_array_equal(
        output_model.var_poisson[good_pix], (input_model.var_poisson * FRINGE_CONSTANT_SQ)[good_pix]
    )
    assert_array_equal(
        output_model.var_rnoise[good_pix], (input_model.var_rnoise * FRINGE_CONSTANT_SQ)[good_pix]
    )
    assert_array_equal(
        output_model.var_flat[good_pix], (input_model.var_flat * FRINGE_CONSTANT_SQ)[good_pix]
    )

    # Check that correction was not done on pixel with NaN values for all SCI,
    # ERR, and VAR arrays (i.e., these pixels have not been corrected).
    assert np.isnan(output_model.data[0, 0])
    assert np.isnan(output_model.err[0, 0])
    assert np.isnan(output_model.var_poisson[0, 0])
    assert np.isnan(output_model.var_rnoise[0, 0])
    assert np.isnan(output_model.var_flat[0, 0])


@pytest.fixture
def setup_inputs():
    """Create input and fringe models."""

    def _setup(shape=(2, 2), bad_corner=False, fringe_constant=1):
        input_data = np.ones(shape) * 6.0
        input_err = input_data * 0.1
        err_poisson = np.sqrt(input_data)
        fake_var = input_err * input_err * 0.3333
        input_model = IFUImageModel(
            data=input_data,
            err=input_err,
            var_poisson=err_poisson * err_poisson,
            var_rnoise=fake_var,
            var_flat=fake_var,
        )

        fringe_model = FringeModel(data=np.ones(shape) / fringe_constant)

        if bad_corner:
            # Make 1 bad pixel.
            input_model.data[0, 0] = np.nan
            input_model.err[0, 0] = np.nan
            input_model.var_poisson[0, 0] = np.nan
            input_model.var_rnoise[0, 0] = np.nan
            input_model.var_flat[0, 0] = np.nan

        return input_model, fringe_model

    return _setup
