"""Unit tests for superbias subtraction."""

import numpy as np
from stdatamodels.jwst.datamodels import dqflags

from jwst.superbias.bias_sub import do_correction


def test_basic_superbias_subtraction(setup_full_cube):
    """Check basic superbias subtraction."""

    # Create inputs, data, and superbiases
    ngroups = 5
    nrows = 2048
    ncols = 2048
    blevel = 2000.0

    data, bias = setup_full_cube(ngroups, nrows, ncols)

    # Add signal values and bias values
    data.data[:] = blevel
    bias.data[:] = blevel

    # Run the pipeline on a copy as this is now being done in superbias_step
    output = do_correction(data.copy(), bias)

    # Make sure that manual superbias subtraction matches pipeline output
    manual = data.data - bias.data
    assert np.array_equal(output.data, manual)


def test_subarray_correction(setup_subarray_cube):
    """Check that the proper subarray is extracted from the full frame
    reference file during subtraction."""

    # Create inputs, subarray SUB320A335R data, and superbiases
    ngroups = 5
    nrows = 320
    ncols = 320
    blevel = 2000.0
    xstart = 486
    ystart = 1508

    data, bias = setup_subarray_cube(xstart, ystart, ngroups, nrows, ncols)
    manualbias = bias.data[ystart - 1 : ystart - 1 + nrows, xstart - 1 : xstart - 1 + ncols]

    # Add signal values and bias values
    data.data[:] = blevel
    manualbias[:] = blevel

    # Run the pipeline on a copy as this is now being done in superbias_step
    output = do_correction(data.copy(), bias)

    # Make sure the subarray was extracted correctly from superbias file
    # Make sure that manual superbias subtraction matches pipeline output
    manualout = data.data - manualbias
    assert np.array_equal(output.data, manualout)


def test_dq_propagation(setup_full_cube):
    """Check that the PIXELDQ array of the science exposure is correctly
    combined with the reference file DQ array."""

    # Create inputs, data, and superbiases
    ngroups = 5
    nrows = 2048
    ncols = 2048
    dqval1 = 5
    dqval2 = 10

    data, bias = setup_full_cube(ngroups, nrows, ncols)

    # Add DQ values to the data and reference file
    data.pixeldq[5, 5] = dqval1
    bias.dq[5, 5] = dqval2

    # Run the pipeline
    output = do_correction(data, bias)

    # Make sure DQ values from data and reference file are added in the output
    assert output.pixeldq[5, 5] == dqval1 + dqval2


def test_nans_in_superbias(setup_full_cube):
    """Check that no superbias subtraction is done for pixels that have a
    value of NaN in the reference file."""

    # Create inputs, data, and superbiases
    ngroups = 5
    nrows = 2048
    ncols = 2048
    blevel = 2000.0

    data, bias = setup_full_cube(ngroups, nrows, ncols)

    # Add signal values and bias values
    data.data[:] = blevel
    bias.data[:] = blevel

    # Set superbias value for pixel to NaN
    bias.data[500, 500] = np.nan

    # Run the pipeline on a copy as this is now being done in superbias_step
    output = do_correction(data.copy(), bias)

    # Check that subtraction was not done on pixel with NaN reference value
    # Check that subtraction was done on other pixels
    assert np.array_equal(output.data[0, :, 500, 500], data.data[0, :, 500, 500])
    assert np.array_equal(output.data[0, :, 50, 50], data.data[0, :, 50, 50] - blevel)


def test_zeroframe(setup_full_cube):
    """ """
    darr1 = [11800.0, 11793.0, 11823.0, 11789.0, 11857.0]
    darr2 = [11800.0, 11793.0, 11823.0, 11789.0, 11857.0]
    darr3 = [10579.0, 10594.0, 10620.0, 10583.0, 10621.0]
    zarr = [0.0, 10500.0, 10579.0]

    ngroups, nrows, ncols = len(darr1), 1, len(zarr)
    ramp, bias = setup_full_cube(ngroups, nrows, ncols)

    cdq = np.array([dqflags.group["SATURATED"]] * ngroups)

    # Set up data array
    ramp.data[0, :, 0, 0] = np.array(darr1)
    ramp.data[0, :, 0, 1] = np.array(darr2)
    ramp.data[0, :, 0, 2] = np.array(darr3)

    ramp.groupdq[0, :, 0, 0] = np.array(cdq)
    ramp.groupdq[0, :, 0, 1] = np.array(cdq)

    nints = ramp.data.shape[0]
    ramp.zeroframe = np.zeros((nints, nrows, ncols), dtype=float)

    ramp.zeroframe[0, 0, :] = np.array(zarr)

    bval = 5000.0
    bias.data[:, :] = bval

    ramp.meta.exposure.zero_frame = True
    output = do_correction(ramp, bias)

    check = np.array([zarr[0], zarr[1] - bval, zarr[2] - bval])
    np.testing.assert_equal(check, output.zeroframe[0, 0, :])
