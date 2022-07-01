"""

Unit tests for superbias subtraction

"""

import pytest
import numpy as np
from jwst.superbias import SuperBiasStep
from jwst.superbias.bias_sub import do_correction
from jwst.datamodels import RampModel, SuperBiasModel
from jwst.datamodels import dqflags


def test_basic_superbias_subtraction(setup_full_cube):
    '''Check basic superbias subtraction.'''

    # Create inputs, data, and superbiases
    ngroups = 5
    nrows = 2048
    ncols = 2048
    blevel = 2000.

    data, bias = setup_full_cube(ngroups, nrows, ncols)

    # Add signal values and bias values
    data.data[:] = blevel
    bias.data[:] = blevel

    # Run the pipeline
    output = do_correction(data, bias)

    # Make sure that manual superbias subtraction matches pipeline output
    manual = data.data - bias.data
    assert np.array_equal(output.data, manual)


def test_subarray_correction(setup_subarray_cube):
    '''Check that the proper subarray is extracted from the full frame
       reference file during subtraction.'''

    # Create inputs, subarray SUB320A335R data, and superbiases
    ngroups = 5
    nrows = 320
    ncols = 320
    blevel = 2000.
    xstart = 486
    ystart = 1508

    data, bias = setup_subarray_cube(xstart, ystart, ngroups, nrows, ncols)
    manualbias = bias.data[ystart - 1:ystart - 1 + nrows, xstart - 1:xstart - 1 + ncols]

    # Add signal values and bias values
    data.data[:] = blevel
    manualbias[:] = blevel

    # Run the pipeline
    output = do_correction(data, bias)

    # Make sure the subarray was extracted correctly from superbias file
    # Make sure that manual superbias subtraction matches pipeline output
    manualout = data.data - manualbias
    assert np.array_equal(output.data, manualout)


def test_dq_propagation(setup_full_cube):
    '''Check that the PIXELDQ array of the science exposure is correctly
       combined with the reference file DQ array.'''

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
    '''Check that no superbias subtraction is done for pixels that have a
       value of NaN in the reference file.'''

    # Create inputs, data, and superbiases
    ngroups = 5
    nrows = 2048
    ncols = 2048
    blevel = 2000.

    data, bias = setup_full_cube(ngroups, nrows, ncols)

    # Add signal values and bias values
    data.data[:] = blevel
    bias.data[:] = blevel

    # Set superbias value for pixel to NaN
    bias.data[500, 500] = np.nan

    # Run the pipeline
    output = do_correction(data, bias)

    # Check that subtraction was not done on pixel with NaN reference value
    # Check that subtraction was done on other pixels
    assert np.array_equal(output.data[0, :, 500, 500], data.data[0, :, 500, 500])
    assert np.array_equal(output.data[0, :, 50, 50], data.data[0, :, 50, 50] - blevel)


def test_full_step(setup_full_cube):
    '''Test full run of the SuperBiasStep.'''

    # Create inputs, data, and superbiases
    ngroups = 5
    nrows = 2048
    ncols = 2048

    data, bias = setup_full_cube(ngroups, nrows, ncols)

    # Add signal values and bias values
    # Use signal = 0 ADU so value will be negative after superbias step
    data.data[0, :, 500, 500] = 0

    # Run the pipeline
    output = SuperBiasStep.call(data)

    # Check that pixel value is negative after bias is subtracted
    assert np.sign(output.data[0, 0, 500, 500]) == -1


def test_zeroframe(setup_full_cube):
    """
    """
    darr1 = [11800., 11793., 11823., 11789., 11857.]
    darr2 = [11800., 11793., 11823., 11789., 11857.]
    darr3 = [10579., 10594., 10620., 10583., 10621.]
    zarr = [0., 10500., 10579.]

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

    bval = 5000.
    bias.data[:, :] = bval

    ramp.meta.exposure.zero_frame = True
    output = do_correction(ramp, bias)

    check = np.array([zarr[0], zarr[1] - bval, zarr[2] - bval])
    np.testing.assert_equal(check, output.zeroframe[0, 0, :])


@pytest.fixture(scope='function')
def setup_full_cube():
    ''' Set up fake NIRCam FULL data to test.'''

    def _cube(ngroups, nrows, ncols):

        nints = 1

        # create a JWST datamodel for NIRCam FULL data
        data_model = RampModel((nints, ngroups, nrows, ncols))
        data_model.meta.subarray.xstart = 1
        data_model.meta.subarray.ystart = 1
        data_model.meta.subarray.xsize = ncols
        data_model.meta.subarray.ysize = nrows
        data_model.meta.exposure.ngroups = ngroups
        data_model.meta.instrument.name = 'NIRCAM'
        data_model.meta.instrument.detector = 'NRCA1'
        data_model.meta.observation.date = '2017-10-01'
        data_model.meta.observation.time = '00:00:00'

        # create a superbias model for the superbias step
        bias_model = SuperBiasModel((2048, 2048))
        bias_model.meta.subarray.xstart = 1
        bias_model.meta.subarray.ystart = 1
        bias_model.meta.subarray.xsize = 2048
        bias_model.meta.subarray.ysize = 2048
        bias_model.meta.instrument.name = 'NIRCAM'
        bias_model.meta.description = 'Fake data.'
        bias_model.meta.telescope = 'JWST'
        bias_model.meta.reftype = 'SuperBiasModel'
        bias_model.meta.author = 'Alicia'
        bias_model.meta.pedigree = 'Dummy'
        bias_model.meta.useafter = '2015-10-01T00:00:00'

        return data_model, bias_model

    return _cube


@pytest.fixture(scope='function')
def setup_subarray_cube():
    ''' Set up fake NIRCam subarray data to test.'''

    def _cube(xstart, ystart, ngroups, nrows, ncols):

        nints = 1

        # create a JWST datamodel for NIRCam SUB320A335R data
        data_model = RampModel((nints, ngroups, nrows, ncols))
        data_model.meta.subarray.name = 'SUB320A335R'
        data_model.meta.subarray.xstart = xstart
        data_model.meta.subarray.ystart = ystart
        data_model.meta.subarray.xsize = ncols
        data_model.meta.subarray.ysize = nrows
        data_model.meta.exposure.ngroups = ngroups
        data_model.meta.instrument.name = 'NIRCAM'
        data_model.meta.instrument.detector = 'NRCALONG'
        data_model.meta.observation.date = '2019-10-14'
        data_model.meta.observation.time = '16:44:12.000'

        # create a superbias model for the superbias step
        bias_model = SuperBiasModel((2048, 2048))
        bias_model.meta.subarray.xstart = 1
        bias_model.meta.subarray.ystart = 1
        bias_model.meta.subarray.xsize = 2048
        bias_model.meta.subarray.ysize = 2048
        bias_model.meta.instrument.name = 'NIRCAM'
        bias_model.meta.description = 'Fake data.'
        bias_model.meta.telescope = 'JWST'
        bias_model.meta.reftype = 'SuperBiasModel'
        bias_model.meta.author = 'Alicia'
        bias_model.meta.pedigree = 'Dummy'
        bias_model.meta.useafter = '2015-10-01T00:00:00'

        return data_model, bias_model

    return _cube
