"""

Unit tests for saturation flagging

"""

import pytest
import numpy as np

from stdatamodels.jwst.datamodels import RampModel, SaturationModel, dqflags

from jwst.saturation import SaturationStep
from jwst.saturation.saturation import flag_saturation, irs2_flag_saturation


def test_basic_saturation_flagging(setup_nrc_cube):
    '''Check that the saturation flag is set when a pixel value is above the
       threshold given by the reference file.'''

    # Create inputs, data, and saturation maps
    ngroups = 5
    nrows = 20
    ncols = 20
    satvalue = 60000

    data, satmap = setup_nrc_cube(ngroups, nrows, ncols)

    # Add ramp values up to the saturation limit
    data.data[0, 0, 5, 5] = 0
    data.data[0, 1, 5, 5] = 20000
    data.data[0, 2, 5, 5] = 40000
    data.data[0, 3, 5, 5] = 60000   # Signal reaches saturation limit
    data.data[0, 4, 5, 5] = 62000

    # Set saturation value in the saturation model
    satmap.data[5, 5] = satvalue

    # Run the pipeline
    output = flag_saturation(data, satmap, n_pix_grow_sat=1)

    # Make sure that groups with signal > saturation limit get flagged
    satindex = np.argmax(output.data[0, :, 4:6, 4:6] == satvalue)
    assert np.all(output.groupdq[0, satindex:, 4:6, 4:6] == dqflags.group['SATURATED'])


def test_nirspec_irs2_saturation_flagging(setup_nrs_irs2_cube):
    data, satmap = setup_nrs_irs2_cube()

    pixx, pixy = 1000, 1000
    data.data[0, 3, pixx, pixy] = 65000  # Signal exceeds saturation limit of 60000
    data.data[0, 4, pixx, pixy] = 67000

    # Run saturation detection
    output = irs2_flag_saturation(data, satmap, n_pix_grow_sat=1)

    # Make sure that groups with signal > saturation limit get flagged
    assert np.all(output.groupdq[0, 3:, pixx - 1: pixx + 1, pixy - 1: pixy + 1] == dqflags.group['SATURATED'])


def test_irs2_zero_frame(setup_nrs_irs2_cube):
    """
    Tests IRS2 ZEROFRAME processing.

    This ensures ZEROFRAME data that is outside saturation and AD floor
    boundaries get zeroed out.
    """
    ramp, sat = setup_nrs_irs2_cube()

    # Setup ZEROFRAME
    nints, ngroups, nrows, ncols = ramp.data.shape
    dims = (nints, nrows, ncols)
    ramp.zeroframe = np.ones(dims, dtype=ramp.data.dtype) * 1000.

    nint, row, col = 0, 1000, 1000
    ramp.meta.exposure.zero_frame = True
    ramp.zeroframe[nint, row, col] = 65000.
    ramp.zeroframe[nint, row + 1, col + 1] = -100.

    output = irs2_flag_saturation(ramp, sat, n_pix_grow_sat=0)

    check_zframe = np.ones(dims, dtype=ramp.data.dtype) * 1000.
    check_zframe[nint, row, col] = 0.
    check_zframe[nint, row + 1, col + 1] = 0.
    tol = 1.e-5
    np.testing.assert_allclose(output.zeroframe, check_zframe, tol)


def test_ad_floor_flagging(setup_nrc_cube):
    """Check that the ad_floor flag is set when a pixel value is zero or
    negative."""

    # Create inputs, data, and saturation maps
    ngroups = 5
    nrows = 20
    ncols = 20
    satvalue = 60000

    data, satmap = setup_nrc_cube(ngroups, nrows, ncols)

    # Add ramp values up to the saturation limit
    data.data[0, 0, 5, 5] = 0  # Signal at bottom rail - low saturation
    data.data[0, 1, 5, 5] = 0  # Signal at bottom rail - low saturation
    data.data[0, 2, 5, 5] = 20
    data.data[0, 3, 5, 5] = 40
    data.data[0, 4, 5, 5] = 60

    # frames that should be flagged as saturation (low)
    satindxs = [0, 1]

    # Set saturation value in the saturation model
    satmap.data[5, 5] = satvalue

    # Run the pipeline
    output = flag_saturation(data, satmap, n_pix_grow_sat=1)

    # Check if the right frames are flagged as saturated
    assert np.all(output.groupdq[0, satindxs, 5, 5]
                  == dqflags.group['DO_NOT_USE'] | dqflags.group['AD_FLOOR'])


def test_ad_floor_and_saturation_flagging(setup_nrc_cube):
    """Check that the ad_floor flag is set when a pixel value is zero or
    negative and the saturation flag when the pixel is above the saturation threshold."""

    # Create inputs, data, and saturation maps
    ngroups = 5
    nrows = 20
    ncols = 20
    satvalue = 60000

    data, satmap = setup_nrc_cube(ngroups, nrows, ncols)

    # Add ramp values up to the saturation limit
    data.data[0, 0, 5, 5] = 0  # Signal at bottom rail - low saturation
    data.data[0, 1, 5, 5] = 0  # Signal at bottom rail - low saturation
    data.data[0, 2, 5, 5] = 20
    data.data[0, 3, 5, 5] = 40
    data.data[0, 4, 5, 5] = 61000  # Signal above the saturation threshold

    # frames that should be flagged as ad_floor
    floorindxs = [0, 1]
    # frames that should be flagged as saturation
    satindxs = [4]

    # Set saturation value in the saturation model
    satmap.data[5, 5] = satvalue

    # Run the pipeline
    output = flag_saturation(data, satmap, n_pix_grow_sat=1)

    # Check if the right frames are flagged as ad_floor
    assert np.all(output.groupdq[0, floorindxs, 5, 5]
                  == dqflags.group['DO_NOT_USE'] | dqflags.group['AD_FLOOR'])
    # Check if the right frames are flagged as saturated
    assert np.all(np.bitwise_and(output.groupdq[0, satindxs, 4:6, 4:6],
                  dqflags.group['SATURATED']) == dqflags.group['SATURATED'])


def test_signal_fluctuation_flagging(setup_nrc_cube):
    '''Check that once a pixel is flagged as saturated in a group, all
       subsequent groups should also be flagged as saturated, even if the
       signal value drops back below saturation.'''

    # Create inputs, data, and saturation maps
    ngroups = 5
    nrows = 20
    ncols = 20
    satvalue = 60000

    data, satmap = setup_nrc_cube(ngroups, nrows, ncols)

    # Add ramp values up to the saturation limit
    data.data[0, 0, 5, 5] = 10
    data.data[0, 1, 5, 5] = 20000
    data.data[0, 2, 5, 5] = 40000
    data.data[0, 3, 5, 5] = 60000   # Signal reaches saturation limit
    data.data[0, 4, 5, 5] = 40000   # Signal drops below saturation limit

    # Set saturation value in the saturation model
    satmap.data[5, 5] = satvalue

    # Run the pipeline
    output = flag_saturation(data, satmap, n_pix_grow_sat=1)

    # Make sure that all groups after first saturated group are flagged
    satindex = np.argmax(output.data[0, :, 5, 5] == satvalue)
    assert np.all(np.bitwise_and(output.groupdq[0, satindex:, 4:6, 4:6],
                                 dqflags.group['SATURATED']) == dqflags.group['SATURATED'])


def test_all_groups_saturated(setup_nrc_cube):
    '''Check case where all groups are saturated.'''

    # Create inputs, data, and saturation maps
    ngroups = 5
    nrows = 20
    ncols = 20
    satvalue = 60000

    data, satmap = setup_nrc_cube(ngroups, nrows, ncols)

    # Add ramp values at or above saturation limit
    data.data[0, 0, 5, 5] = 60000
    data.data[0, 1, 5, 5] = 62000
    data.data[0, 2, 5, 5] = 62000
    data.data[0, 3, 5, 5] = 60000
    data.data[0, 4, 5, 5] = 62000

    # Set saturation value in the saturation model
    satmap.data[5, 5] = satvalue

    # Run the pipeline
    output = flag_saturation(data, satmap, n_pix_grow_sat=1)

    # Make sure all groups are flagged
    assert np.all(np.bitwise_and(output.groupdq[0, :, 4:6, 4:6],
                                 dqflags.group['SATURATED']) == dqflags.group['SATURATED'])


def test_subarray_extraction(setup_miri_cube):
    '''Check the step correctly handles subarrays.'''

    # Create input data
    # Create model of data with 0 value array
    ngroups = 50
    nrows = 224
    ncols = 288

    data, satmap = setup_miri_cube(1, 467, ngroups, nrows, ncols)

    # Place DQ flags in DQ array that would be in subarray
    # MASK1550 file has colstart=1, rowstart=467
    satmap.dq[542, 100:105] = dqflags.pixel['DO_NOT_USE']

    # Test a value of NaN in the reference file with an existing DQ flag
    satmap.data[550, 100] = np.nan
    satmap.dq[550, 100] = dqflags.pixel['NONLINEAR']

    # Run the pipeline
    output = flag_saturation(data, satmap, n_pix_grow_sat=1)

    # Check for DQ flag in PIXELDQ of subarray image
    assert output.pixeldq[76, 100] == dqflags.pixel['DO_NOT_USE']
    assert output.pixeldq[76, 104] == dqflags.pixel['DO_NOT_USE']

    # Pixel 84, 100 in subarray maps to 550, 100 in reference file
    # Check that pixel was flagged 'NO_SAT_CHECK' and that original
    # DQ flag persists (i.e. did not get overwritten)
    assert (output.pixeldq[84, 100] ==
           dqflags.pixel['NO_SAT_CHECK'] + dqflags.pixel['NONLINEAR'])


def test_dq_propagation(setup_nrc_cube):
    '''Check PIXELDQ propagation.'''

    # Create inputs, data, and saturation maps
    ngroups = 5
    nrows = 20
    ncols = 20
    dqval1 = 5
    dqval2 = 10

    data, satmap = setup_nrc_cube(ngroups, nrows, ncols)

    # Add DQ values to the data and reference file
    data.pixeldq[5, 5] = dqval1
    satmap.dq[5, 5] = dqval2

    # Run the pipeline
    output = flag_saturation(data, satmap, n_pix_grow_sat=1)

    # Make sure DQ values from data and reference file are added in the output
    assert output.pixeldq[5, 5] == dqval1 + dqval2


def test_no_sat_check(setup_nrc_cube):
    '''Check that pixels flagged with NO_SAT_CHECK in the reference file get
       added to the DQ mask and are not flagged as saturated.'''

    # Create inputs, data, and saturation maps
    ngroups = 5
    nrows = 20
    ncols = 20
    satvalue = 60000

    data, satmap = setup_nrc_cube(ngroups, nrows, ncols)

    # Add ramp values up to the saturation limit
    data.data[0, 0, 5, 5] = 10
    data.data[0, 1, 5, 5] = 20000
    data.data[0, 2, 5, 5] = 40000
    data.data[0, 3, 5, 5] = 60000
    data.data[0, 4, 5, 5] = 62000   # Signal reaches saturation limit

    # Set saturation value in the saturation model & DQ value for NO_SAT_CHECK
    satmap.data[5, 5] = satvalue
    satmap.dq[5, 5] = dqflags.pixel['NO_SAT_CHECK']

    # Also set an existing DQ flag in input science data
    data.pixeldq[5, 5] = dqflags.pixel['RC']

    # Run the pipeline
    output = flag_saturation(data, satmap, n_pix_grow_sat=1)

    # Make sure output GROUPDQ does not get flagged as saturated
    # Make sure PIXELDQ is set to NO_SAT_CHECK and original flag
    assert np.all(output.groupdq[0, :, 4:6, 4:6] != dqflags.group['SATURATED'])
    assert output.pixeldq[5, 5] == (dqflags.pixel['NO_SAT_CHECK'] +
                                    dqflags.pixel['RC'])


def test_nans_in_mask(setup_nrc_cube):
    '''Check that pixels in the reference files that have value NaN are not
       flagged as saturated in the data and that in the PIXELDQ array the
       pixel is set to NO_SAT_CHECK.'''

    # Create inputs, data, and saturation maps
    ngroups = 5
    nrows = 20
    ncols = 20

    data, satmap = setup_nrc_cube(ngroups, nrows, ncols)

    # Add ramp values up to the saturation limit
    data.data[0, 0, 5, 5] = 10
    data.data[0, 1, 5, 5] = 20000
    data.data[0, 2, 5, 5] = 40000
    data.data[0, 3, 5, 5] = 60000
    data.data[0, 4, 5, 5] = 62000

    # Set saturation value for pixel to NaN
    satmap.data[5, 5] = np.nan

    # Run the pipeline
    output = flag_saturation(data, satmap, n_pix_grow_sat=1)

    # Check that output GROUPDQ is not flagged as saturated
    assert np.all(output.groupdq[0, :, 4:6, 4:6] != dqflags.group['SATURATED'])
    # Check that output PIXELDQ is set to NO_SAT_CHECK
    assert output.pixeldq[5, 5] == dqflags.pixel['NO_SAT_CHECK']


def test_full_step(setup_nrc_cube):
    '''Test full run of the SaturationStep.'''

    # Create inputs, data, and saturation maps
    ngroups = 5
    nrows = 20
    ncols = 20

    data, satmap = setup_nrc_cube(ngroups, nrows, ncols)

    # set the entire array to a small non-zero value to avoid labeling
    # almost everything as low saturated
    data.data[:, :, :, :] = 1

    # Add ramp values up to the saturation limit
    data.data[0, 0, 5, 5] = 10
    data.data[0, 1, 5, 5] = 20000
    data.data[0, 2, 5, 5] = 40000
    data.data[0, 3, 5, 5] = 70000   # Signal reaches saturation limit
    data.data[0, 4, 5, 5] = 73000

    # Run the pipeline
    output = SaturationStep.call(data)

    # Check that correct pixel and group 3+ are flagged as saturated
    assert dqflags.group['SATURATED'] == np.max(output.groupdq[0, :, 5, 5])
    assert np.all(output.groupdq[0, 3:, 5, 5] == dqflags.pixel['SATURATED'])
    # Check that other pixel and groups are not flagged
    # assert np.all(output.groupdq[0, :3, 5, 5] != dqflags.group['SATURATED'])
    # assert np.all(output.groupdq[0, :, 10, 10] != dqflags.group['SATURATED'])


@pytest.fixture(scope='function')
def setup_nrc_cube():
    ''' Set up fake NIRCam data to test.'''

    def _cube(ngroups, nrows, ncols):

        nints = 1

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

        saturation_model = SaturationModel((2048, 2048))
        saturation_model.meta.subarray.xstart = 1
        saturation_model.meta.subarray.ystart = 1
        saturation_model.meta.subarray.xsize = 2048
        saturation_model.meta.subarray.ysize = 2048
        saturation_model.meta.instrument.name = 'NIRCAM'
        saturation_model.meta.description = 'Fake data.'
        saturation_model.meta.telescope = 'JWST'
        saturation_model.meta.reftype = 'SaturationModel'
        saturation_model.meta.author = 'Alicia'
        saturation_model.meta.pedigree = 'Dummy'
        saturation_model.meta.useafter = '2015-10-01T00:00:00'

        return data_model, saturation_model

    return _cube


@pytest.fixture(scope='function')
def setup_miri_cube():
    ''' Set up fake MIRI data to test.'''

    def _cube(xstart, ystart, ngroups, nrows, ncols):

        nints = 1

        # create a JWST datamodel for MIRI data
        data_model = RampModel((nints, ngroups, nrows, ncols))
        data_model.data += 1
        data_model.meta.instrument.name = 'MIRI'
        data_model.meta.instrument.detector = 'MIRIMAGE'
        data_model.meta.instrument.filter = 'F1500W'
        data_model.meta.instrument.band = 'N/A'
        data_model.meta.observation.date = '2016-06-01'
        data_model.meta.observation.time = '00:00:00'
        data_model.meta.exposure.type = 'MIR_IMAGE'
        data_model.meta.subarray.name = 'MASK1550'
        data_model.meta.subarray.xstart = xstart
        data_model.meta.subarray.xsize = ncols
        data_model.meta.subarray.ystart = ystart
        data_model.meta.subarray.ysize = nrows

        # create a saturation model for the saturation step
        saturation_model = SaturationModel((1032, 1024))
        saturation_model.meta.description = 'Fake data.'
        saturation_model.meta.telescope = 'JWST'
        saturation_model.meta.reftype = 'SaturationModel'
        saturation_model.meta.author = 'Alicia'
        saturation_model.meta.pedigree = 'Dummy'
        saturation_model.meta.useafter = '2015-10-01T00:00:00'
        saturation_model.meta.instrument.name = 'MIRI'
        saturation_model.meta.instrument.detector = 'MIRIMAGE'
        saturation_model.meta.subarray.xstart = 1
        saturation_model.meta.subarray.xsize = 1024
        saturation_model.meta.subarray.ystart = 1
        saturation_model.meta.subarray.ysize = 1032

        return data_model, saturation_model

    return _cube


@pytest.fixture(scope='function')
def setup_nrs_irs2_cube():

    def _cube():

        # create a JWST datamodel for NIRSPEC IRS2 data
        data_model = RampModel((1, 5, 3200, 2048))
        data_model.data = np.ones(((1, 5, 3200, 2048)))
        data_model.groupdq = np.zeros(((1, 5, 3200, 2048)))
        data_model.pixeldq = np.zeros(((3200, 2048)))
        data_model.meta.instrument.name = 'NIRSPEC'
        data_model.meta.instrument.detector = 'NRS1'
        data_model.meta.instrument.filter = 'F100LP'
        data_model.meta.observation.date = '2019-07-19'
        data_model.meta.observation.time = '23:23:30.912'
        data_model.meta.exposure.type = 'NRS_LAMP'
        data_model.meta.subarray.name = 'FULL'
        data_model.meta.subarray.xstart = 1
        data_model.meta.subarray.xsize = 2048
        data_model.meta.subarray.ystart = 1
        data_model.meta.subarray.ysize = 2048
        data_model.meta.exposure.nrs_normal = 16
        data_model.meta.exposure.nrs_reference = 4
        data_model.meta.exposure.readpatt = 'NRSIRS2RAPID'

        # create a saturation model for the saturation step
        saturation_model = SaturationModel((2048, 2048))
        saturation_model.data = np.ones((2048, 2048)) * 60000  # saturation limit for every pixel is 60000
        saturation_model.meta.description = 'Fake data.'
        saturation_model.meta.telescope = 'JWST'
        saturation_model.meta.reftype = 'SaturationModel'
        saturation_model.meta.useafter = '2015-10-01T00:00:00'
        saturation_model.meta.instrument.name = 'NIRSPEC'
        saturation_model.meta.instrument.detector = 'NRS1'
        saturation_model.meta.author = 'Clare'
        saturation_model.meta.pedigree = 'Dummy'
        saturation_model.meta.subarray.xstart = 1
        saturation_model.meta.subarray.xsize = 2048
        saturation_model.meta.subarray.ystart = 1
        saturation_model.meta.subarray.ysize = 2048

        return data_model, saturation_model
    return _cube
