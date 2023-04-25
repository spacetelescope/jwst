"""
Unit tests for reset anomaly correction
"""

import pytest
import numpy as np

from stdatamodels.jwst.datamodels import RampModel, ResetModel, dqflags

from jwst.reset.reset_sub import (
    do_correction as resetcorr
)


def test_correction(make_rampmodel, make_resetmodel):
    '''Check that data is unchanged for science data for groups > number of groups of reset '''

    # size of integration
    nints = 1
    ngroups = 15
    xsize = 200
    ysize = 200

    # create raw input data for step: ramp = 0 to 14
    dm_ramp = make_rampmodel(nints, ngroups, ysize, xsize)
    for i in range(0, ngroups):
        dm_ramp.data[0, i] = i

    refgroups = 12
    # create reset reference file model with frames less than science data
    reset = make_resetmodel(nints, refgroups, ysize, xsize)

    # populate data array of reference file
    for i in range(0, refgroups):
        reset.data[0, i] = i

    # set up test array
    test = make_rampmodel(nints, ngroups, ysize, xsize)
    # populate test data  = first 12 groups = 0, 11
    for i in range(0, refgroups):
        test.data[0, i, :, :] = 0

    # next groups =  12,13,14,15
    for i in range(refgroups, ngroups):
        test.data[0, i, :, :] = i

    # apply correction
    outfile = resetcorr(dm_ramp, reset)

    # test that the science data  are corrected for this first refgroups - should be 0
    # refgroups to ngroups should be = group # -1
    np.testing.assert_array_equal(outfile.data, test.data)


def test_nan(make_rampmodel, make_resetmodel):
    '''Verify that when a reset  has NaNs, these are correctly assumed as zero and the PIXELDQ is set properly'''

    # size of integration
    nints = 1
    ngroups = 10
    xsize = 200
    ysize = 200

    # create raw input data for step
    dm_ramp = make_rampmodel(nints, ngroups, ysize, xsize)

    # populate data array of science cube
    for i in range(0, ngroups - 1):
        dm_ramp.data[0, i, :, :] = i

    # create reset reference file model with more frames than science data
    refgroups = 15
    reset = make_resetmodel(nints, refgroups, ysize, xsize)

    # populate data array of reference file
    for i in range(0, refgroups - 1):
        reset.data[0, i] = i * 0.1

    # set NaN in reset data
    reset.data[0, 5, 100, 100] = np.nan

    # apply correction
    outfile = resetcorr(dm_ramp, reset)

    # test that the NaN reset reference pixel was set to 0 (nothing subtracted)

    assert outfile.data[0, 5, 100, 100] == 5.0


def test_dq_combine(make_rampmodel, make_resetmodel):
    '''Verify that the DQ array of the reset is correctly combined with the PIXELDQ array of the science data.'''

    # size of integration
    nints = 1
    ngroups = 5
    xsize = 200
    ysize = 200

    # create raw input data for step
    dm_ramp = make_rampmodel(nints, ngroups, ysize, xsize)

    # populate data array of science cube
    for i in range(1, ngroups - 1):
        dm_ramp.data[0, i, :, :] = i

    # create reset reference file model with more frames than science data
    refgroups = 12
    reset = make_resetmodel(nints, refgroups, ysize, xsize)

    jump_det = dqflags.pixel['JUMP_DET']
    saturated = dqflags.pixel['SATURATED']
    do_not_use = dqflags.pixel['DO_NOT_USE']
    nonlinear = dqflags.pixel['NONLINEAR']

    # populate dq flags of sci pixeldq and reference dq
    dm_ramp.pixeldq[50, 50] = jump_det
    dm_ramp.pixeldq[50, 51] = saturated

    reset.dq[50, 50] = np.bitwise_or(do_not_use, nonlinear)
    reset.dq[50, 51] = nonlinear

    # run correction step
    outfile = resetcorr(dm_ramp, reset)

    t50_50 = jump_det | do_not_use | nonlinear
    t50_51 = saturated | nonlinear
    # check that dq flags were correctly added
    assert outfile.pixeldq[50, 50] == t50_50
    assert outfile.pixeldq[50, 51] == t50_51


def test_2_int(make_rampmodel, make_resetmodel):
    '''Verify the reset correction is done by integration for MIRI observations'''

    # size of integration
    nints = 2
    ngroups = 10
    xsize = 200
    ysize = 200

    # create raw input data for step
    dm_ramp = make_rampmodel(nints, ngroups, ysize, xsize)

    # populate data array of science cube
    for i in range(0, ngroups - 1):
        dm_ramp.data[:, i] = i

    # create reset reference file model with more frames than science data
    refgroups = 10
    reset = make_resetmodel(nints, refgroups, ysize, xsize)

    # populate data array of reference file
    for i in range(0, refgroups - 1):
        reset.data[0, i] = i * 0.1
        reset.data[1, i] = i * 0.2

    # run correction
    outfile = resetcorr(dm_ramp, reset)

    # check that the reset file is subtracted frame by frame from the science data
    diff = dm_ramp.data[0] - reset.data[0, :ngroups]
    diff_int2 = dm_ramp.data[1] - reset.data[1, :ngroups]

    # test that the output data file is equal to the difference found when subtracting ref file from sci file
    np.testing.assert_array_equal(outfile.data[0], diff)
    np.testing.assert_array_equal(outfile.data[1], diff_int2)


@pytest.fixture(scope='function')
def make_rampmodel():
    '''Make MIRI Ramp model for testing'''
    def _ramp(nints, ngroups, ysize, xsize):
        # create the data and groupdq arrays
        csize = (nints, ngroups, ysize, xsize)
        data = np.full(csize, 1.0)  # default = 1.0

        # create a JWST datamodel for MIRI data
        dm_ramp = RampModel(data=data)

        dm_ramp.meta.instrument.name = 'MIRI'
        dm_ramp.meta.observation.date = '2018-01-01'
        dm_ramp.meta.observation.time = '00:00:00'
        dm_ramp.meta.subarray.xstart = 1
        dm_ramp.meta.subarray.xsize = xsize
        dm_ramp.meta.subarray.ystart = 1
        dm_ramp.meta.subarray.ysize = ysize
        dm_ramp.meta.exposure.nints = nints
        dm_ramp.meta.exposure.ngroups = ngroups
        dm_ramp.meta.description = 'Fake data.'

        return dm_ramp

    return _ramp


@pytest.fixture(scope='function')
def make_resetmodel():
    '''Make MIRI Reset model for testing'''
    def _reset(nints, ngroups, ysize, xsize):
        # create the data and groupdq arrays
        nints = 2
        csize = (nints, ngroups, ysize, xsize)
        data = np.full(csize, 1.0)  # default = 1.0

        # create a JWST datamodel for MIRI data
        reset = ResetModel(data=data)
        reset.meta.exposure.nints = nints
        reset.meta.exposure.ngroups = ngroups
        reset.meta.instrument.name = 'MIRI'
        reset.meta.date = '2018-01-01'
        reset.meta.time = '00:00:00'
        reset.meta.description = 'Fake data.'
        reset.meta.reftype = 'ResetModel'
        reset.meta.author = 'Jane Morrison'
        reset.meta.pedigree = 'Dummy'
        reset.meta.useafter = '2015-10-01T00:00:00'
        return reset

    return _reset
