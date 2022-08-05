"""
Unit tests for dark current correction
"""

import pytest
import numpy as np
from numpy.testing import assert_allclose

from stcal.dark_current.dark_sub import average_dark_frames_3d as average_dark_frames
from stcal.dark_current.dark_sub import do_correction as darkcorr

from stcal.dark_current.dark_class import DarkData

from jwst.dark_current.dark_current_step import DarkCurrentStep

from jwst.datamodels import RampModel, DarkModel, DarkMIRIModel, dqflags


# Define frame_time and number of groups in the generated dark reffile
TFRAME = 10.73677
NGROUPS_DARK = 10

DELIM = "-" * 80


@pytest.fixture(scope='function')
def setup_nrc_cube():
    '''Set up fake NIRCam data to test.'''

    def _cube(readpatt, ngroups, nframes, groupgap, nrows, ncols):

        nints = 1

        data_model = RampModel((nints, ngroups, nrows, ncols))
        data_model.meta.subarray.xstart = 1
        data_model.meta.subarray.ystart = 1
        data_model.meta.subarray.xsize = ncols
        data_model.meta.subarray.ysize = nrows
        data_model.meta.exposure.ngroups = ngroups
        data_model.meta.exposure.groupgap = groupgap
        data_model.meta.exposure.nframes = nframes
        data_model.meta.exposure.frame_time = TFRAME
        data_model.meta.exposure.group_time = (nframes + groupgap) * TFRAME
        data_model.meta.instrument.name = 'NIRCAM'
        data_model.meta.instrument.detector = 'NRCA1'
        data_model.meta.observation.date = '2017-10-01'
        data_model.meta.observation.time = '00:00:00'

        dark_model = DarkModel((NGROUPS_DARK, 2048, 2048))
        dark_model.meta.subarray.xstart = 1
        dark_model.meta.subarray.ystart = 1
        dark_model.meta.subarray.xsize = 2048
        dark_model.meta.subarray.ysize = 2048
        dark_model.meta.exposure.ngroups = NGROUPS_DARK
        dark_model.meta.exposure.groupgap = 0
        dark_model.meta.exposure.nframes = 1
        dark_model.meta.instrument.name = 'NIRCAM'
        dark_model.meta.description = 'Fake data.'
        dark_model.meta.telescope = 'JWST'
        dark_model.meta.reftype = 'DarkModel'
        dark_model.meta.author = 'Alicia'
        dark_model.meta.pedigree = 'Dummy'
        dark_model.meta.useafter = '2015-10-01T00:00:00'

        return data_model, dark_model

    return _cube


def _params():
    """Returns list of tuples, one for each readpatt, generating parameters for
    test_frame_averaging. Parameters are the following:

        (readpatt, ngroups, nframes, groupgap, nrows, ncols)

    Note groupgap = nskip
    """

    # Dictionary of NIRCam readout patterns
    readpatterns = dict(
        DEEP8=dict(ngroups=20, nframes=8, nskip=12),
        DEEP2=dict(ngroups=20, nframes=2, nskip=18),
        MEDIUM8=dict(ngroups=10, nframes=8, nskip=2),
        MEDIUM2=dict(ngroups=10, nframes=2, nskip=8),
        SHALLOW4=dict(ngroups=10, nframes=4, nskip=1),
        SHALLOW2=dict(ngroups=10, nframes=2, nskip=3),
        BRIGHT2=dict(ngroups=10, nframes=2, nskip=0),
        BRIGHT1=dict(ngroups=10, nframes=1, nskip=1),
        RAPID=dict(ngroups=10, nframes=1, nskip=0),
    )

    params = []
    ngroups = 3
    # NIRCam is 2048x2048, but we reduce the size to 20x20 for speed/memory
    nrows = 20
    ncols = 20
    for readpatt, values in readpatterns.items():
        params.append((readpatt, ngroups, values['nframes'], values['nskip'], nrows, ncols))

    return params


# Refac done
@pytest.mark.parametrize('readpatt, ngroups, nframes, groupgap, nrows, ncols', _params())
def test_frame_averaging(setup_nrc_cube, readpatt, ngroups, nframes, groupgap, nrows, ncols):
    '''Check that if nframes>1 or groupgap>0, then the pipeline reconstructs
       the dark reference file to match the frame averaging and groupgap
       settings of the exposure.'''

    # Create data and dark model
    data, dark_model = setup_nrc_cube(readpatt, ngroups, nframes, groupgap, nrows, ncols)
    dark = DarkData(dark_model=dark_model)

    # Add ramp values to dark model data array
    dark.data[:, 10, 10] = np.arange(0, NGROUPS_DARK)
    dark.err[:, 10, 10] = np.arange(10, NGROUPS_DARK + 10)

    # Run the pipeline's averaging function
    avg_dark = average_dark_frames(dark, ngroups, nframes, groupgap)

    # Group input groups into collections of frames which will be averaged
    total_frames = (nframes * ngroups) + (groupgap * (ngroups - 1))

    # Get starting/ending indexes of the input groups to be averaged
    gstrt_ind = np.arange(0, total_frames, nframes + groupgap)
    gend_ind = gstrt_ind + nframes

    # Prepare arrays to hold results of averaging
    manual_avg = np.zeros((ngroups), dtype=np.float32)
    manual_errs = np.zeros((ngroups), dtype=np.float32)

    # Manually average the input data to compare with pipeline output
    for newgp, gstart, gend in zip(range(ngroups), gstrt_ind, gend_ind):

        # Average the data frames
        newframe = np.mean(dark.data[gstart:gend, 10, 10])
        manual_avg[newgp] = newframe

        # ERR arrays will be quadratic sum of error values
        manual_errs[newgp] = np.sqrt(np.sum(dark.err[gstart:gend, 10, 10]**2)) / (gend - gstart)

    # Check that pipeline output matches manual averaging results
    assert_allclose(manual_avg, avg_dark.data[:, 10, 10], rtol=1e-5)
    assert_allclose(manual_errs, avg_dark.err[:, 10, 10], rtol=1e-5)

    # Check that meta data was properly updated
    assert avg_dark.exp_nframes == nframes
    assert avg_dark.exp_ngroups == ngroups
    assert avg_dark.exp_groupgap == groupgap


# Refac done
def test_more_sci_frames(make_rampmodel, make_darkmodel):
    '''Check that data is unchanged if there are more frames in the science
    data is than in the dark reference file and verify that when the dark is not applied,
    the data is correctly flagged as such'''

    # size of integration
    nints = 1
    ngroups = 7
    xsize = 200
    ysize = 200

    # create raw input data for step
    dm_ramp = make_rampmodel(nints, ngroups, ysize, xsize)
    dm_ramp.meta.exposure.nframes = 1
    dm_ramp.meta.exposure.groupgap = 0

    # populate data array of science cube
    for i in range(0, ngroups - 1):
        dm_ramp.data[0, i] = i

    refgroups = 5
    # create dark reference file model with fewer frames than science data
    dark = make_darkmodel(refgroups, ysize, xsize)

    # populate data array of reference file
    for i in range(0, refgroups - 1):
        dark.data[0, i] = i * 0.1

    # apply correction
    out_data, avg_data = darkcorr(dm_ramp, dark)

    # test that the science data are not changed; input file = output file
    np.testing.assert_array_equal(dm_ramp.data, out_data.data)

    # get dark correction status from header
    # darkstatus = outfile.meta.cal_step.dark_sub
    darkstatus = out_data.cal_step

    assert darkstatus == 'SKIPPED'


# Refac done
def test_sub_by_frame(make_rampmodel, make_darkmodel):
    '''Check that if NFRAMES=1 and GROUPGAP=0 for the science data, the dark reference data are
    directly subtracted frame by frame'''

    # size of integration
    nints = 1
    ngroups = 10
    xsize = 200
    ysize = 200

    # create raw input data for step
    dm_ramp = make_rampmodel(nints, ngroups, ysize, xsize)
    dm_ramp.meta.exposure.nframes = 1
    dm_ramp.meta.exposure.groupgap = 0

    # populate data array of science cube
    for i in range(0, ngroups - 1):
        dm_ramp.data[0, i] = i

    # create dark reference file model with more frames than science data
    refgroups = 15
    dark = make_darkmodel(refgroups, ysize, xsize)

    # populate data array of reference file
    for i in range(0, refgroups - 1):
        dark.data[0, i] = i * 0.1

    # apply correction
    outfile, avg_dark = darkcorr(dm_ramp, dark)

    # remove the single dimension at start of file (1, 30, 1032, 1024) so comparison in assert works
    outdata = np.squeeze(outfile.data)

    # check that the dark file is subtracted frame by frame from the science data
    diff = dm_ramp.data[0] - dark.data[0, :ngroups]

    # test that the output data file is equal to the difference found when subtracting ref file from sci file
    np.testing.assert_array_equal(outdata, diff, err_msg='dark file should be subtracted from sci file ')


# Refac done
def test_nan(make_rampmodel, make_darkmodel):
    '''Verify that when a dark has NaNs, these are correctly assumed as zero and the PIXELDQ is set properly'''

    # size of integration
    nints = 1
    ngroups = 10
    xsize = 200
    ysize = 200

    # create raw input data for step
    dm_ramp = make_rampmodel(nints, ngroups, ysize, xsize)
    dm_ramp.meta.exposure.nframes = 1
    dm_ramp.meta.exposure.groupgap = 0

    # populate data array of science cube
    for i in range(0, ngroups - 1):
        dm_ramp.data[0, i, :, :] = i

    # create dark reference file model with more frames than science data
    refgroups = 15
    dark = make_darkmodel(refgroups, ysize, xsize)

    # populate data array of reference file
    for i in range(0, refgroups - 1):
        dark.data[0, i] = i * 0.1

    # set NaN in dark file

    dark.data[0, 5, 100, 100] = np.nan

    # apply correction
    outfile, avg_dark = darkcorr(dm_ramp, dark)

    # test that the NaN dark reference pixel was set to 0 (nothing subtracted)
    assert outfile.data[0, 5, 100, 100] == 5.0


# Refac done
def test_dq_combine(make_rampmodel, make_darkmodel):
    '''Verify that the DQ array of the dark is correctly combined with the PIXELDQ array of the science data.'''

    # size of integration
    nints = 1
    ngroups = 5
    xsize = 200
    ysize = 200

    # create raw input data for step
    dm_ramp = make_rampmodel(nints, ngroups, ysize, xsize)
    dm_ramp.meta.exposure.nframes = 1
    dm_ramp.meta.exposure.groupgap = 0

    # populate data array of science cube
    for i in range(1, ngroups - 1):
        dm_ramp.data[0, i, :, :] = i

    # create dark reference file model with more frames than science data
    refgroups = 7
    dark = make_darkmodel(refgroups, ysize, xsize)

    jump_det = dqflags.pixel['JUMP_DET']
    saturated = dqflags.pixel['SATURATED']
    do_not_use = dqflags.pixel['DO_NOT_USE']

    # populate dq flags of sci pixeldq and reference dq
    dm_ramp.pixeldq[50, 50] = jump_det
    dm_ramp.pixeldq[50, 51] = saturated

    dark.dq[0, 0, 50, 50] = do_not_use
    dark.dq[0, 0, 50, 51] = do_not_use

    # run correction step
    outfile, avg_dark = darkcorr(dm_ramp, dark)

    # check that dq flags were correctly added
    assert outfile.pixeldq[50, 50] == np.bitwise_or(jump_det, do_not_use)
    assert outfile.pixeldq[50, 51] == np.bitwise_or(saturated, do_not_use)


# Refac done
def test_2_int(make_rampmodel, make_darkmodel):
    '''Verify the dark correction is done by integration for MIRI observations'''

    # size of integration
    nints = 2
    ngroups = 10
    xsize = 200
    ysize = 200

    # create raw input data for step
    dm_ramp = make_rampmodel(nints, ngroups, ysize, xsize)
    dm_ramp.meta.exposure.nframes = 1
    dm_ramp.meta.exposure.groupgap = 0

    # populate data array of science cube
    for i in range(0, ngroups - 1):
        dm_ramp.data[:, i] = i

    # create dark reference file model with more frames than science data
    refgroups = 15
    dark = make_darkmodel(refgroups, ysize, xsize)

    # populate data array of reference file
    for i in range(0, refgroups - 1):
        dark.data[0, i] = i * 0.1
        dark.data[1, i] = i * 0.2

    # run correction
    outfile, avg_dark = darkcorr(dm_ramp, dark)

    # check that the dark file is subtracted frame by frame from the science data
    diff = dm_ramp.data[0] - dark.data[0, :ngroups]
    diff_int2 = dm_ramp.data[1] - dark.data[1, :ngroups]

    # test that the output data file is equal to the difference found when subtracting ref file from sci file
    np.testing.assert_array_equal(outfile.data[0], diff)
    np.testing.assert_array_equal(outfile.data[1], diff_int2)


# Refac done
def test_frame_avg(make_rampmodel, make_darkmodel):
    '''Check that if NFRAMES>1 or GROUPGAP>0, the frame-averaged dark data are
    subtracted group-by-group from science data groups and the ERR arrays are not modified'''

    # size of integration
    nints = 1
    ngroups = 5
    xsize = 1032
    ysize = 1024

    # create raw input data for step
    dm_ramp = make_rampmodel(nints, ngroups, ysize, xsize)
    dm_ramp.meta.exposure.nframes = 4
    dm_ramp.meta.exposure.groupgap = 0

    # populate data array of science cube
    for i in range(0, ngroups - 1):
        dm_ramp.data[:, i] = i + 1

    # create dark reference file model

    refgroups = 20  # This needs to be 20 groups for the calculations to work
    dark = make_darkmodel(refgroups, ysize, xsize)

    # populate data array of reference file
    for i in range(0, refgroups - 1):
        dark.data[0, i] = i * 0.1

    # apply correction
    outfile, avg_dark = darkcorr(dm_ramp, dark)

    # dark frames should be averaged in groups of 4 frames

    # this will result in average values of 0.15, 0.55, 0.95, and 1.35
    # these values are then subtracted from frame values of 1, 2, 3 and 4

    assert outfile.data[0, 0, 500, 500] == pytest.approx(0.85)
    assert outfile.data[0, 1, 500, 500] == pytest.approx(1.45)
    assert outfile.data[0, 2, 500, 500] == pytest.approx(2.05)
    assert outfile.data[0, 3, 500, 500] == pytest.approx(2.65)

    # check that the error array is not modified.
    np.testing.assert_array_equal(outfile.err[:, :], 0)


# ------------------------------------------------------------------------------
def test_basic_step(make_rampmodel, make_darkmodel):
    """
    Same as test_more_sci_frames above, but done calling the step code.
    """
    # size of integration
    nints, ngroups, nrows, ncols = 1, 10, 200, 200

    # create raw input data for step
    dm_ramp = make_rampmodel(nints, ngroups, nrows, ncols)
    dm_ramp.meta.exposure.nframes = 1
    dm_ramp.meta.exposure.groupgap = 0

    # populate data array of science cube
    for i in range(0, ngroups - 1):
        dm_ramp.data[0, i] = i

    # create dark reference file model with more frames than science data
    refgroups = 15
    dark = make_darkmodel(refgroups, nrows, ncols)

    # populate data array of reference file
    for i in range(0, refgroups - 1):
        dark.data[0, i] = i * 0.1

    dark_model = DarkCurrentStep.call(dm_ramp, override_dark=dark)

    assert dark_model.meta.cal_step.dark_sub == "COMPLETE"

    outdata = np.squeeze(dark_model.data)

    # check that the dark file is subtracted frame by frame from the science data
    diff = dm_ramp.data[0] - dark.data[0, :ngroups]

    # test that the output data file is equal to the difference found when subtracting ref file from sci file
    np.testing.assert_array_equal(outdata, diff, err_msg='dark file should be subtracted from sci file ')


@pytest.fixture(scope='function')
def make_rampmodel():
    '''Make MIRI Ramp model for testing'''

    def _ramp(nints, ngroups, ysize, xsize):

        # create the data and groupdq arrays
        csize = (nints, ngroups, ysize, xsize)
        data = np.full(csize, 1.0)

        # create a JWST datamodel for MIRI data
        dm_ramp = RampModel(data=data)

        dm_ramp.meta.instrument.name = 'MIRI'
        dm_ramp.meta.observation.date = '2018-01-01'
        dm_ramp.meta.observation.time = '00:00:00'
        dm_ramp.meta.subarray.xstart = 1
        dm_ramp.meta.subarray.xsize = xsize
        dm_ramp.meta.subarray.ystart = 1
        dm_ramp.meta.subarray.ysize = ysize
        dm_ramp.meta.description = 'Fake data.'

        return dm_ramp

    return _ramp


@pytest.fixture(scope='function')
def make_darkmodel():
    '''Make MIRI dark model for testing'''

    def _dark(ngroups, ysize, xsize):
        # create the data and groupdq arrays
        nints = 2
        csize = (nints, ngroups, ysize, xsize)
        data = np.full(csize, 1.0)

        # create a JWST datamodel for MIRI data
        dark = DarkMIRIModel(data=data)

        dark.meta.instrument.name = 'MIRI'
        dark.meta.date = '2018-01-01'
        dark.meta.time = '00:00:00'
        dark.meta.subarray.xstart = 1
        dark.meta.subarray.xsize = xsize
        dark.meta.subarray.ystart = 1
        dark.meta.subarray.ysize = ysize
        dark.meta.exposure.nframes = 1
        dark.meta.exposure.groupgap = 0
        dark.meta.description = 'Fake data.'
        dark.meta.reftype = 'DarkModel'
        dark.meta.author = 'Alicia'
        dark.meta.pedigree = 'Dummy'
        dark.meta.useafter = '2015-10-01T00:00:00'

        return dark

    return _dark
