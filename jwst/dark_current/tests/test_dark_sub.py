"""
Unit tests for dark current correction
"""

import pytest
import numpy as np

from jwst.dark_current.dark_sub import average_dark_frames
from jwst.datamodels import (
    RampModel, 
    DarkModel, 
    MIRIRampModel, 
    DarkMIRIModel, 
    dqflags,
)
from jwst.dark_current.dark_sub import do_correction as darkcorr
from jwst.dark_current import DarkCurrentStep
from astropy.io import fits

# Dictionary of NIRCam readout patterns
DEEP8 = {}
DEEP8['ngroups'] = 20.
DEEP8['nframes'] = 8
DEEP8['nskip'] = 12

DEEP2 = {}
DEEP2['ngroups'] = 20.
DEEP2['nframes'] = 2
DEEP2['nskip'] = 18

MEDIUM8 = {}
MEDIUM8['ngroups'] = 10.
MEDIUM8['nframes'] = 8
MEDIUM8['nskip'] = 2

MEDIUM2 = {}
MEDIUM2['ngroups'] = 10.
MEDIUM2['nframes'] = 2
MEDIUM2['nskip'] = 8

SHALLOW4 = {}
SHALLOW4['ngroups'] = 10.
SHALLOW4['nframes'] = 4
SHALLOW4['nskip'] = 1

SHALLOW2 = {}
SHALLOW2['ngroups'] = 10.
SHALLOW2['nframes'] = 2
SHALLOW2['nskip'] = 3

BRIGHT2 = {}
BRIGHT2['ngroups'] = 10.
BRIGHT2['nframes'] = 2
BRIGHT2['nskip'] = 0

BRIGHT1 = {}
BRIGHT1['ngroups'] = 10.
BRIGHT1['nframes'] = 1
BRIGHT1['nskip'] = 1

RAPID = {}
RAPID['ngroups'] = 10.
RAPID['nframes'] = 1
RAPID['nskip'] = 0

READPATTERNS = {}
READPATTERNS['DEEP8'] = DEEP8
READPATTERNS['DEEP2'] = DEEP2
READPATTERNS['MEDIUM8'] = MEDIUM8
READPATTERNS['MEDIUM2'] = MEDIUM2
READPATTERNS['SHALLOW4'] = SHALLOW4
READPATTERNS['SHALLOW2'] = SHALLOW2
READPATTERNS['BRIGHT2'] = BRIGHT2
READPATTERNS['BRIGHT1'] = BRIGHT1
READPATTERNS['RAPID'] = RAPID

TFRAME = 10.73677


def test_frame_averaging(setup_nrc_cube):
    '''Check that if nframes>1 or groupgap>0, then the pipeline reconstructs
       the dark reference file to match the frame averaging and groupgap
       settings of the exposure.'''

    # Values to build the fake data arrays
    ngroups = 5
    nrows = 2048
    ncols = 2048

    # Loop over the NIRCam readout patterns:
    for name in READPATTERNS:

        # Get the configuration for the readout pattern
        readpatt = name
        nframes = READPATTERNS[readpatt]['nframes']
        groupgap = READPATTERNS[readpatt]['nskip']

        # Create data and dark model
        data, dark = setup_nrc_cube(readpatt, ngroups, nrows, ncols)

        # Add ramp values to dark model data array
        dark.data[:, 500, 500] = np.arange(0, 100)
        dark.err[:, 500, 500] = np.arange(100, 200)

        # Run the pipeline's averaging function
        avg_dark = average_dark_frames(dark, ngroups, nframes, groupgap)

        # Group input groups into collections of frames which will be averaged
        total_frames = (nframes * ngroups) + (groupgap * (ngroups-1))

        # Get starting/ending indexes of the input groups to be averaged
        gstrt_ind = np.arange(0, total_frames, nframes + groupgap)
        gend_ind = gstrt_ind + nframes

        # Prepare arrays to hold results of averaging
        manual_avg = np.zeros((ngroups))
        manual_errs = np.zeros((ngroups))

        # Manually average the input data to compare with pipeline output
        for newgp, gstart, gend in zip(range(ngroups), gstrt_ind, gend_ind):

            # Average the data frames
            newframe = np.mean(dark.data[gstart:gend, 500, 500])
            manual_avg[newgp] = newframe

            # ERR arrays will be quadratic sum of error values
            manual_errs[newgp] = np.sqrt(np.sum(dark.err[gstart:gend, 500, 500]**2)) / (gend - gstart)

        # Check that pipeline output matches manual averaging results
        assert np.all(manual_avg == avg_dark.data[:, 500, 500])
        assert np.all(manual_errs == avg_dark.err[:, 500, 500])

        # Check that meta data was properly updated
        assert avg_dark.meta.exposure.nframes == nframes
        assert avg_dark.meta.exposure.ngroups == ngroups
        assert avg_dark.meta.exposure.groupgap == groupgap



def test_more_sci_frames():
    '''Check that data is unchanged if there are more frames in the science data is than in the reference file'''

    # size of integration
    nints = 1
    ngroups = 30
    xsize = 1032
    ysize = 1024

    # create raw input data for step
    dm_ramp = make_rampmodel(nints, ngroups, ysize, xsize)
    dm_ramp.meta.exposure.nframes = 1
    dm_ramp.meta.exposure.groupgap = 0

    # populate data array of science cube
    for i in range(0, ngroups-1):
        dm_ramp.data[0, i, :, :] = i

    refgroups = 20
    # create dark reference file model with fewer frames than science data
    dark = make_darkmodel(refgroups, ysize, xsize)

    # populate data array of reference file
    for i in range(0, refgroups - 1):
        dark.data[0, i, :, :] = i * 0.1

    # apply correction
    outfile = darkcorr(dm_ramp, dark)

    # check that no correction/subtraction was applied; input file = output file
    diff = dm_ramp.data[:, :, :, :] - outfile.data[:, :, :, :]

    # test that the science data are not changed

    np.testing.assert_array_equal(np.full((nints, ngroups, ysize, xsize), 0.0, dtype=float),
                                  diff, err_msg='no changes should be seen in array ')


def test_sub_by_frame():
    '''Check that if NFRAMES=1 and GROUPGAP=0 for the science data, the dark reference data are
    directly subtracted frame by frame'''

    # size of integration
    nints = 1
    ngroups = 30
    xsize = 1032
    ysize = 1024

    # create raw input data for step
    dm_ramp = make_rampmodel(nints, ngroups, ysize, xsize)
    dm_ramp.meta.exposure.nframes = 1
    dm_ramp.meta.exposure.groupgap = 0

    # populate data array of science cube
    for i in range(0, ngroups-1):
        dm_ramp.data[0, i, :, :] = i

    refgroups = 50
    # create dark reference file model with fewer frames than science data
    dark = make_darkmodel(refgroups, ysize, xsize)

    # populate data array of reference file
    for i in range(0, refgroups - 1):
        dark.data[0, i, :, :] = i * 0.1

    # apply correction
    outfile = darkcorr(dm_ramp, dark)

    # remove the single dimension at start of file (1, 30, 1032, 1024) so comparison in assert works
    outdata = np.squeeze(outfile.data)

    # check that the dark file is subtracted frame by frame from the science data
    diff = dm_ramp.data[0, :, :, :] - dark.data[0, :ngroups, :, :]

    # test that the output data file is equal to the difference found when subtracting ref file from sci file

    np.testing.assert_array_equal(outdata, diff, err_msg='dark file should be subtracted from sci file ')


def test_nan():
    '''Verify that when a dark has NaNs, these are correctly assumed as zero and the PIXELDQ is set properly'''

    # size of integration
    nints = 1
    ngroups = 10
    xsize = 1032
    ysize = 1024

    # create raw input data for step
    dm_ramp = make_rampmodel(nints, ngroups, ysize, xsize)
    dm_ramp.meta.exposure.nframes = 1
    dm_ramp.meta.exposure.groupgap = 0

    # populate data array of science cube
    for i in range(0, ngroups-1):
        dm_ramp.data[0, i, :, :] = i

    refgroups = 15
    # create dark reference file model with fewer frames than science data
    dark = make_darkmodel(refgroups, ysize, xsize)

    # populate data array of reference file
    for i in range(0, refgroups - 1):
        dark.data[0, i, :, :] = i * 0.1

    # set NaN in dark file
    dark.data[0, 5, 500, 500] = np.nan

    # apply correction
    outfile = darkcorr(dm_ramp, dark)

    print(outfile.pixeldq[500, 500])
    print(outfile.groupdq[0, 5, 500, 500])

    # test that the NaN dark reference pixel was set to 0 (nothing subtracted)
    assert outfile.data[0, 5, 500, 500] == 5.0

    # test that the output dq file is flagged (with what)


def test_dq_combine():
    '''Verify that the DQ array of the dark is correctly combined with the PIXELDQ array of the science data.'''

    # size of integration
    nints = 1
    ngroups = 5
    xsize = 1032
    ysize = 1024

    # create raw input data for step
    dm_ramp = make_rampmodel(nints, ngroups, ysize, xsize)
    dm_ramp.meta.exposure.nframes = 1
    dm_ramp.meta.exposure.groupgap = 0

    # populate data array of science cube
    for i in range(0, ngroups-1):
        dm_ramp.data[0, i, :, :] = i

    refgroups = 15
    # create dark reference file model with fewer frames than science data
    dark = make_darkmodel(refgroups, ysize, xsize)

    # populate dq flags of sci pixeldq and reference dq
    dm_ramp.pixeldq[500, 500] = 4
    dm_ramp.pixeldq[500, 501] = 2

    dark.dq[0, 0, 500, 500] = 1
    dark.dq[0, 0, 500, 501] = 1

    # run correction step
    outfile = darkcorr(dm_ramp, dark)

    # check that dq flags were correctly added
    assert outfile.pixeldq[500, 500] == 5
    assert outfile.pixeldq[500, 501] == 3


def test_2_int():
    '''Verify the dark correction is done by integration for MIRI observations'''

    # size of integration
    nints = 2
    ngroups = 10
    xsize = 1032
    ysize = 1024

    # create raw input data for step
    dm_ramp = make_rampmodel(nints, ngroups, ysize, xsize)
    dm_ramp.meta.exposure.nframes = 1
    dm_ramp.meta.exposure.groupgap = 0

    # populate data array of science cube
    for i in range(0, ngroups-1):
        dm_ramp.data[:, i, :, :] = i

    refgroups = 15
    # create dark reference file model with fewer frames than science data
    dark = make_darkmodel(refgroups, ysize, xsize)

    # populate data array of reference file
    for i in range(0, refgroups - 1):
        dark.data[0, i, :, :] = i * 0.1
        dark.data[1, i, :, :] = i * 0.2

    # run correction
    outfile = darkcorr(dm_ramp, dark)

    # perform subtractions manually

    # check that the dark file is subtracted frame by frame from the science data
    diff = dm_ramp.data[0, :, :, :] - dark.data[0, :ngroups, :, :]
    diff_int2 = dm_ramp.data[1, :, :, :] - dark.data[1, :ngroups, :, :]

    # test that the output data file is equal to the difference found when subtracting ref file from sci file

    np.testing.assert_array_equal(outfile.data[0, :, :, :], diff,
                                  err_msg='dark file should be subtracted from sci file ')
    np.testing.assert_array_equal(outfile.data[1, :, :, :], diff_int2,
                                  err_msg='dark file should be subtracted from sci file ')


def test_dark_skipped():
    '''Verify that when the dark is not applied, the data is correctly flagged as such.'''

    # size of integration
    nints = 1
    ngroups = 30
    xsize = 1032
    ysize = 1024

    # create raw input data for step
    dm_ramp = make_rampmodel(nints, ngroups, ysize, xsize)
    dm_ramp.meta.exposure.nframes = 1
    dm_ramp.meta.exposure.groupgap = 0

    # populate data array of science cube
    for i in range(0, ngroups-1):
        dm_ramp.data[0, i, :, :] = i

    refgroups = 20
    # create dark reference file model with fewer frames than science data
    dark = make_darkmodel(refgroups, ysize, xsize)

    # populate data array of reference file
    for i in range(0, refgroups - 1):
        dark.data[0, i, :, :] = i * 0.1

    # apply correction
    outfile = darkcorr(dm_ramp, dark)

    # get dark correction status from header
    darkstatus = outfile.meta.cal_step.dark_sub
    print('Dark status', darkstatus)

    assert darkstatus == 'SKIPPED'


def test_frame_avg():
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
    for i in range(0, ngroups-1):
        dm_ramp.data[:, i, :, :] = i + 1

    refgroups = 20
    # create dark reference file model
    dark = make_darkmodel(refgroups, ysize, xsize)

    # populate data array of reference file
    for i in range(0, refgroups - 1):
        dark.data[0, i, :, :] = i * 0.1

    # apply correction
    outfile = darkcorr(dm_ramp, dark)

    # dark frames should be averaged in groups of 4 frames
    # this will result in average values of 0.15, 0.55, .095. 1.35

    assert outfile.data[0, 0, 500, 500] == pytest.approx(0.85)
    assert outfile.data[0, 1, 500, 500] == pytest.approx(1.45)
    assert outfile.data[0, 2, 500, 500] == pytest.approx(2.05)
    assert outfile.data[0, 3, 500, 500] == pytest.approx(2.65)

    # check that the error array is not modified.
    np.testing.assert_array_equal(outfile.err[:, :], 0,
                                  err_msg='error array should remain 0 ')


def test_sub_masklyot():
    '''make sure correct subarray (masklyot) reference file is called'''

    # size of integration
    nints = 1
    ngroups = 10
    xsize = 320
    ysize = 304

    # create raw input data for step
    dm_ramp = make_rampmodel(nints, ngroups, ysize, xsize)
    dm_ramp.meta.exposure.nframes = 1
    dm_ramp.meta.exposure.groupgap = 0
    dm_ramp.meta.exposure.readpatt = 'FAST'
    dm_ramp.meta.exposure.type = 'MIR_IMAGE'
    dm_ramp.meta.instrument.detector = 'MIRIMAGE'
    dm_ramp.meta.instrument.filter = 'F2300C'
    dm_ramp.meta.instrument.band = 'N/A'
    dm_ramp.meta.subarray.name = 'MASKLYOT'

    # populate data array of science cube
    for i in range(0, ngroups-1):
        dm_ramp.data[:, i, :, :] = i

    # set subarray
    sub = 'MASKLYOT'

    # run pipeline
    outfile = DarkCurrentStep.call(dm_ramp)

    # get dark reference file from header
    darkfile = outfile.meta.ref_file.dark.name
    print('Darkfile', darkfile)

    # parse darkfile name to discard crds://
    dfile = darkfile.split("/")[2]
    print('split file', dfile)

    # read in reference file
    filepath = '/grp/crds/jwst/references/jwst/' + dfile
    reffile = fits.open(filepath)
    refhead = reffile[0].header
    sub_file = refhead['SUBARRAY']
    refdata = reffile[1].data

    # do manual subtraction of a single frame
    diff = dm_ramp.data[0, 3, :, :] - refdata[0, 3, :, :]

    # check that correction equals manual subtraction
    np.testing.assert_array_equal(outfile.data[0, 3, :, :], diff,
                                  err_msg='dark file should be subtracted from sci file ')

    # test if filters match
    assert (sub == sub_file)


def test_sub_mask1550():
    '''make sure correct subarray (mask1550) reference file is called'''

    # size of integration
    nints = 1
    ngroups = 10
    xsize = 288
    ysize = 224

    # create raw input data for step
    dm_ramp = make_rampmodel(nints, ngroups, ysize, xsize)
    dm_ramp.meta.exposure.nframes = 1
    dm_ramp.meta.exposure.groupgap = 0
    dm_ramp.meta.exposure.readpatt = 'FAST'
    dm_ramp.meta.exposure.type = 'MIR_IMAGE'
    dm_ramp.meta.instrument.detector = 'MIRIMAGE'
    dm_ramp.meta.instrument.filter = 'F1550C'
    dm_ramp.meta.instrument.band = 'N/A'
    dm_ramp.meta.subarray.name = 'MASK1550'

    # populate data array of science cube
    for i in range(0, ngroups-1):
        dm_ramp.data[:, i, :, :] = i

    # set subarray
    sub = 'MASK1550'

    # run pipeline
    outfile = DarkCurrentStep.call(dm_ramp)

    # get dark reference file from header
    darkfile = outfile.meta.ref_file.dark.name
    print('Darkfile', darkfile)

    # parse darkfile name to discard crds://
    dfile = darkfile.split("/")[2]
    print('split file', dfile)

    # read in reference file
    filepath = '/grp/crds/jwst/references/jwst/' + dfile
    reffile = fits.open(filepath)
    refhead = reffile[0].header
    sub_file = refhead['SUBARRAY']

    refdata = reffile[1].data

    # do manual subtraction of a single frame
    diff = dm_ramp.data[0, 3, :, :] - refdata[0, 3, :, :]

    # check that correction equals manual subtraction
    np.testing.assert_array_equal(outfile.data[0, 3, :, :], diff,
                                  err_msg='dark file should be subtracted from sci file ')

    # test if filters match
    assert (sub == sub_file)


def test_sub_mask1140():
    '''make sure correct subarray (mask1140) reference file is called'''

    # size of integration
    nints = 1
    ngroups = 10
    xsize = 288
    ysize = 224

    # create raw input data for step
    dm_ramp = make_rampmodel(nints, ngroups, ysize, xsize)
    dm_ramp.meta.exposure.nframes = 1
    dm_ramp.meta.exposure.groupgap = 0
    dm_ramp.meta.exposure.readpatt = 'FAST'
    dm_ramp.meta.exposure.type = 'MIR_IMAGE'
    dm_ramp.meta.instrument.detector = 'MIRIMAGE'
    dm_ramp.meta.instrument.filter = 'F1140C'
    dm_ramp.meta.instrument.band = 'N/A'
    dm_ramp.meta.subarray.name = 'MASK1140'

    # populate data array of science cube
    for i in range(0, ngroups-1):
        dm_ramp.data[:, i, :, :] = i

    # set subarray
    sub = 'MASK1140'

    # run pipeline
    outfile = DarkCurrentStep.call(dm_ramp)

    # get dark reference file from header
    darkfile = outfile.meta.ref_file.dark.name
    print('Darkfile', darkfile)

    # parse darkfile name to discard crds://
    dfile = darkfile.split("/")[2]
    print('split file', dfile)

    # read in reference file
    filepath = '/grp/crds/jwst/references/jwst/' + dfile
    reffile = fits.open(filepath)
    refhead = reffile[0].header
    sub_file = refhead['SUBARRAY']

    refdata = reffile[1].data

    # do manual subtraction of a single frame
    diff = dm_ramp.data[0, 3, :, :] - refdata[0, 3, :, :]

    # check that correction equals manual subtraction
    np.testing.assert_array_equal(outfile.data[0, 3, :, :], diff,
                                  err_msg='dark file should be subtracted from sci file ')

    # test if filters match
    assert (sub == sub_file)


def test_sub_mask1065():
    '''make sure correct subarray (mask1065) reference file is called'''

    # size of integration
    nints = 1
    ngroups = 10
    xsize = 288
    ysize = 224

    # create raw input data for step
    dm_ramp = make_rampmodel(nints, ngroups, ysize, xsize)
    dm_ramp.meta.exposure.nframes = 1
    dm_ramp.meta.exposure.groupgap = 0
    dm_ramp.meta.exposure.readpatt = 'FAST'
    dm_ramp.meta.exposure.type = 'MIR_IMAGE'
    dm_ramp.meta.instrument.detector = 'MIRIMAGE'
    dm_ramp.meta.instrument.filter = 'F1065C'
    dm_ramp.meta.instrument.band = 'N/A'
    dm_ramp.meta.subarray.name = 'MASK1065'

    # populate data array of science cube
    for i in range(0, ngroups-1):
        dm_ramp.data[:, i, :, :] = i

    # set subarray
    sub = 'MASK1065'

    # run pipeline
    outfile = DarkCurrentStep.call(dm_ramp)

    # get dark reference file from header
    darkfile = outfile.meta.ref_file.dark.name
    print('Darkfile', darkfile)

    # parse darkfile name to discard crds://
    dfile = darkfile.split("/")[2]
    print('split file', dfile)

    # read in reference file
    filepath = '/grp/crds/jwst/references/jwst/' + dfile
    reffile = fits.open(filepath)
    refhead = reffile[0].header
    sub_file = refhead['SUBARRAY']

    refdata = reffile[1].data

    # do manual subtraction of a single frame
    diff = dm_ramp.data[0, 3, :, :] - refdata[0, 3, :, :]

    # check that correction equals manual subtraction
    np.testing.assert_array_equal(outfile.data[0, 3, :, :], diff,
                                  err_msg='dark file should be subtracted from sci file ')

    # test if filters match
    assert (sub == sub_file)


def test_sub_sub256():
    '''make sure correct subarray (sub256) reference file is called'''

    # size of integration
    nints = 1
    ngroups = 10
    xsize = 256
    ysize = 256

    # create raw input data for step
    dm_ramp = make_rampmodel(nints, ngroups, ysize, xsize)
    dm_ramp.meta.exposure.nframes = 1
    dm_ramp.meta.exposure.groupgap = 0
    dm_ramp.meta.exposure.readpatt = 'FAST'
    dm_ramp.meta.exposure.type = 'MIR_IMAGE'
    dm_ramp.meta.instrument.detector = 'MIRIMAGE'
    dm_ramp.meta.instrument.filter = 'F560W'
    dm_ramp.meta.instrument.band = 'N/A'
    dm_ramp.meta.subarray.name = 'SUB256'

    # populate data array of science cube
    for i in range(0, ngroups-1):
        dm_ramp.data[:, i, :, :] = i

    # set subarray
    sub = 'SUB256'

    # run pipeline
    outfile = DarkCurrentStep.call(dm_ramp)

    # get dark reference file from header
    darkfile = outfile.meta.ref_file.dark.name
    print('Darkfile', darkfile)

    # parse darkfile name to discard crds://
    dfile = darkfile.split("/")[2]
    print('split file', dfile)

    # read in reference file
    filepath = '/grp/crds/jwst/references/jwst/' + dfile
    reffile = fits.open(filepath)
    refhead = reffile[0].header
    sub_file = refhead['SUBARRAY']

    refdata = reffile[1].data

    # do manual subtraction of a single frame
    diff = dm_ramp.data[0, 3, :, :] - refdata[0, 3, :, :]

    # check that correction equals manual subtraction
    np.testing.assert_array_equal(outfile.data[0, 3, :, :], diff,
                                  err_msg='dark file should be subtracted from sci file ')

    # test if filters match
    assert (sub == sub_file)


def test_sub_sub128():
    '''make sure correct subarray (sub128) reference file is called'''

    # size of integration
    nints = 1
    ngroups = 10
    xsize = 136
    ysize = 128

    # create raw input data for step
    dm_ramp = make_rampmodel(nints, ngroups, ysize, xsize)
    dm_ramp.meta.exposure.nframes = 1
    dm_ramp.meta.exposure.groupgap = 0
    dm_ramp.meta.exposure.readpatt = 'FAST'
    dm_ramp.meta.exposure.type = 'MIR_IMAGE'
    dm_ramp.meta.instrument.detector = 'MIRIMAGE'
    dm_ramp.meta.instrument.filter = 'F770W'
    dm_ramp.meta.instrument.band = 'N/A'
    dm_ramp.meta.subarray.name = 'SUB128'

    # populate data array of science cube
    for i in range(0, ngroups-1):
        dm_ramp.data[:, i, :, :] = i

    # set subarray
    sub = 'SUB128'

    # run pipeline
    outfile = DarkCurrentStep.call(dm_ramp)

    # get dark reference file from header
    darkfile = outfile.meta.ref_file.dark.name
    print('Darkfile', darkfile)

    # parse darkfile name to discard crds://
    dfile = darkfile.split("/")[2]
    print('split file', dfile)

    # read in reference file
    filepath = '/grp/crds/jwst/references/jwst/' + dfile
    reffile = fits.open(filepath)
    refhead = reffile[0].header
    sub_file = refhead['SUBARRAY']

    refdata = reffile[1].data

    # do manual subtraction of a single frame
    diff = dm_ramp.data[0, 3, :, :] - refdata[0, 3, :, :]

    # check that correction equals manual subtraction
    np.testing.assert_array_equal(outfile.data[0, 3, :, :], diff,
                                  err_msg='dark file should be subtracted from sci file ')

    # test if filters match
    assert (sub == sub_file)


def test_sub_sub64():
    '''make sure correct subarray (sub64) reference file is called'''

    # size of integration
    nints = 1
    ngroups = 10
    xsize = 72
    ysize = 64

    # create raw input data for step
    dm_ramp = make_rampmodel(nints, ngroups, ysize, xsize)
    dm_ramp.meta.exposure.nframes = 1
    dm_ramp.meta.exposure.groupgap = 0
    dm_ramp.meta.exposure.readpatt = 'FAST'
    dm_ramp.meta.exposure.type = 'MIR_IMAGE'
    dm_ramp.meta.instrument.detector = 'MIRIMAGE'
    dm_ramp.meta.instrument.filter = 'F770W'
    dm_ramp.meta.instrument.band = 'N/A'
    dm_ramp.meta.subarray.name = 'SUB64'

    # populate data array of science cube
    for i in range(0, ngroups-1):
        dm_ramp.data[:, i, :, :] = i

    # set subarray
    sub = 'SUB64'

    # run pipeline
    outfile = DarkCurrentStep.call(dm_ramp)

    # get dark reference file from header
    darkfile = outfile.meta.ref_file.dark.name
    print('Darkfile', darkfile)

    # parse darkfile name to discard crds://
    dfile = darkfile.split("/")[2]
    print('split file', dfile)

    # read in reference file
    filepath = '/grp/crds/jwst/references/jwst/' + dfile
    reffile = fits.open(filepath)
    refhead = reffile[0].header
    sub_file = refhead['SUBARRAY']

    refdata = reffile[1].data

    # do manual subtraction of a single frame
    diff = dm_ramp.data[0, 3, :, :] - refdata[0, 3, :, :]

    # check that correction equals manual subtraction
    np.testing.assert_array_equal(outfile.data[0, 3, :, :], diff,
                                  err_msg='dark file should be subtracted from sci file ')

    # test if filters match
    assert (sub == sub_file)


def test_sub_brightsky():
    '''make sure correct subarray (brightsky) reference file is called'''

    # size of integration
    nints = 1
    ngroups = 10
    xsize = 512
    ysize = 512

    # create raw input data for step
    dm_ramp = make_rampmodel(nints, ngroups, ysize, xsize)
    dm_ramp.meta.exposure.nframes = 1
    dm_ramp.meta.exposure.groupgap = 0
    dm_ramp.meta.exposure.readpatt = 'FAST'
    dm_ramp.meta.exposure.type = 'MIR_IMAGE'
    dm_ramp.meta.instrument.detector = 'MIRIMAGE'
    dm_ramp.meta.instrument.filter = 'F770W'
    dm_ramp.meta.instrument.band = 'N/A'
    dm_ramp.meta.subarray.name = 'BRIGHTSKY'

    # populate data array of science cube
    for i in range(0, ngroups-1):
        dm_ramp.data[:, i, :, :] = i

    # set subarray
    sub = 'BRIGHTSKY'

    # run pipeline
    outfile = DarkCurrentStep.call(dm_ramp)

    # get dark reference file from header
    darkfile = outfile.meta.ref_file.dark.name
    print('Darkfile', darkfile)

    # parse darkfile name to discard crds://
    dfile = darkfile.split("/")[2]
    print('split file', dfile)

    # read in reference file
    filepath = '/grp/crds/jwst/references/jwst/' + dfile
    reffile = fits.open(filepath)
    refhead = reffile[0].header
    sub_file = refhead['SUBARRAY']

    refdata = reffile[1].data

    # do manual subtraction of a single frame
    diff = dm_ramp.data[0, 3, :, :] - refdata[0, 3, :, :]

    # check that correction equals manual subtraction
    np.testing.assert_array_equal(outfile.data[0, 3, :, :], diff,
                                  err_msg='dark file should be subtracted from sci file ')

    # test if filters match
    assert (sub == sub_file)


@pytest.mark.xfail
def test_sub_slitlessprism():
    '''make sure correct subarray reference file is called'''
    # we have issues with mismatch between SLITLESSPRISM and SUBPRISM in the headers of the files

    # size of integration
    nints = 1
    ngroups = 10
    xsize = 72
    ysize = 416

    # create raw input data for step
    dm_ramp = make_rampmodel(nints, ngroups, ysize, xsize)
    dm_ramp.meta.exposure.nframes = 1
    dm_ramp.meta.exposure.groupgap = 0
    dm_ramp.meta.exposure.readpatt = 'FAST'
    dm_ramp.meta.exposure.type = 'MIR_IMAGE'
    dm_ramp.meta.instrument.detector = 'MIRIMAGE'
    dm_ramp.meta.instrument.filter = 'F770W'
    dm_ramp.meta.instrument.band = 'N/A'
    dm_ramp.meta.subarray.name = 'SUBPRISM'

    # populate data array of science cube
    for i in range(0, ngroups-1):
        dm_ramp.data[:, i, :, :] = i

    # set subarray
    sub = 'SUBPRISM'

    # run pipeline
    outfile = DarkCurrentStep.call(dm_ramp)

    # get dark reference file from header
    darkfile = outfile.meta.ref_file.dark.name
    print('Darkfile', darkfile)

    # parse darkfile name to discard crds://
    dfile = darkfile.split("/")[2]
    print('split file', dfile)

    # read in reference file
    filepath = '/grp/crds/jwst/references/jwst/' + dfile
    reffile = fits.open(filepath)
    refhead = reffile[0].header
    sub_file = refhead['SUBARRAY']

    refdata = reffile[1].data

    # do manual subtraction of a single frame
    diff = dm_ramp.data[0, 3, :, :] - refdata[0, 3, :, :]

    # check that correction equals manual subtraction
    np.testing.assert_array_equal(outfile.data[0, 3, :, :], diff,
                                  err_msg='dark file should be subtracted from sci file ')

    print(sub, sub_file)
    # test if filters match
    assert (sub == sub_file)  # tests show sub=SUBPRISM and sub_file=SLITLESSPRISM


def test_slowmode():
    '''make sure correct slow mode reference file is called'''

    # size of integration
    nints = 1
    ngroups = 10
    xsize = 1032
    ysize = 1024

    # create raw input data for step
    dm_ramp = make_rampmodel(nints, ngroups, ysize, xsize)
    dm_ramp.meta.exposure.nframes = 1
    dm_ramp.meta.exposure.groupgap = 0
    dm_ramp.meta.exposure.readpatt = 'SLOW'
    dm_ramp.meta.exposure.type = 'MIR_IMAGE'
    dm_ramp.meta.instrument.detector = 'MIRIMAGE'
    dm_ramp.meta.instrument.filter = 'F770W'
    dm_ramp.meta.instrument.band = 'N/A'
    dm_ramp.meta.subarray.name = 'FULL'

    # populate data array of science cube
    for i in range(0, ngroups-1):
        dm_ramp.data[:, i, :, :] = i

    # set mode
    readout = 'SLOW'

    # run pipeline
    outfile = DarkCurrentStep.call(dm_ramp)

    # get dark reference file from header
    darkfile = outfile.meta.ref_file.dark.name
    print('Darkfile', darkfile)

    # parse darkfile name to discard crds://
    dfile = darkfile.split("/")[2]
    print('split file', dfile)

    # read in reference file
    filepath = '/grp/crds/jwst/references/jwst/' + dfile
    reffile = fits.open(filepath)
    refhead = reffile[0].header
    rmode = refhead['READPATT']
    refdata = reffile[1].data

    # do manual subtraction of a single frame
    diff = dm_ramp.data[0, 3, :, :] - refdata[0, 3, :, :]

    # check that correction equals manual subtraction
    np.testing.assert_array_equal(outfile.data[0, 3, :, :], diff,
                                  err_msg='dark file should be subtracted from sci file ')

    # test if filters match
    assert (readout == rmode)


def test_ifulong():
    '''make sure correct mrs (IFULONG) mode reference file is called'''

    # size of integration
    nints = 1
    ngroups = 10
    xsize = 1032
    ysize = 1024

    # create raw input data for step
    dm_ramp = make_rampmodel(nints, ngroups, ysize, xsize)
    dm_ramp.meta.exposure.nframes = 1
    dm_ramp.meta.exposure.groupgap = 0
    dm_ramp.meta.exposure.readpatt = 'FAST'
    dm_ramp.meta.exposure.type = 'MIR_MRS'
    dm_ramp.meta.instrument.detector = 'MIRIFULONG'
    dm_ramp.meta.instrument.filter = 'F770W'
    dm_ramp.meta.instrument.band = 'N/A'
    dm_ramp.meta.subarray.name = 'FULL'

    # populate data array of science cube
    for i in range(0, ngroups-1):
        dm_ramp.data[:, i, :, :] = i

    # set mode
    det = 'MIRIFULONG'

    # run pipeline
    outfile = DarkCurrentStep.call(dm_ramp)

    # get dark reference file from header
    darkfile = outfile.meta.ref_file.dark.name
    print('Darkfile', darkfile)

    # parse darkfile name to discard crds://
    dfile = darkfile.split("/")[2]
    print('split file', dfile)

    # read in reference file
    filepath = '/grp/crds/jwst/references/jwst/' + dfile
    reffile = fits.open(filepath)
    refhead = reffile[0].header
    det_file = refhead['DETECTOR']
    refdata = reffile[1].data

    # do manual subtraction of a single frame
    diff = dm_ramp.data[0, 3, :, :] - refdata[0, 3, :, :]

    # check that correction equals manual subtraction
    np.testing.assert_array_equal(outfile.data[0, 3, :, :], diff,
                                  err_msg='dark file should be subtracted from sci file ')

    # test if filters match
    assert (det == det_file)


def test_ifushort():
    '''make sure correct mrs (IFUSHORT) mode reference file is called'''

    # size of integration
    nints = 1
    ngroups = 10
    xsize = 1032
    ysize = 1024

    # create raw input data for step
    dm_ramp = make_rampmodel(nints, ngroups, ysize, xsize)
    dm_ramp.meta.exposure.nframes = 1
    dm_ramp.meta.exposure.groupgap = 0
    dm_ramp.meta.exposure.readpatt = 'FAST'
    dm_ramp.meta.exposure.type = 'MIR_MRS'
    dm_ramp.meta.instrument.detector = 'MIRIFUSHORT'
    dm_ramp.meta.instrument.filter = 'F770W'
    dm_ramp.meta.instrument.band = 'N/A'
    dm_ramp.meta.subarray.name = 'FULL'

    # populate data array of science cube
    for i in range(0, ngroups-1):
        dm_ramp.data[:, i, :, :] = i

    # set mode
    det = 'MIRIFUSHORT'

    # run pipeline
    outfile = DarkCurrentStep.call(dm_ramp)

    # get dark reference file from header
    darkfile = outfile.meta.ref_file.dark.name
    print('Darkfile', darkfile)

    # parse darkfile name to discard crds://
    dfile = darkfile.split("/")[2]
    print('split file', dfile)

    # read in reference file
    filepath = '/grp/crds/jwst/references/jwst/' + dfile
    reffile = fits.open(filepath)
    refhead = reffile[0].header
    det_file = refhead['DETECTOR']
    refdata = reffile[1].data

    # do manual subtraction of a single frame
    diff = dm_ramp.data[0, 3, :, :] - refdata[0, 3, :, :]

    # check that correction equals manual subtraction
    np.testing.assert_array_equal(outfile.data[0, 3, :, :], diff,
                                  err_msg='dark file should be subtracted from sci file ')

    # test if filters match
    assert (det == det_file)


def make_rampmodel(nints, ngroups, ysize, xsize):
    '''Make MIRI Ramp model for testing'''
    # create the data and groupdq arrays
    csize = (nints, ngroups, ysize, xsize)
    data = np.full(csize, 1.0)
    pixeldq = np.zeros((ysize, xsize), dtype=int)
    groupdq = np.zeros(csize, dtype=int)
    err = np.zeros((ysize, xsize), dtype=int)

    # create a JWST datamodel for MIRI data
    dm_ramp = MIRIRampModel(data=data, pixeldq=pixeldq, groupdq=groupdq, err=err)

    dm_ramp.meta.instrument.name = 'MIRI'
    dm_ramp.meta.observation.date = '2018-01-01'
    dm_ramp.meta.observation.time = '00:00:00'
    dm_ramp.meta.subarray.xstart = 1
    dm_ramp.meta.subarray.xsize = xsize
    dm_ramp.meta.subarray.ystart = 1
    dm_ramp.meta.subarray.ysize = ysize
    dm_ramp.meta.description = 'Fake data.'


    return dm_ramp


def make_darkmodel(ngroups, ysize, xsize):
    '''Make MIRI dark model for testing'''
    # create the data and groupdq arrays
    nints = 2
    csize = (nints, ngroups, ysize, xsize)
    data = np.full(csize, 1.0)
    dq = np.zeros((nints, 1, ysize, xsize), dtype=int)

    # create a JWST datamodel for MIRI data
    dark = DarkMIRIModel(data=data, dq=dq)

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


@pytest.fixture(scope='function')
def setup_nrc_cube():
    '''Set up fake NIRCam data to test.'''

    def _cube(readpatt, ngroups, nrows, ncols):

        nints = 1
        groupgap = READPATTERNS[readpatt.upper()]['nskip']
        nframes = READPATTERNS[readpatt.upper()]['nframes']

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

        dark_model = DarkModel((100, 2048, 2048))
        dark_model.meta.subarray.xstart = 1
        dark_model.meta.subarray.ystart = 1
        dark_model.meta.subarray.xsize = 2048
        dark_model.meta.subarray.ysize = 2048
        dark_model.meta.exposure.ngroups = 100
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