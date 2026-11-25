"""
Unit tests for dark current correction
"""

import warnings

import numpy as np
import pytest
from numpy.testing import assert_allclose

from jwst.dark_current.dark_class import DarkData, ScienceData
from jwst.dark_current.dark_sub import average_dark_frames_3d as average_dark_frames
from jwst.dark_current.dark_sub import do_correction_data as darkcorr

dqflags = {
    "DO_NOT_USE": 2**0,  # Bad pixel. Do not use.
    "SATURATED": 2**1,  # Pixel saturated during exposure
    "JUMP_DET": 2**2,  # Jump detected during exposure
    "DROPOUT": 2**3,  # Data lost in transmission
    "AD_FLOOR": 2**6,  # Below A/D floor (0 DN, was RESERVED_3)
}


# Define frame_time and number of groups in the generated dark reffile
TFRAME = 10.73677
NGROUPS_DARK = 10

DELIM = "-" * 80


@pytest.fixture()
def make_rampmodel():
    """Make MIRI Ramp model for testing"""

    def _ramp(nints, ngroups, nrows, ncols):
        ramp_data = ScienceData()

        dims = (nints, ngroups, nrows, ncols)
        imshape = (nrows, ncols)

        ramp_data.data = np.full(dims, 1.0, dtype=np.float32)
        ramp_data.groupdq = np.zeros(dims, dtype=np.uint32)
        ramp_data.pixeldq = np.zeros(imshape, dtype=np.uint32)

        ramp_data.instrument_name = "MIRI"

        ramp_data.cal_step = None

        return ramp_data

    return _ramp


@pytest.fixture()
def make_darkmodel():
    """Make MIRI dark model for testing"""

    def _dark(ngroups, nrows, ncols):
        # create the data and groupdq arrays
        nints = 2
        csize = (nints, ngroups, nrows, ncols)

        dark = DarkData()
        dark.data = np.full(csize, 1.0)
        dark.groupdq = np.zeros(csize, dtype=np.uint32)

        dark.exp_nframes = 1
        dark.exp_ngroups = ngroups
        dark.exp_groupgap = 0

        return dark

    return _dark


@pytest.fixture()
def setup_nrc_cube():
    """Set up fake NIRCam data to test."""

    def _cube(readpatt, ngroups, nframes, groupgap, nrows, ncols):
        nints = 1
        dims = (nints, ngroups, nrows, ncols)
        ramp_data = ScienceData()
        ramp_data.data = np.zeros(dims, dtype=np.float32)
        ramp_data.groupdq = np.zeros(dims, dtype=np.uint32)
        ramp_data.pixeldq = np.zeros(dims[-2:], dtype=np.uint32)

        ramp_data.instrument_name = "NIRCAM"
        ramp_data.exp_nframes = nframes
        ramp_data.exp_groupgap = groupgap

        ramp_data.cal_step = None

        dark_dims = (NGROUPS_DARK, 2048, 2048)
        dark_data = DarkData(dark_dims)
        dark_data.exp_nframes = nframes
        dark_data.exp_ngroups = NGROUPS_DARK
        dark_data.exp_groupgap = 0

        return ramp_data, dark_data

    return _cube


def _params():
    """Returns list of tuples, one for each readpatt, generating parameters for
    test_frame_averaging. Parameters are the following:

        (readpatt, ngroups, nframes, groupgap, nrows, ncols)

    Note groupgap = nskip
    """

    # Dictionary of NIRCam readout patterns
    readpatterns = {
        "DEEP8": {"ngroups": 20, "nframes": 8, "nskip": 12},
        "DEEP2": {"ngroups": 20, "nframes": 2, "nskip": 18},
        "MEDIUM8": {"ngroups": 10, "nframes": 8, "nskip": 2},
        "MEDIUM2": {"ngroups": 10, "nframes": 2, "nskip": 8},
        "SHALLOW4": {"ngroups": 10, "nframes": 4, "nskip": 1},
        "SHALLOW2": {"ngroups": 10, "nframes": 2, "nskip": 3},
        "BRIGHT2": {"ngroups": 10, "nframes": 2, "nskip": 0},
        "BRIGHT1": {"ngroups": 10, "nframes": 1, "nskip": 1},
        "RAPID": {"ngroups": 10, "nframes": 1, "nskip": 0},
    }

    params = []
    ngroups = 3
    # NIRCam is 2048x2048, but we reduce the ncols to 20x20 for speed/memory
    nrows = 20
    ncols = 20
    for readpatt, values in readpatterns.items():
        params.append((readpatt, ngroups, values["nframes"], values["nskip"], nrows, ncols))

    return params


@pytest.mark.parametrize(
    ("readpatt", "ngroups", "nframes", "groupgap", "nrows", "ncols"), _params()
)
def test_frame_averaging(setup_nrc_cube, readpatt, ngroups, nframes, groupgap, nrows, ncols):
    """Check that if nframes>1 or groupgap>0, then the pipeline reconstructs
    the dark reference file to match the frame averaging and groupgap
    settings of the exposure."""

    # Create data and dark model
    data, dark = setup_nrc_cube(readpatt, ngroups, nframes, groupgap, nrows, ncols)

    # Add ramp values to dark model data array
    dark.data[:, 10, 10] = np.arange(0, NGROUPS_DARK)

    # Run the pipeline's averaging function
    avg_dark = average_dark_frames(dark, ngroups, nframes, groupgap)

    # Group input groups into collections of frames which will be averaged
    total_frames = (nframes * ngroups) + (groupgap * (ngroups - 1))

    # Get starting/ending indexes of the input groups to be averaged
    gstrt_ind = np.arange(0, total_frames, nframes + groupgap)
    gend_ind = gstrt_ind + nframes

    # Prepare arrays to hold results of averaging
    manual_avg = np.zeros((ngroups), dtype=np.float32)

    # Manually average the input data to compare with pipeline output
    with warnings.catch_warnings(), np.errstate(invalid="ignore"):
        warnings.filterwarnings("ignore", message="Mean of empty slice", category=RuntimeWarning)
        for newgp, gstart, gend in zip(range(ngroups), gstrt_ind, gend_ind, strict=True):
            # Average the data frames
            newframe = np.mean(dark.data[gstart:gend, 10, 10])
            manual_avg[newgp] = newframe

    # Check that pipeline output matches manual averaging results
    assert_allclose(manual_avg, avg_dark.data[:, 10, 10], rtol=1e-5)

    # Check that meta data was properly updated
    assert avg_dark.exp_nframes == nframes
    assert avg_dark.exp_ngroups == ngroups
    assert avg_dark.exp_groupgap == groupgap


def test_sub_by_frame(make_rampmodel, make_darkmodel):
    """
    Check that if NFRAMES=1 and GROUPGAP=0 for the science data,
    the dark reference data are directly subtracted frame by frame
    """

    # size of integration
    nints, ngroups, nrows, ncols = 1, 10, 200, 200

    # create raw input data for step
    dm_ramp = make_rampmodel(nints, ngroups, nrows, ncols)
    dm_ramp.exp_nframes = 1
    dm_ramp.exp_groupgap = 0

    # populate data array of science cube
    for i in range(ngroups - 1):
        dm_ramp.data[0, i] = i

    # create dark reference file model with more frames than science data
    refgroups = 15
    dark = make_darkmodel(refgroups, nrows, ncols)

    # populate data array of reference file
    for i in range(refgroups - 1):
        dark.data[0, i] = i * 0.1

    # apply correction
    outfile, avg_dark = darkcorr(dm_ramp, dark)

    assert outfile.cal_step == "COMPLETE"

    # remove the single dimension at start of file (1, 30, 1032, 1024)
    # so comparison in assert works
    outdata = np.squeeze(outfile.data)

    # check that the dark file is subtracted frame by frame
    # from the science data
    diff = dm_ramp.data[0] - dark.data[0, :ngroups]

    # test that the output data file is equal to the difference
    # found when subtracting ref file from sci file
    tol = 1.0e-6
    np.testing.assert_allclose(outdata, diff, tol)


def test_nan(make_rampmodel, make_darkmodel):
    """
    Verify that when a dark has NaNs, these are correctly
    assumed as zero and the PIXELDQ is set properly
    """

    # size of integration
    nints, ngroups, nrows, ncols = 1, 10, 200, 200

    # create raw input data for step
    dm_ramp = make_rampmodel(nints, ngroups, nrows, ncols)
    dm_ramp.exp_nframes = 1
    dm_ramp.exp_groupgap = 0

    # populate data array of science cube
    for i in range(ngroups - 1):
        dm_ramp.data[0, i, :, :] = i

    # create dark reference file model with more frames than science data
    refgroups = 15
    dark = make_darkmodel(refgroups, nrows, ncols)

    # populate data array of reference file
    for i in range(refgroups - 1):
        dark.data[0, i] = i * 0.1

    # set NaN in dark file

    dark.data[0, 5, 100, 100] = np.nan

    # apply correction
    outfile, avg_dark = darkcorr(dm_ramp, dark)

    # test that the NaN dark reference pixel was set to 0 (nothing subtracted)
    assert outfile.data[0, 5, 100, 100] == 5.0


def test_dq_combine(make_rampmodel, make_darkmodel):
    """
    Verify that the DQ array of the dark is correctly
    combined with the PIXELDQ array of the science data.
    """

    # size of integration
    nints, ngroups, nrows, ncols = 1, 5, 200, 200

    # create raw input data for step
    dm_ramp = make_rampmodel(nints, ngroups, nrows, ncols)
    dm_ramp.exp_nframes = 1
    dm_ramp.exp_groupgap = 0

    # populate data array of science cube
    for i in range(1, ngroups - 1):
        dm_ramp.data[0, i, :, :] = i

    # create dark reference file model with more frames than science data
    refgroups = 7
    dark = make_darkmodel(refgroups, nrows, ncols)

    # populate dq flags of sci pixeldq and reference dq
    dm_ramp.pixeldq[50, 50] = dqflags["JUMP_DET"]
    dm_ramp.pixeldq[50, 51] = dqflags["SATURATED"]

    dark.groupdq[0, 0, 50, 50] = dqflags["DO_NOT_USE"]
    dark.groupdq[0, 0, 50, 51] = dqflags["DO_NOT_USE"]

    # run correction step
    outfile, avg_dark = darkcorr(dm_ramp, dark)

    # check that dq flags were correctly added
    assert outfile.pixeldq[50, 50] == np.bitwise_or(dqflags["JUMP_DET"], dqflags["DO_NOT_USE"])

    assert outfile.pixeldq[50, 51] == np.bitwise_or(dqflags["SATURATED"], dqflags["DO_NOT_USE"])


def test_frame_avg(make_rampmodel, make_darkmodel):
    """
    Check that if NFRAMES>1 or GROUPGAP>0, the frame-averaged dark data are
    subtracted group-by-group from science data groups.
    """

    # size of integration
    nints, ngroups, nrows, ncols = 1, 5, 1024, 1032

    # create raw input data for step
    dm_ramp = make_rampmodel(nints, ngroups, nrows, ncols)
    dm_ramp.exp_nframes = 4
    dm_ramp.exp_groupgap = 0

    # populate data array of science cube
    for i in range(ngroups - 1):
        dm_ramp.data[:, i] = i + 1

    # create dark reference file model

    refgroups = 20  # This needs to be 20 groups for the calculations to work
    dark = make_darkmodel(refgroups, nrows, ncols)

    # populate data array of reference file
    for i in range(refgroups - 1):
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


def test_dark_extrapolation(make_rampmodel, make_darkmodel, setup_nrc_cube):
    """
    Check that the dark is extrapolated when it has insufficient frames to cover the science input.

    MIRI uses multi-integration 4-D darks, while NIR instruments use 3-D single-int darks.
    Extrapolation code branches depending on dark shape, so test both.
    """

    # size of integration
    nints, ngroups, nrows, ncols = 2, 20, 1024, 1032

    # create raw input data for step
    dm_ramp = make_rampmodel(nints, ngroups, nrows, ncols)
    dm_ramp.exp_nframes = 1
    dm_ramp.exp_groupgap = 0

    # Science array will have rate of 1, starting at 1.
    for i in range(ngroups):
        dm_ramp.data[:, i] = i + 1

    # create dark reference file model

    refgroups = 10  # This needs to be <20 groups for the extrapolation to occur.
    dark = make_darkmodel(refgroups, nrows, ncols)

    # Int 1 will have dark current of 0.1, starting at 0.
    # Int 2 will have dark current of 0.3, starting at 0.
    for i in range(refgroups):
        dark.data[0, i] = i * 0.1
        dark.data[1, i] = i * 0.3
    # apply correction
    outfile, avg_dark = darkcorr(dm_ramp, dark)

    assert_allclose(outfile.data[0, :, 500, 500], np.linspace(1.0, 18.1, ngroups), rtol=1.0e-5)
    assert_allclose(dark.data[0, :, 500, 500], np.linspace(0.0, 1.9, ngroups), rtol=1.0e-5)
    assert_allclose(outfile.data[1, :, 500, 500], np.linspace(1.0, 14.3, ngroups), rtol=1.0e-5)
    assert_allclose(dark.data[1, :, 500, 500], np.linspace(0.0, 5.7, ngroups), rtol=1.0e-5)

    nrc_ngroups = 40
    data, dark = setup_nrc_cube("rp", nrc_ngroups, nframes=1, groupgap=0, nrows=2048, ncols=2048)

    # Add ramp values to dark model data array
    dark.data[:, 10, 10] = np.arange(0, NGROUPS_DARK) * 0.2

    data.data[:, :, 10, 10] = np.arange(1, nrc_ngroups + 1)

    outfile, avg_dark = darkcorr(data, dark)

    assert_allclose(outfile.data[0, :, 10, 10], np.linspace(1, 32.2, nrc_ngroups), rtol=1.0e-5)
    assert_allclose(dark.data[:, 10, 10], np.linspace(0, 7.8, nrc_ngroups), rtol=1.0e-5)
