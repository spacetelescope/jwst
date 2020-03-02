import numpy as np
import pytest

from jwst.datamodels import GainModel, ReadnoiseModel
from jwst.datamodels import RampModel
from jwst.jump.jump import detect_jumps
import multiprocessing
from jwst.datamodels import dqflags

def test_nocrs_noflux(setup_inputs):
    """"
    All pixel values are zero. So slope should be zero
    """
    model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=5)
    out_model = detect_jumps(model1, gain, rnModel, 4.0,  1, 200, 4, True)
    assert (0 == np.max(out_model.groupdq))

def test_nocrs_noflux_badgain_pixel(setup_inputs):
    """"
    all pixel values are zero. So slope should be zero, pixel with bad gain should
    have pixel dq set to 'NO_GAIN_VALUE' and 'DO_NOT_USE'
    """
    model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=5)
    gain.data[7, 7] = -10 #bad gain
    gain.data[17, 17] = np.nan  # bad gain
    out_model = detect_jumps(model1, gain, rnModel, 4.0,  1, 200, 4, True)
    assert(np.bitwise_and(out_model.pixeldq[7, 7], dqflags.pixel['NO_GAIN_VALUE']))
    assert (np.bitwise_and(out_model.pixeldq[7, 7], dqflags.pixel['DO_NOT_USE']))
    assert (np.bitwise_and(out_model.pixeldq[17, 17], dqflags.pixel['NO_GAIN_VALUE']))
    assert (np.bitwise_and(out_model.pixeldq[17, 17], dqflags.pixel['DO_NOT_USE']))


def test_nocrs_noflux_subarray(setup_inputs):
    """"
    All pixel values are zero. This shows that the subarray reference files get extracted from the full frame
    versions.
    """
    model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=5, subarray=True)
    out_model = detect_jumps(model1, gain, rnModel, 4.0,  1, 200, 4, True)
    assert (0 == np.max(out_model.groupdq))

def test_onecr_10_groups_neighbors_flagged(setup_inputs):
    """"
    A single CR in a 10 group exposure
    """
    grouptime = 3.0
    ingain = 200
    inreadnoise = np.float64(7)
    ngroups = 10
    model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=ngroups,
                                                          gain=ingain, readnoise=inreadnoise, deltatime=grouptime)
    # two segments perfect fit, second segment has twice the slope
    model1.data[0, 0, 5, 5] = 15.0
    model1.data[0, 1, 5, 5] = 20.0
    model1.data[0, 2, 5, 5] = 25.0
    model1.data[0, 3, 5, 5] = 30.0
    model1.data[0, 4, 5, 5] = 35.0
    model1.data[0, 5, 5, 5] = 140.0
    model1.data[0, 6, 5, 5] = 150.0
    model1.data[0, 7, 5, 5] = 160.0
    model1.data[0, 8, 5, 5] = 170.0
    model1.data[0, 9, 5, 5] = 180.0
    out_model = detect_jumps(model1, gain, rnModel, 4.0,  1, 200, 4, True)
    assert (4 == np.max(out_model.groupdq[0, 5, 5, 5]))
    assert (4 == out_model.groupdq[0, 5, 5, 6])
    assert (4 == out_model.groupdq[0, 5, 5, 4])
    assert (4 == out_model.groupdq[0, 5, 6, 5])
    assert (4 == out_model.groupdq[0, 5, 4, 5])

def test_nocr_100_groups_nframes1(setup_inputs):
    """"
    NO CR in a 100 group exposure to make sure that frames_per_group is passed correctly to
    twopoint_difference. This test recreates the problem found in issue #4571.
    """
    grouptime = 3.0
    ingain = 1 #to make the noise calculation simple
    inreadnoise = np.float64(7)
    ngroups = 100
    model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=ngroups, nrows=100, ncols=100,
                                                          gain=ingain, readnoise=inreadnoise,
                                                          deltatime=grouptime)
    model1.meta.exposure.nframes = 1
    # two segments perfect fit, second segment has twice the slope
    model1.data[0, 0, 5, 5] = 14.0
    model1.data[0, 1, 5, 5] = 20.0
    model1.data[0, 2, 5, 5] = 27.0
    model1.data[0, 3, 5, 5] = 30.0
    model1.data[0, 4, 5, 5] = 38.0
    model1.data[0, 5, 5, 5] = 40.0
    model1.data[0, 6, 5, 5] = 50.0
    model1.data[0, 7, 5, 5] = 52.0
    model1.data[0, 8, 5, 5] = 63.0
    model1.data[0, 9, 5, 5] = 68.0
    for i in range(10,100):
        model1.data[0,i,5,5] = i * 5
    out_model = detect_jumps(model1, gain, rnModel, 4.0,  1, 200, 4, True)
    assert (0 == np.max(out_model.groupdq))

def test_twoints_onecr_each_10_groups_neighbors_flagged(setup_inputs):
    """"
    Two integrations with CRs in different locations. This makes sure we are correctly
    dealing with integrations.
    """
    grouptime = 3.0
    ingain = 200
    inreadnoise = np.float64(7)
    ngroups = 10
    model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=ngroups, nints=2,
                                                          gain=ingain, readnoise=inreadnoise, deltatime=grouptime)
    # two segments perfect fit, second segment has twice the slope
    model1.data[0, 0, 5, 5] = 15.0
    model1.data[0, 1, 5, 5] = 20.0
    model1.data[0, 2, 5, 5] = 25.0
    model1.data[0, 3, 5, 5] = 30.0
    model1.data[0, 4, 5, 5] = 35.0
    model1.data[0, 5, 5, 5] = 140.0
    model1.data[0, 6, 5, 5] = 150.0
    model1.data[0, 7, 5, 5] = 160.0
    model1.data[0, 8, 5, 5] = 170.0
    model1.data[0, 9, 5, 5] = 180.0
    model1.data[1, 0, 15, 5] = 15.0
    model1.data[1, 1, 15, 5] = 20.0
    model1.data[1, 2, 15, 5] = 25.0
    model1.data[1, 3, 15, 5] = 30.0
    model1.data[1, 4, 15, 5] = 35.0
    model1.data[1, 5, 15, 5] = 40.0
    model1.data[1, 6, 15, 5] = 45.0
    model1.data[1, 7, 15, 5] = 160.0
    model1.data[1, 8, 15, 5] = 170.0
    model1.data[1, 9, 15, 5] = 180.0
    out_model = detect_jumps(model1, gain, rnModel, 4.0,  1, 200, 4, True)
    assert (4 == np.max(out_model.groupdq[0, 5, 5, 5]))
    assert (4 == out_model.groupdq[0, 5, 5, 6])
    assert (4 == out_model.groupdq[0, 5, 5, 4])
    assert (4 == out_model.groupdq[0, 5, 6, 5])
    assert (4 == out_model.groupdq[0, 5, 4, 5])
    assert (4 == out_model.groupdq[1, 7, 15, 5])
    assert (4 == out_model.groupdq[1, 7, 15, 6])
    assert (4 == out_model.groupdq[1, 7, 15, 4])
    assert (4 == out_model.groupdq[1, 7, 16, 5])
    assert (4 == out_model.groupdq[1, 7, 14, 5])

def test_flagging_of_CRs_across_slice_boundaries(setup_inputs):
    """"
    A multiprocessing test that has two CRs on the boundary between two slices.
    This makes sure that we are correctly flagging neighbors in different  slices.
    """
    grouptime = 3.0
    ingain = 200
    inreadnoise = np.float64(7)
    ngroups = 10

    model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=ngroups, nints=2,
                                                          gain=ingain, readnoise=inreadnoise,
                                                          deltatime=grouptime)
    nrows = model1.data.shape[3]
    num_cores = multiprocessing.cpu_count()
    max_cores = 'half'
    numslices = num_cores // 2
    if numslices > 1:
        yincrement = int(nrows / numslices)
        # two segments perfect fit, second segment has twice the slope
        #add a CR on the last row of the first slice
        model1.data[0, 0, yincrement-1, 5] = 15.0
        model1.data[0, 1, yincrement-1, 5] = 20.0
        model1.data[0, 2, yincrement-1, 5] = 25.0
        model1.data[0, 3, yincrement-1, 5] = 30.0
        model1.data[0, 4, yincrement-1, 5] = 35.0
        model1.data[0, 5, yincrement-1, 5] = 140.0
        model1.data[0, 6, yincrement-1, 5] = 150.0
        model1.data[0, 7, yincrement-1, 5] = 160.0
        model1.data[0, 8, yincrement-1, 5] = 170.0
        model1.data[0, 9, yincrement-1, 5] = 180.0
        #add a CR on the first row of the second slice
        model1.data[1, 0, yincrement, 25] = 15.0
        model1.data[1, 1, yincrement, 25] = 20.0
        model1.data[1, 2, yincrement, 25] = 25.0
        model1.data[1, 3, yincrement, 25] = 30.0
        model1.data[1, 4, yincrement, 25] = 35.0
        model1.data[1, 5, yincrement, 25] = 40.0
        model1.data[1, 6, yincrement, 25] = 50.0
        model1.data[1, 7, yincrement, 25] = 160.0
        model1.data[1, 8, yincrement, 25] = 170.0
        model1.data[1, 9, yincrement, 25] = 180.0
        out_model = detect_jumps(model1, gain, rnModel, 4.0,  max_cores, 200, 4, True)
        #check that the neighbors of the CR on the last row were flagged
        assert (4 == out_model.groupdq[0, 5, yincrement-1, 5])
        assert (4 == out_model.groupdq[0, 5, yincrement-1, 6])
        assert (4 == out_model.groupdq[0, 5, yincrement-1, 4])
        assert (4 == out_model.groupdq[0, 5, yincrement, 5])
        assert (4 == out_model.groupdq[0, 5, yincrement-2, 5])
        # check that the neighbors of the CR on the first row were flagged
        assert (4 == out_model.groupdq[1, 7, yincrement, 25])
        assert (4 == out_model.groupdq[1, 7, yincrement, 26])
        assert (4 == out_model.groupdq[1, 7, yincrement, 24])
        assert (4 == out_model.groupdq[1, 7, yincrement+1, 25])
        assert (4 == out_model.groupdq[1, 7, yincrement-1, 25])


def test_twoints_onecr_10_groups_neighbors_flagged_multi(setup_inputs):
    """"
    A multiprocessing test that has two CRs on the boundary between two slices
    in different integrations. This makes sure that we are correctly flagging
    neighbors in different slices and that we are parsing the integrations correctly.
    """
    grouptime = 3.0
    ingain = 200
    inreadnoise = np.float64(7)
    ngroups = 10
    model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=ngroups, nints=2,
                                                          gain=ingain, readnoise=inreadnoise, deltatime=grouptime)
    # two segments perfect fit, second segment has twice the slope
    model1.data[0, 0, 5, 5] = 15.0
    model1.data[0, 1, 5, 5] = 20.0
    model1.data[0, 2, 5, 5] = 25.0
    model1.data[0, 3, 5, 5] = 30.0
    model1.data[0, 4, 5, 5] = 35.0
    model1.data[0, 5, 5, 5] = 140.0
    model1.data[0, 6, 5, 5] = 150.0
    model1.data[0, 7, 5, 5] = 160.0
    model1.data[0, 8, 5, 5] = 170.0
    model1.data[0, 9, 5, 5] = 180.0
    model1.data[1, 0, 15, 5] = 15.0
    model1.data[1, 1, 15, 5] = 20.0
    model1.data[1, 2, 15, 5] = 25.0
    model1.data[1, 3, 15, 5] = 30.0
    model1.data[1, 4, 15, 5] = 35.0
    model1.data[1, 5, 15, 5] = 40.0
    model1.data[1, 6, 15, 5] = 45.0
    model1.data[1, 7, 15, 5] = 160.0
    model1.data[1, 8, 15, 5] = 170.0
    model1.data[1, 9, 15, 5] = 180.0
    out_model = detect_jumps(model1, gain, rnModel, 4.0,  'half', 200, 4, True)
    assert (4 == np.max(out_model.groupdq[0, 5, 5, 5]))
    assert (4 == out_model.groupdq[0, 5, 5, 6])
    assert (4 == out_model.groupdq[0, 5, 5, 4])
    assert (4 == out_model.groupdq[0, 5, 6, 5])
    assert (4 == out_model.groupdq[0, 5, 4, 5])
    assert (4 == out_model.groupdq[1, 7, 15, 5])
    assert (4 == out_model.groupdq[1, 7, 15, 6])
    assert (4 == out_model.groupdq[1, 7, 15, 4])
    assert (4 == out_model.groupdq[1, 7, 16, 5])
    assert (4 == out_model.groupdq[1, 7, 14, 5])

@pytest.mark.skip(reason="Test is only used to test performance issue. No need to run every time.")
def test_every_pixel_CR_neighbors_flagged(setup_inputs):
    """"
    A multiprocessing test that has a jump in every pixel. This is used
    to test the performance gain from multiprocessing.
    """
    grouptime = 3.0
    ingain = 200
    inreadnoise = np.float64(7)
    ngroups = 10
    model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=ngroups,
                                                          gain=ingain, readnoise=inreadnoise, deltatime=grouptime)
    # two segments perfect fit, second segment has twice the slope
    model1.data[0, 0, :, :] = 15.0
    model1.data[0, 1, :, :] = 20.0
    model1.data[0, 2, :, :] = 25.0
    model1.data[0, 3, :, :] = 30.0
    model1.data[0, 4, :, :] = 35.0
    model1.data[0, 5, :, :] = 140.0
    model1.data[0, 6, :, :] = 150.0
    model1.data[0, 7, :, :] = 160.0
    model1.data[0, 8, :, :] = 170.0
    model1.data[0, 9, :, :] = 180.0
    out_model = detect_jumps(model1, gain, rnModel, 4.0,  'half', 200, 4, True)
    assert (4 == np.max(out_model.groupdq[0, 5, 5, 5]))
    assert (4 == out_model.groupdq[0, 5, 5, 6])
    assert (4 == out_model.groupdq[0, 5, 5, 4])
    assert (4 == out_model.groupdq[0, 5, 6, 5])
    assert (4 == out_model.groupdq[0, 5, 4, 5])

def test_crs_on_edge_with_neighbor_flagging(setup_inputs):
    """"
    A test to make sure that the neighbors of CRs on the edges of the
    array are flagged correctly.
    """
    grouptime = 3.0
    ingain = 200
    inreadnoise = np.float64(7)
    ngroups = 10
    model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=ngroups,
                                                          gain=ingain, readnoise=inreadnoise,
                                                          deltatime=grouptime)
    # two segments perfect fit, second segment has twice the slope
    # CR on 1st row
    model1.data[0, 0, 0, 15] = 15.0
    model1.data[0, 1, 0, 15] = 20.0
    model1.data[0, 2, 0, 15] = 25.0
    model1.data[0, 3, 0, 15] = 30.0
    model1.data[0, 4, 0, 15] = 35.0
    model1.data[0, 5, 0, 15] = 140.0
    model1.data[0, 6, 0, 15] = 150.0
    model1.data[0, 7, 0, 15] = 160.0
    model1.data[0, 8, 0, 15] = 170.0
    model1.data[0, 9, 0, 15] = 180.0
    # CR on last row
    model1.data[0, 0, 1023, 5] = 15.0
    model1.data[0, 1, 1023, 5] = 20.0
    model1.data[0, 2, 1023, 5] = 25.0
    model1.data[0, 3, 1023, 5] = 30.0
    model1.data[0, 4, 1023, 5] = 35.0
    model1.data[0, 5, 1023, 5] = 140.0
    model1.data[0, 6, 1023, 5] = 150.0
    model1.data[0, 7, 1023, 5] = 160.0
    model1.data[0, 8, 1023, 5] = 170.0
    model1.data[0, 9, 1023, 5] = 180.0
    # CR on 1st column
    model1.data[0, 0, 5, 0] = 15.0
    model1.data[0, 1, 5, 0] = 20.0
    model1.data[0, 2, 5, 0] = 25.0
    model1.data[0, 3, 5, 0] = 30.0
    model1.data[0, 4, 5, 0] = 35.0
    model1.data[0, 5, 5, 0] = 140.0
    model1.data[0, 6, 5, 0] = 150.0
    model1.data[0, 7, 5, 0] = 160.0
    model1.data[0, 8, 5, 0] = 170.0
    model1.data[0, 9, 5, 0] = 180.0
    # CR on last column
    model1.data[0, 0, 15, 1027] = 15.0
    model1.data[0, 1, 15, 1027] = 20.0
    model1.data[0, 2, 15, 1027] = 25.0
    model1.data[0, 3, 15, 1027] = 30.0
    model1.data[0, 4, 15, 1027] = 35.0
    model1.data[0, 5, 15, 1027] = 140.0
    model1.data[0, 6, 15, 1027] = 150.0
    model1.data[0, 7, 15, 1027] = 160.0
    model1.data[0, 8, 15, 1027] = 170.0
    model1.data[0, 9, 15, 1027] = 180.0
    out_model = detect_jumps(model1, gain, rnModel, 4.0,  1, 200, 10, True)
    # flag CR and three neighbors of first row CR
    assert (4 == out_model.groupdq[0, 5, 0, 15])
    assert (4 == out_model.groupdq[0, 5, 1, 15])
    assert (4 == out_model.groupdq[0, 5, 0, 14])
    assert (4 == out_model.groupdq[0, 5, 0, 16])
    assert (out_model.groupdq[0, 5, -1, 15] == 0) # The one not to flag
    # flag CR and three neighbors of last row CR
    assert (4 == out_model.groupdq[0, 5, 1023, 5])
    assert (4 == out_model.groupdq[0, 5, 1022, 5])
    assert (4 == out_model.groupdq[0, 5, 1023, 4])
    assert (4 == out_model.groupdq[0, 5, 1023, 6])
    # flag CR and three neighbors of first column CR
    assert (4 == out_model.groupdq[0, 5, 5, 0])
    assert (4 == out_model.groupdq[0, 5, 6, 0])
    assert (4 == out_model.groupdq[0, 5, 4, 0])
    assert (4 == out_model.groupdq[0, 5, 5, 1])
    assert (out_model.groupdq[0, 5, 5, -1] == 0)# The one not to flag
    # flag CR and three neighbors of last column CR
    assert (4 == out_model.groupdq[0, 5, 15, 1027])
    assert (4 == out_model.groupdq[0, 5, 15, 1026])
    assert (4 == out_model.groupdq[0, 5, 16, 1027])
    assert (4 == out_model.groupdq[0, 5, 14, 1027])


def test_onecr_10_groups(setup_inputs):
    """"
    A test to make sure that neighbors are not flagged when they are not requested to be flagged.
    """
    grouptime = 3.0
    ingain = 200
    inreadnoise = np.float64(7)
    ngroups = 10
    model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=ngroups,
                                                          gain=ingain, readnoise=inreadnoise, deltatime=grouptime)
    # two segments perfect fit, second segment has twice the slope
    model1.data[0, 0, 5, 5] = 15.0
    model1.data[0, 1, 5, 5] = 20.0
    model1.data[0, 2, 5, 5] = 25.0
    model1.data[0, 3, 5, 5] = 30.0
    model1.data[0, 4, 5, 5] = 35.0
    model1.data[0, 5, 5, 5] = 140.0
    model1.data[0, 6, 5, 5] = 150.0
    model1.data[0, 7, 5, 5] = 160.0
    model1.data[0, 8, 5, 5] = 170.0
    model1.data[0, 9, 5, 5] = 180.0
    out_model = detect_jumps(model1, gain, rnModel, 4.0,  1, 200, 10, False)
    assert (out_model.groupdq[0, 5, 5, 5] == 4)
    assert (out_model.groupdq[0, 5, 4, 5] == 0)
    assert (out_model.groupdq[0, 5, 6, 5] == 0)
    assert (out_model.groupdq[0, 5, 5, 6] == 0)
    assert (out_model.groupdq[0, 5, 5, 4] == 0)

def test_onecr_10_groups_fullarray(setup_inputs):
    """"
    A test that has a cosmic ray in the 5th group for all pixels except column 10. In column
    10 the jump is in the 7th group.
    """
    grouptime = 3.0
    ingain = 5
    inreadnoise = np.float64(7)
    ngroups = 10
    model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=ngroups,
                                                          gain=ingain, readnoise=inreadnoise, deltatime=grouptime)
    #
    model1.data[0, 0, 5, :] = 15.0
    model1.data[0, 1, 5, :] = 20.0
    model1.data[0, 2, 5, :] = 25.0
    model1.data[0, 3, 5, :] = 30.0
    model1.data[0, 4, 5, :] = 35.0
    model1.data[0, 5, 5, :] = 140.0
    model1.data[0, 6, 5, :] = 150.0
    model1.data[0, 7, 5, :] = 160.0
    model1.data[0, 8, 5, :] = 170.0
    model1.data[0, 9, 5, :] = 180.0
    # move the CR to group 7 for row 10 and make difference be 300
    model1.data[0, 3, 5, 10] = 100
    model1.data[0, 4, 5, 10] = 130
    model1.data[0, 5, 5, 10] = 160
    model1.data[0, 6, 5, 10] = 190
    model1.data[0, 7, 5, 10] = 400
    model1.data[0, 8, 5, 10] = 410
    model1.data[0, 9, 5, 10] = 420
    out_model = detect_jumps(model1, gain, rnModel, 4.0,  1, 200, 10, False)
    assert (np.all(out_model.groupdq[0, 5, 5, 0:10] == 4)) # The jump is in group 5 for columns 0-9
    assert (out_model.groupdq[0, 7, 5, 10] == 4)  # The jump is in group 7 for column 10
    assert (np.all(out_model.groupdq[0, 5, 5, 11:] == 4)) # The jump is in group 5 for columns 11+


def test_onecr_50_groups(setup_inputs):
    """"
    A test with a fifty group integration. There are two jumps in pixel 5,5. One in group 5 and
    one in group 30.
    """
    grouptime = 3.0
    ingain = 5
    inreadnoise = np.float64(7)
    ngroups = 50
    model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=ngroups,
                                                          gain=ingain, readnoise=inreadnoise, deltatime=grouptime)

    model1.data[0, 0, 5, 5] = 15.0
    model1.data[0, 1, 5, 5] = 20.0
    model1.data[0, 2, 5, 5] = 25.0
    model1.data[0, 3, 5, 5] = 30.0
    model1.data[0, 4, 5, 5] = 35.0
    model1.data[0, 5, 5, 5] = 140.0
    model1.data[0, 6, 5, 5] = 150.0
    model1.data[0, 7, 5, 5] = 160.0
    model1.data[0, 8, 5, 5] = 170.0
    model1.data[0, 9, 5, 5] = 180.0
    model1.data[0, 10:30, 5, 5] = np.arange(190, 290, 5)
    model1.data[0, 30:50, 5, 5] = np.arange(500, 600, 5)
    out_model = detect_jumps(model1, gain, rnModel, 4.0,  1, 200, 10, False)
    assert (out_model.groupdq[0, 5, 5, 5] == 4) # CR in group 5
    assert (out_model.groupdq[0, 30, 5, 5] == 4) # CR in group 30
    assert (np.all(out_model.groupdq[0, 6:30, 5, 5] == 0)) # groups in between are not flagged


def test_single_CR_neighbor_flag( setup_inputs):
    """"
    A single CR in a 10 group exposure. Tests that:
    - if neighbor-flagging is set, the 4 neighboring pixels *ARE* flagged, and
    - if neighbor-flagging is *NOT* set, the 4 neighboring pixels are *NOT* flagged
    """
    grouptime = 3.0
    ingain = 5
    inreadnoise = np.float64(7)
    ngroups = 10

    model1, gdq, rnModel, pixdq, err, gain = \
        setup_inputs( ngroups=ngroups, nrows=5, ncols=6, gain=ingain, readnoise=inreadnoise,
                      deltatime=grouptime )

    # two segments perfect fit, second segment has twice the slope
    model1.data[0, 0, 3, 3] = 15.0
    model1.data[0, 1, 3, 3] = 20.0
    model1.data[0, 2, 3, 3] = 25.0
    model1.data[0, 3, 3, 3] = 30.0
    model1.data[0, 4, 3, 3] = 35.0
    model1.data[0, 5, 3, 3] = 140.0
    model1.data[0, 6, 3, 3] = 150.0
    model1.data[0, 7, 3, 3] = 160.0
    model1.data[0, 8, 3, 3] = 170.0
    model1.data[0, 9, 3, 3] = 180.0

    # Flag neighbors
    out_model = detect_jumps( model1, gain, rnModel, 4.0,  1, 200, 4, True )

    assert (4 == np.max(out_model.groupdq[0, 5, 3, 3]))
    assert (4 == out_model.groupdq[0, 5, 3, 4])
    assert (4 == out_model.groupdq[0, 5, 3, 2])
    assert (4 == out_model.groupdq[0, 5, 2, 3])
    assert (4 == out_model.groupdq[0, 5, 4, 3])

    # Do not flag neighbors
    out_model = detect_jumps( model1, gain, rnModel, 4.0,  1, 200, 4, False )

    assert (4 == np.max(out_model.groupdq[0, 5, 3, 3]))
    assert (0 == out_model.groupdq[0, 5, 3, 4])
    assert (0 == out_model.groupdq[0, 5, 3, 2])
    assert (0 == out_model.groupdq[0, 5, 2, 3])
    assert (0 == out_model.groupdq[0, 5, 4, 3])


def test_proc(setup_inputs):
    """"
    A single CR in a 10 group exposure. Verify that the pixels flagged using
    multiprocessing are identical to the pixels flagged when no
    multiprocessing is done.
    """
    grouptime = 3.0
    ingain = 5
    inreadnoise = np.float64(7)
    ngroups = 10

    model1, gdq, rnModel, pixdq, err, gain = \
        setup_inputs( ngroups=ngroups, nrows=5, ncols=6, nints=2, gain=ingain, readnoise=inreadnoise,
                      deltatime=grouptime )

    model1.data[0, 0, 2, 3] = 15.0
    model1.data[0, 1, 2, 3] = 21.0
    model1.data[0, 2, 2, 3] = 25.0
    model1.data[0, 3, 2, 3] = 30.2
    model1.data[0, 4, 2, 3] = 35.0
    model1.data[0, 5, 2, 3] = 140.0
    model1.data[0, 6, 2, 3] = 151.0
    model1.data[0, 7, 2, 3] = 160.0
    model1.data[0, 8, 2, 3] = 170.0
    model1.data[0, 9, 2, 3] = 180.0

    out_model_a = detect_jumps( model1, gain, rnModel, 4.0, 'half', 200, 4, True )
    out_model_b = detect_jumps( model1, gain, rnModel, 4.0, None, 200, 4, True )
    assert( out_model_a.groupdq == out_model_b.groupdq ).all()

    out_model_c = detect_jumps( model1, gain, rnModel, 4.0, 'All', 200, 4, True )
    assert( out_model_a.groupdq == out_model_c.groupdq ).all()


def test_adjacent_CRs( setup_inputs ):
    """
    Three CRs in a 10 group exposure; the CRs have overlapping neighboring
    pixels. This test makes sure that the correct pixels are flagged.
    """
    grouptime = 3.0
    ingain = 5
    inreadnoise = np.float64(7)
    ngroups = 10
    model1, gdq, rnModel, pixdq, err, gain = \
        setup_inputs( ngroups=ngroups, nrows=5, ncols=6, gain=ingain,
                      readnoise=inreadnoise, deltatime=grouptime )

    # Populate arrays for 1st CR, centered at (x=2, y=3)
    x=2; y=3
    model1.data[0, 0, y, x] = 15.0
    model1.data[0, 1, y, x] = 20.0
    model1.data[0, 2, y, x] = 26.0
    model1.data[0, 3, y, x] = 30.0
    model1.data[0, 4, y, x] = 35.0
    model1.data[0, 5, y, x] = 140.0
    model1.data[0, 6, y, x] = 150.0
    model1.data[0, 7, y, x] = 161.0
    model1.data[0, 8, y, x] = 170.0
    model1.data[0, 9, y, x] = 180.0

    # Populate arrays for 2nd CR, centered at (x=2, y=2)
    x=2; y=2
    model1.data[0, 0, y, x] = 20.0
    model1.data[0, 1, y, x] = 30.0
    model1.data[0, 2, y, x] = 41.0
    model1.data[0, 3, y, x] = 51.0
    model1.data[0, 4, y, x] = 62.0
    model1.data[0, 5, y, x] = 170.0
    model1.data[0, 6, y, x] = 200.0
    model1.data[0, 7, y, x] = 231.0
    model1.data[0, 8, y, x] = 260.0
    model1.data[0, 9, y, x] = 290.0

    # Populate arrays for 3rd CR, centered at (x=3, y=2)
    x=3; y=2
    model1.data[0, 0, y, x] = 120.0
    model1.data[0, 1, y, x] = 140.0
    model1.data[0, 2, y, x] = 161.0
    model1.data[0, 3, y, x] = 181.0
    model1.data[0, 4, y, x] = 202.0
    model1.data[0, 5, y, x] = 70.0
    model1.data[0, 6, y, x] = 100.0
    model1.data[0, 7, y, x] = 131.0
    model1.data[0, 8, y, x] = 160.0
    model1.data[0, 9, y, x] = 190.0

    out_model = detect_jumps(model1, gain, rnModel, 4.0, 'half', 200, 4, True)

    # 1st CR (centered at x=2, y=3)
    assert (4 == out_model.groupdq[ 0,5,2,2 ])
    assert (4 == out_model.groupdq[ 0,5,3,1 ])
    assert (4 == out_model.groupdq[ 0,5,3,2 ])
    assert (4 == out_model.groupdq[ 0,5,3,3 ])
    assert (4 == out_model.groupdq[ 0,5,4,2 ])

    # 2nd CR (centered at x=2, y=2)
    assert (4 == out_model.groupdq[ 0,5,1,2 ])
    assert (4 == out_model.groupdq[ 0,5,2,1 ])
    assert (4 == out_model.groupdq[ 0,5,2,3 ])

    # 3rd CR (centered at x=3, y=2)
    assert (4 == out_model.groupdq[ 0,5,1,3 ])
    assert (4 == out_model.groupdq[ 0,5,2,4 ])

# Need test for multi-ints near zero with positive and negative slopes

@pytest.fixture
def setup_inputs():
    def _setup(ngroups=10, readnoise=10, nints=1,
               nrows=1024, ncols=1032, nframes=1, grouptime=1.0, gain=1, deltatime=1,
               gain_subarray = False, readnoise_subarray = False, subarray = False):
        times = np.array(list(range(ngroups)), dtype=np.float64) * deltatime
        gain = np.ones(shape=(nrows, ncols), dtype=np.float64) * gain


        pixdq = np.zeros(shape=(nrows, ncols), dtype=np.float64)
        read_noise = np.full((nrows, ncols), readnoise, dtype=np.float64)
        gdq = np.zeros(shape=(nints, ngroups, nrows, ncols), dtype=np.uint32)
        if subarray:
            data = np.zeros(shape=(nints, ngroups, 20, 20), dtype=np.float64)
            err = np.ones(shape=(nints, ngroups, 20, 20), dtype=np.float64)
        else:
            data = np.zeros(shape=(nints, ngroups, nrows, ncols), dtype=np.float64)
            err = np.ones(shape=(nints, ngroups, nrows, ncols), dtype=np.float64)
        model1 = RampModel(data=data, err=err, pixeldq=pixdq, groupdq=gdq, times=times)
        model1.meta.instrument.name = 'MIRI'
        model1.meta.instrument.detector = 'MIRIMAGE'
        model1.meta.instrument.filter = 'F480M'
        model1.meta.observation.date = '2015-10-13'
        model1.meta.exposure.type = 'MIR_IMAGE'
        model1.meta.exposure.group_time = deltatime
        model1.meta.subarray.name = 'FULL'
        model1.meta.subarray.xstart = 1
        model1.meta.subarray.ystart = 1
        if subarray:
            model1.meta.subarray.xsize = 20
            model1.meta.subarray.ysize = 20
        else:
            model1.meta.subarray.xsize = ncols
            model1.meta.subarray.ysize = nrows
        model1.meta.exposure.frame_time = deltatime
        model1.meta.exposure.ngroups = ngroups
        model1.meta.exposure.group_time = deltatime
        model1.meta.exposure.nframes = 1
        model1.meta.exposure.groupgap = 0
        gain = GainModel(data=gain)
        gain.meta.instrument.name = 'MIRI'
        gain.meta.subarray.xstart = 1
        gain.meta.subarray.ystart = 1
        gain.meta.subarray.xsize = ncols
        gain.meta.subarray.ysize = nrows
        rnModel = ReadnoiseModel(data=read_noise)
        rnModel.meta.instrument.name = 'MIRI'
        rnModel.meta.subarray.xstart = 1
        rnModel.meta.subarray.ystart = 1
        rnModel.meta.subarray.xsize = ncols
        rnModel.meta.subarray.ysize = nrows
        return model1, gdq, rnModel, pixdq, err, gain

    return _setup
