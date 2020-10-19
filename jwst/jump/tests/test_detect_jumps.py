import numpy as np
from numpy.testing import assert_array_equal
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
    model1, rnModel, gain = setup_inputs(ngroups=5)
    out_model = detect_jumps(model1, gain, rnModel, 4.0,  1, 200, 4, True)
    assert (0 == np.max(out_model.groupdq))


def test_nocrs_noflux_badgain_pixel(setup_inputs):
    """"
    all pixel values are zero. So slope should be zero, pixel with bad gain should
    have pixel dq set to 'NO_GAIN_VALUE' and 'DO_NOT_USE'
    """
    model1, rnModel, gain = setup_inputs(ngroups=5, nrows=20, ncols=20)
    gain.data[7, 7] = -10  # bad gain
    gain.data[17, 17] = np.nan  # bad gain
    out_model = detect_jumps(model1, gain, rnModel, 4.0,  1, 200, 4, True)
    assert(np.bitwise_and(out_model.pixeldq[7, 7], dqflags.pixel['NO_GAIN_VALUE']))
    assert (np.bitwise_and(out_model.pixeldq[7, 7], dqflags.pixel['DO_NOT_USE']))
    assert (np.bitwise_and(out_model.pixeldq[17, 17], dqflags.pixel['NO_GAIN_VALUE']))
    assert (np.bitwise_and(out_model.pixeldq[17, 17], dqflags.pixel['DO_NOT_USE']))


def test_nocrs_noflux_subarray(setup_inputs):
    """"
    All pixel values are zero. This shows that the subarray reference files get
    extracted from the full frame versions.
    """
    model1, rnModel, gain = setup_inputs(ngroups=5, subarray=True)
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
    model1, rnModel, gain = setup_inputs(ngroups=ngroups, gain=ingain, nrows=10, ncols=10,
                                         readnoise=inreadnoise, deltatime=grouptime)

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
    ingain = 1  # to make the noise calculation simple
    inreadnoise = np.float64(7)
    ngroups = 100
    model1, rnModel, gain = setup_inputs(ngroups=ngroups, nrows=10, ncols=10,
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
    model1, rnModel, gain = setup_inputs(ngroups=ngroups, nints=2, nrows=20, ncols=20,
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


def test_multiple_neighbor_jumps_firstlastbad(setup_inputs):
    """
    This test is based on actual MIRI data that was having the incorrect
    group flagged with JUMP_DET (it was flagging group 2 instead of group 5).
    This makes sure that group 5 is getting flagged.
    Note that the first and last frames/groups are all flagged with DO_NOT_USE,
    due to the application of the first/last frame steps.
    """
    grouptime = 3.0
    ingain = 5.5
    inreadnoise = 6.5
    ngroups = 10
    nrows = 10
    ncols = 10

    model1, rnModel, gain = setup_inputs(
        ngroups=ngroups, nints=1, nrows=nrows, ncols=ncols,
        gain=ingain, readnoise=inreadnoise, deltatime=grouptime)

    # Setup the desired pixel values
    model1.data[0,:,1,1] = [10019.966, 10057.298, 10078.248, 10096.01,
       20241.627, 20248.752, 20268.047, 20284.895, 20298.705, 20314.25]
    model1.data[0,:,1,2] = [10016.457, 10053.907, 10063.568, 10076.166,
       11655.773, 11654.063, 11681.795, 11693.763, 11712.788, 11736.994]
    model1.data[0,:,1,3] = [10013.259, 10050.348, 10070.398, 10097.658, 10766.534, 10787.84 ,
       10802.418, 10818.872, 10832.695, 10861.175]
    model1.data[0,:,1,4] = [10016.422, 10053.959, 10070.934, 10090.381, 10104.014, 10127.665,
       10143.687, 10172.227, 10178.138, 10199.59]
    model1.data[0,:,2,1] = [10021.067, 10042.973, 10059.062, 10069.323, 18732.406, 18749.602,
       18771.908, 18794.695, 18803.223, 18819.523]
    model1.data[0,:,2,2] = [10019.651, 10043.371, 10056.423, 10085.121, 40584.703, 40606.08 ,
       40619.51 , 40629.574, 40641.9  , 40660.145]
    model1.data[0,:,2,3] = [10021.223, 10042.112, 10052.958, 10067.142, 28188.316, 28202.922,
       28225.557, 28243.79 , 28253.883, 28273.586]
    model1.data[0,:,2,4] = [10022.608, 10037.174, 10069.476, 10081.729, 11173.748, 11177.344,
       11201.127, 11219.607, 11229.468, 11243.174]
    model1.data[0,:,2,5] = [10011.095, 10047.422, 10061.066, 10079.375, 10106.405, 10116.071,
       10129.348, 10136.305, 10161.373, 10181.479]
    model1.data[0,:,3,1] = [10011.877, 10052.809, 10075.108, 10085.111, 10397.106, 10409.291,
       10430.475, 10445.3  , 10462.004, 10484.906]
    model1.data[0,:,3,2] = [10012.124 , 10059.202, 10078.984, 10092.74, 11939.488,
       11958.45, 11977.5625, 11991.776, 12025.897, 12027.326]
    model1.data[0,:,3,3] = [10013.282, 10046.887, 10062.308, 10085.447, 28308.426, 28318.957,
       28335.55 , 28353.832, 28371.746, 28388.848]
    model1.data[0,:,3,4] = [10016.784, 10048.249, 10060.097, 10074.606, 21506.082, 21522.027,
       21542.309, 21558.34 , 21576.365, 21595.58]
    model1.data[0,:,3,5] = [10014.916 , 10052.995 , 10063.7705, 10092.866, 10538.075,
       10558.318, 10570.754, 10597.343, 10608.488, 10628.104]
    model1.data[0,:,4,1] = [10017.438, 10038.94 , 10057.657, 10069.987, 10090.22 , 10114.296,
       10133.543, 10148.657, 10158.109, 10172.842]
    model1.data[0,:,4,2] = [10011.129, 10037.982, 10054.445, 10079.703, 10097.964, 10110.593,
       10135.701, 10149.448, 10171.771, 10185.874]
    model1.data[0,:,4,3] = [10021.109, 10043.658, 10063.909, 10072.364, 10766.232, 10774.402,
       10790.677, 10809.337, 10833.65 , 10849.55]
    model1.data[0,:,4,4] = [10023.877, 10035.997, 10052.321, 10077.937, 10529.645, 10541.947,
       10571.127, 10577.249, 10599.716, 10609.544]

    model1.groupdq[0,0,:,:] = 1  # Flag first frame as DO_NOT_USE
    model1.groupdq[0,-1,:,:] = 1  # Flag last frame as DO_NOT_USE

    # run jump detection
    out_model = detect_jumps(model1, gain, rnModel, rejection_threshold=200.0,
                             max_cores=None, max_jump_to_flag_neighbors=200,
                             min_jump_to_flag_neighbors=10, flag_4_neighbors=True)

    # Verify that the correct groups have been flagged. The entries for pixels
    # 2,2 and 3,3 are the ones that had previously been flagged in group 2 instead
    # of group 5.
    assert_array_equal(out_model.groupdq[0,:,1,1], [1, 0, 0, 0, 4, 0, 0, 0, 0, 1])
    assert_array_equal(out_model.groupdq[0,:,1,2], [1, 0, 0, 0, 4, 0, 0, 0, 0, 1])
    assert_array_equal(out_model.groupdq[0,:,1,3], [1, 0, 0, 0, 0, 0, 0, 0, 0, 1])
    assert_array_equal(out_model.groupdq[0,:,1,4], [1, 0, 0, 0, 0, 0, 0, 0, 0, 1])
    assert_array_equal(out_model.groupdq[0,:,2,1], [1, 0, 0, 0, 4, 0, 0, 0, 0, 1])
    assert_array_equal(out_model.groupdq[0,:,2,2], [1, 0, 0, 0, 4, 0, 0, 0, 0, 1])  # <--
    assert_array_equal(out_model.groupdq[0,:,2,3], [1, 0, 0, 0, 4, 0, 0, 0, 0, 1])
    assert_array_equal(out_model.groupdq[0,:,2,4], [1, 0, 0, 0, 0, 0, 0, 0, 0, 1])
    assert_array_equal(out_model.groupdq[0,:,3,1], [1, 0, 0, 0, 0, 0, 0, 0, 0, 1])
    assert_array_equal(out_model.groupdq[0,:,3,2], [1, 0, 0, 0, 4, 0, 0, 0, 0, 1])
    assert_array_equal(out_model.groupdq[0,:,3,3], [1, 0, 0, 0, 4, 0, 0, 0, 0, 1])  # <--
    assert_array_equal(out_model.groupdq[0,:,3,4], [1, 0, 0, 0, 4, 0, 0, 0, 0, 1])
    assert_array_equal(out_model.groupdq[0,:,4,1], [1, 0, 0, 0, 0, 0, 0, 0, 0, 1])
    assert_array_equal(out_model.groupdq[0,:,4,2], [1, 0, 0, 0, 0, 0, 0, 0, 0, 1])
    assert_array_equal(out_model.groupdq[0,:,4,3], [1, 0, 0, 0, 0, 0, 0, 0, 0, 1])
    assert_array_equal(out_model.groupdq[0,:,4,4], [1, 0, 0, 0, 0, 0, 0, 0, 0, 1])


def test_nirspec_saturated_pix(setup_inputs):
    """
    This test is based on an actual NIRSpec exposure that has some pixels
    flagged as saturated in one or more groups, which the jump step is
    supposed to ignore, but an old version of the code was setting JUMP flags
    for some of the saturated groups. This is to verify that the saturated
    groups are no longer flagged with jumps.
    """
    grouptime = 3.0
    ingain = 1.0
    inreadnoise = 10.7
    ngroups = 7
    nrows = 6
    ncols = 6

    model1, rnModel, gain = setup_inputs(ngroups=ngroups, nints=1, nrows=nrows, ncols=ncols,
                                         gain=ingain, readnoise=inreadnoise, deltatime=grouptime)

    # Setup the needed input pixel and DQ values
    model1.data[0,:,1,1] = [639854.75, 4872.451, -17861.791, 14022.15, 22320.176,
                            1116.3828, 1936.9746]
    model1.groupdq[0,:,1,1] = [0, 0, 0, 0, 0, 2, 2]
    model1.data[0,:,2,2] = [8.25666812e+05, -1.10471914e+05, 1.95755371e+02,  1.83118457e+03,
                            1.72250879e+03,  1.81733496e+03, 1.65188281e+03]
    model1.groupdq[0,:,2,2] = [0, 0, 2, 2, 2, 2, 2]
    model1.data[0,:,3,3] = [1228767., 46392.234, -3245.6553, 7762.413,
                            37190.76, 266611.62,  5072.4434]
    model1.groupdq[0,:,3,3] = [0, 0, 0, 0, 0, 0, 2]
    model1.data[0,:,4,4] = [7.5306038e+05, 1.8269953e+04, 1.8352356e+02, 2.1245061e+03,
                            2.0628525e+03, 2.1039399e+03, 2.0069873e+03]
    model1.groupdq[0,:,4,4] = [0, 0, 2, 2, 2, 2, 2]

    # run jump detection
    out_model = detect_jumps(model1, gain, rnModel, rejection_threshold=200.0,
                             max_cores=None, max_jump_to_flag_neighbors=200,
                             min_jump_to_flag_neighbors=10, flag_4_neighbors=True)

    # Check the results. There should not be any pixels with DQ values of 6, which
    # is saturated (2) plus jump (4). All the DQ's should be either just 2 or just 4.
    assert_array_equal(out_model.groupdq[0,:,1,1], [0, 4, 4, 4, 0, 2, 2])
    assert_array_equal(out_model.groupdq[0,:,2,2], [0, 4, 2, 2, 2, 2, 2])
    assert_array_equal(out_model.groupdq[0,:,3,3], [0, 4, 4, 0, 0, 4, 2])
    assert_array_equal(out_model.groupdq[0,:,4,4], [0, 4, 2, 2, 2, 2, 2])


def test_flagging_of_CRs_across_slice_boundaries(setup_inputs):
    """"
    A multiprocessing test that has two CRs on the boundary between two slices.
    This makes sure that we are correctly flagging neighbors in different  slices.
    """
    grouptime = 3.0
    ingain = 200
    inreadnoise = np.float64(7)
    ngroups = 10

    model1, rnModel, gain = setup_inputs(ngroups=ngroups, nints=2,
                                         gain=ingain, readnoise=inreadnoise,
                                         deltatime=grouptime)
    nrows = model1.data.shape[3]
    num_cores = multiprocessing.cpu_count()
    max_cores = 'half'
    numslices = num_cores // 2
    if numslices > 1:
        yincrement = int(nrows / numslices)
        # two segments perfect fit, second segment has twice the slope
        # add a CR on the last row of the first slice
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
        # add a CR on the first row of the second slice
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

        # check that the neighbors of the CR on the last row were flagged
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
    model1, rnModel, gain = setup_inputs(ngroups=ngroups, nints=2, nrows=40, ncols=10,
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
    model1, rnModel, gain = setup_inputs(ngroups=ngroups,
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
    model1, rnModel, gain = setup_inputs(ngroups=ngroups, nrows=20, ncols=20,
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
    model1.data[0, 0, 19, 5] = 15.0
    model1.data[0, 1, 19, 5] = 20.0
    model1.data[0, 2, 19, 5] = 25.0
    model1.data[0, 3, 19, 5] = 30.0
    model1.data[0, 4, 19, 5] = 35.0
    model1.data[0, 5, 19, 5] = 140.0
    model1.data[0, 6, 19, 5] = 150.0
    model1.data[0, 7, 19, 5] = 160.0
    model1.data[0, 8, 19, 5] = 170.0
    model1.data[0, 9, 19, 5] = 180.0
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
    model1.data[0, 0, 15, 19] = 15.0
    model1.data[0, 1, 15, 19] = 20.0
    model1.data[0, 2, 15, 19] = 25.0
    model1.data[0, 3, 15, 19] = 30.0
    model1.data[0, 4, 15, 19] = 35.0
    model1.data[0, 5, 15, 19] = 140.0
    model1.data[0, 6, 15, 19] = 150.0
    model1.data[0, 7, 15, 19] = 160.0
    model1.data[0, 8, 15, 19] = 170.0
    model1.data[0, 9, 15, 19] = 180.0

    out_model = detect_jumps(model1, gain, rnModel, 4.0,  1, 200, 10, True)

    # flag CR and three neighbors of first row CR
    assert (4 == out_model.groupdq[0, 5, 0, 15])
    assert (4 == out_model.groupdq[0, 5, 1, 15])
    assert (4 == out_model.groupdq[0, 5, 0, 14])
    assert (4 == out_model.groupdq[0, 5, 0, 16])
    assert (out_model.groupdq[0, 5, -1, 15] == 0) # The one not to flag
    # flag CR and three neighbors of last row CR
    assert (4 == out_model.groupdq[0, 5, 19, 5])
    assert (4 == out_model.groupdq[0, 5, 18, 5])
    assert (4 == out_model.groupdq[0, 5, 19, 4])
    assert (4 == out_model.groupdq[0, 5, 19, 6])
    # flag CR and three neighbors of first column CR
    assert (4 == out_model.groupdq[0, 5, 5, 0])
    assert (4 == out_model.groupdq[0, 5, 6, 0])
    assert (4 == out_model.groupdq[0, 5, 4, 0])
    assert (4 == out_model.groupdq[0, 5, 5, 1])
    assert (out_model.groupdq[0, 5, 5, -1] == 0)# The one not to flag
    # flag CR and three neighbors of last column CR
    assert (4 == out_model.groupdq[0, 5, 15, 19])
    assert (4 == out_model.groupdq[0, 5, 15, 18])
    assert (4 == out_model.groupdq[0, 5, 16, 19])
    assert (4 == out_model.groupdq[0, 5, 14, 19])


def test_onecr_10_groups(setup_inputs):
    """"
    A test to make sure that neighbors are not flagged when they are not requested to be flagged.
    """
    grouptime = 3.0
    ingain = 200
    inreadnoise = np.float64(7)
    ngroups = 10
    model1, rnModel, gain = setup_inputs(ngroups=ngroups, nrows=20, ncols=20,
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
    model1, rnModel, gain = setup_inputs(ngroups=ngroups, gain=ingain, nrows=20, ncols=20,
                                         readnoise=inreadnoise, deltatime=grouptime)

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

    assert (np.all(out_model.groupdq[0, 5, 5, 0:10] == 4))  # The jump is in group 5 for columns 0-9
    assert (out_model.groupdq[0, 7, 5, 10] == 4)  # The jump is in group 7 for column 10
    assert (np.all(out_model.groupdq[0, 5, 5, 11:] == 4))  # The jump is in group 5 for columns 11+


def test_onecr_50_groups(setup_inputs):
    """"
    A test with a fifty group integration. There are two jumps in pixel 5,5. One in group 5 and
    one in group 30.
    """
    grouptime = 3.0
    ingain = 5
    inreadnoise = np.float64(7)
    ngroups = 50
    model1, rnModel, gain = setup_inputs(ngroups=ngroups, nrows=10, ncols=10,
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

    assert (out_model.groupdq[0, 5, 5, 5] == 4)  # CR in group 5
    assert (out_model.groupdq[0, 30, 5, 5] == 4)  # CR in group 30
    assert (np.all(out_model.groupdq[0, 6:30, 5, 5] == 0))  # groups in between are not flagged


def test_single_CR_neighbor_flag(setup_inputs):
    """"
    A single CR in a 10 group exposure. Tests that:
    - if neighbor-flagging is set, the 4 neighboring pixels *ARE* flagged, and
    - if neighbor-flagging is *NOT* set, the 4 neighboring pixels are *NOT* flagged
    """
    grouptime = 3.0
    ingain = 5
    inreadnoise = np.float64(7)
    ngroups = 10

    model1, rnModel, gain = setup_inputs(ngroups=ngroups, nrows=5, ncols=6,
                                         gain=ingain, readnoise=inreadnoise,
                                         deltatime=grouptime)

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

    model1, rnModel, gain = setup_inputs(ngroups=ngroups, nrows=25, ncols=6,
                                         nints=2, gain=ingain, readnoise=inreadnoise,
                                         deltatime=grouptime)

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

    out_model_a = detect_jumps( model1, gain, rnModel, 4.0, None, 200, 4, True )
    out_model_b = detect_jumps( model1, gain, rnModel, 4.0, 'half', 200, 4, True )
    assert( out_model_a.groupdq == out_model_b.groupdq ).all()

    out_model_c = detect_jumps( model1, gain, rnModel, 4.0, 'all', 200, 4, True )
    assert( out_model_a.groupdq == out_model_c.groupdq ).all()


def test_adjacent_CRs(setup_inputs ):
    """
    Three CRs in a 10 group exposure; the CRs have overlapping neighboring
    pixels. This test makes sure that the correct pixels are flagged.
    """
    grouptime = 3.0
    ingain = 5
    inreadnoise = np.float64(7)
    ngroups = 10
    model1, rnModel, gain = setup_inputs(ngroups=ngroups, nrows=15, ncols=6, gain=ingain,
                                         readnoise=inreadnoise, deltatime=grouptime)

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
    assert (4 == out_model.groupdq[0, 5, 2, 2])
    assert (4 == out_model.groupdq[0, 5, 3, 1])
    assert (4 == out_model.groupdq[0, 5, 3, 2])
    assert (4 == out_model.groupdq[0, 5, 3, 3])
    assert (4 == out_model.groupdq[0, 5, 4, 2])

    # 2nd CR (centered at x=2, y=2)
    assert (4 == out_model.groupdq[0, 5, 1, 2])
    assert (4 == out_model.groupdq[0, 5, 2, 1])
    assert (4 == out_model.groupdq[0, 5, 2, 3])

    # 3rd CR (centered at x=3, y=2)
    assert (4 == out_model.groupdq[0, 5, 1, 3])
    assert (4 == out_model.groupdq[0, 5, 2, 4])

# Need test for multi-ints near zero with positive and negative slopes


@pytest.fixture
def setup_inputs():
    def _setup(ngroups=10, readnoise=10, nints=1,
               nrows=1024, ncols=1032, nframes=1, grouptime=1.0, gain=1, deltatime=1,
               gain_subarray = False, readnoise_subarray = False, subarray = False):

        # Populate data arrays for gain and readnoise ref files
        gain = np.ones(shape=(nrows, ncols), dtype=np.float64) * gain
        read_noise = np.full((nrows, ncols), readnoise, dtype=np.float64)

        # Create data array for science RampModel
        if subarray:
            data = np.zeros(shape=(nints, ngroups, 20, 20), dtype=np.float64)
        else:
            data = np.zeros(shape=(nints, ngroups, nrows, ncols), dtype=np.float64)

        # Create the science RampModel and populate meta data
        model1 = RampModel(data=data)
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

        # Create gain datamodel and populate meta
        gain = GainModel(data=gain)
        gain.meta.instrument.name = 'MIRI'
        gain.meta.subarray.xstart = 1
        gain.meta.subarray.ystart = 1
        gain.meta.subarray.xsize = ncols
        gain.meta.subarray.ysize = nrows

        # Create readnoise datamodel and populate meta
        rnModel = ReadnoiseModel(data=read_noise)
        rnModel.meta.instrument.name = 'MIRI'
        rnModel.meta.subarray.xstart = 1
        rnModel.meta.subarray.ystart = 1
        rnModel.meta.subarray.xsize = ncols
        rnModel.meta.subarray.ysize = nrows

        return model1, rnModel, gain

    return _setup
