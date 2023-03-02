import numpy as np
from numpy.testing import assert_array_equal
import pytest

from stdatamodels.jwst.datamodels import GainModel, ReadnoiseModel, RampModel, dqflags

from jwst.jump.jump import run_detect_jumps
import multiprocessing

import time
import random

JUMP_DET = dqflags.group["JUMP_DET"]
DO_NOT_USE = dqflags.group["DO_NOT_USE"]
GOOD = dqflags.group["GOOD"]
SATURATED = dqflags.group["SATURATED"]
NO_GAIN_VALUE = dqflags.pixel["NO_GAIN_VALUE"]


def test_exec_time_0_crs(setup_inputs):
    """"
    Set up with dimension similar to simulated MIRI datasets, Dataset has no
    cosmic rays. Test only the execution time of run_detect_jumps() for
    comparison with nominal time; hopefully indicative of faults with newly
    added code.
    """
    model, rnoise, gain = setup_inputs(ngroups=10, nrows=1024, ncols=1032,
                                       nints=2, readnoise=6.5, gain=5.5,
                                       grouptime=2.775, deltatime=2.775)

    tstart = time.time()
    # using dummy variable in next to prevent "F841-variable is assigned to but never used"
    _ = run_detect_jumps(model, gain, rnoise, 4.0, 5.0, 6.0, 'none', 200, 4, True)
    tstop = time.time()

    t_elapsed = tstop - tstart
    MAX_TIME = 10  # takes 1.6 sec on my Mac

    assert t_elapsed < MAX_TIME


def test_exec_time_many_crs(setup_inputs):
    """"
    Set up with dimension similar to simulated MIRI datasets, Dataset has
    many cosmic rays; approximately one CR per 4 groups. Test only the execution
    time of run_detect_jumps() for comparison with nominal time; hopefully
    indicative of faults with newly added code.
    """
    nrows = 350
    ncols = 400

    model, rnoise, gain = setup_inputs(ngroups=10, nrows=nrows, ncols=ncols,
                                       nints=2, readnoise=6.5, gain=5.5,
                                       grouptime=2.775, deltatime=2.775)

    crs_frac = 0.25  # fraction of groups having a CR
    model = add_crs(model, crs_frac)  # add desired fraction of CRs

    tstart = time.time()
    # using dummy variable in next to prevent "F841-variable is assigned to but never used"
    _ = run_detect_jumps(model, gain, rnoise, 4.0, 5.0, 6.0, 'none', 200, 4, True)
    tstop = time.time()

    t_elapsed = tstop - tstart
    MAX_TIME = 600  # takes ~100 sec on my Mac

    assert t_elapsed < MAX_TIME


def test_nocrs_noflux(setup_inputs):
    """"
    All pixel values are zero. So slope should be zero
    """
    model, rnoise, gain = setup_inputs(ngroups=5)
    out_model = run_detect_jumps(model, gain, rnoise, 4.0, 5.0, 6.0, 'none', 200, 4, True)
    assert np.max(out_model.groupdq) == GOOD


def test_nocrs_noflux_badgain_pixel(setup_inputs):
    """"
    all pixel values are zero. So slope should be zero, pixel with bad gain should
    have pixel dq set to 'NO_GAIN_VALUE' and 'DO_NOT_USE'
    """
    model, rnoise, gain = setup_inputs(ngroups=5, nrows=20, ncols=20)
    gain.data[7, 7] = -10  # bad gain
    gain.data[17, 17] = np.nan  # bad gain
    out_model = run_detect_jumps(model, gain, rnoise, 4.0, 5.0, 6.0, 'none', 200, 4, True)

    # 2 bits are set for each pixel, so use bitwise_and to check is set
    assert np.bitwise_and(out_model.pixeldq[7, 7], NO_GAIN_VALUE) == NO_GAIN_VALUE
    assert np.bitwise_and(out_model.pixeldq[7, 7], DO_NOT_USE) == DO_NOT_USE
    assert np.bitwise_and(out_model.pixeldq[17, 17], NO_GAIN_VALUE) == NO_GAIN_VALUE
    assert np.bitwise_and(out_model.pixeldq[17, 17], DO_NOT_USE) == DO_NOT_USE


def test_nocrs_noflux_subarray(setup_inputs):
    """"
    All pixel values are zero. This shows that the subarray reference files get
    extracted from the full frame versions.
    """
    model, rnoise, gain = setup_inputs(ngroups=5, subarray=True)
    out_model = run_detect_jumps(model, gain, rnoise, 4.0, 5.0, 6.0, 'none', 200, 4, True)
    assert np.max(out_model.groupdq) == GOOD


def test_onecr_10_groups_neighbors_flagged(setup_inputs):
    """"
    A single CR in a 10 group exposure
    """
    grouptime = 3.0
    ingain = 200
    inreadnoise = 7.0
    ngroups = 10
    model, rnoise, gain = setup_inputs(ngroups=ngroups, gain=ingain, nrows=10, ncols=10,
                                       readnoise=inreadnoise, deltatime=grouptime)

    # two segments perfect fit, second segment has twice the slope
    model.data[0, 0, 5, 5] = 15.0
    model.data[0, 1, 5, 5] = 20.0
    model.data[0, 2, 5, 5] = 25.0
    model.data[0, 3, 5, 5] = 30.0
    model.data[0, 4, 5, 5] = 35.0
    model.data[0, 5, 5, 5] = 140.0
    model.data[0, 6, 5, 5] = 150.0
    model.data[0, 7, 5, 5] = 160.0
    model.data[0, 8, 5, 5] = 170.0
    model.data[0, 9, 5, 5] = 180.0
    out_model = run_detect_jumps(model, gain, rnoise, 4.0, 5.0, 6.0, 'none', 200, 4, True)

    assert np.max(out_model.groupdq[0, 5, 5, 5]) == JUMP_DET
    assert out_model.groupdq[0, 5, 5, 6] == JUMP_DET
    assert out_model.groupdq[0, 5, 5, 4] == JUMP_DET
    assert out_model.groupdq[0, 5, 6, 5] == JUMP_DET
    assert out_model.groupdq[0, 5, 4, 5] == JUMP_DET


def test_nocr_100_groups_nframes1(setup_inputs):
    """"
    NO CR in a 100 group exposure to make sure that frames_per_group is passed correctly to
    twopoint_difference. This test recreates the problem found in issue #4571.
    """
    grouptime = 3.0
    ingain = 1
    inreadnoise = 7.0
    ngroups = 100
    model, rnoise, gain = setup_inputs(ngroups=ngroups, gain=ingain, nframes=1,
                                       nrows=10, ncols=10, readnoise=inreadnoise,
                                       deltatime=grouptime)

    # two segments perfect fit, second segment has twice the slope
    model.data[0, 0, 5, 5] = 14.0
    model.data[0, 1, 5, 5] = 20.0
    model.data[0, 2, 5, 5] = 27.0
    model.data[0, 3, 5, 5] = 30.0
    model.data[0, 4, 5, 5] = 38.0
    model.data[0, 5, 5, 5] = 40.0
    model.data[0, 6, 5, 5] = 50.0
    model.data[0, 7, 5, 5] = 52.0
    model.data[0, 8, 5, 5] = 63.0
    model.data[0, 9, 5, 5] = 68.0
    for i in range(10, 100):
        model.data[0, i, 5, 5] = i * 5
    out_model = run_detect_jumps(model, gain, rnoise, 4.0, 5.0, 6.0, 'none', 200, 4, True)
    assert np.max(out_model.groupdq) == GOOD


def test_twoints_onecr_each_10_groups_neighbors_flagged(setup_inputs):
    """"
    Two integrations with CRs in different locations. This makes sure we are correctly
    dealing with integrations.
    """
    grouptime = 3.0
    ingain = 200
    inreadnoise = 7.0
    ngroups = 10
    model, rnoise, gain = setup_inputs(ngroups=ngroups, nints=2, gain=ingain, nrows=20, ncols=20,
                                       readnoise=inreadnoise, deltatime=grouptime)

    # two segments perfect fit, second segment has twice the slope
    model.data[0, 0, 5, 5] = 15.0
    model.data[0, 1, 5, 5] = 20.0
    model.data[0, 2, 5, 5] = 25.0
    model.data[0, 3, 5, 5] = 30.0
    model.data[0, 4, 5, 5] = 35.0
    model.data[0, 5, 5, 5] = 140.0
    model.data[0, 6, 5, 5] = 150.0
    model.data[0, 7, 5, 5] = 160.0
    model.data[0, 8, 5, 5] = 170.0
    model.data[0, 9, 5, 5] = 180.0
    model.data[1, 0, 15, 5] = 15.0
    model.data[1, 1, 15, 5] = 20.0
    model.data[1, 2, 15, 5] = 25.0
    model.data[1, 3, 15, 5] = 30.0
    model.data[1, 4, 15, 5] = 35.0
    model.data[1, 5, 15, 5] = 40.0
    model.data[1, 6, 15, 5] = 45.0
    model.data[1, 7, 15, 5] = 160.0
    model.data[1, 8, 15, 5] = 170.0
    model.data[1, 9, 15, 5] = 180.0
    out_model = run_detect_jumps(model, gain, rnoise, 4.0, 5.0, 6.0, 'none', 200, 4, True)

    assert np.max(out_model.groupdq[0, 5, 5, 5]) == JUMP_DET
    assert out_model.groupdq[0, 5, 5, 6] == JUMP_DET
    assert out_model.groupdq[0, 5, 5, 4] == JUMP_DET
    assert out_model.groupdq[0, 5, 6, 5] == JUMP_DET
    assert out_model.groupdq[0, 5, 4, 5] == JUMP_DET
    assert out_model.groupdq[1, 7, 15, 5] == JUMP_DET
    assert out_model.groupdq[1, 7, 15, 6] == JUMP_DET
    assert out_model.groupdq[1, 7, 15, 4] == JUMP_DET
    assert out_model.groupdq[1, 7, 16, 5] == JUMP_DET
    assert out_model.groupdq[1, 7, 14, 5] == JUMP_DET


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

    model, rnoise, gain = setup_inputs(ngroups=ngroups, nints=1, nrows=nrows,
                                       ncols=ncols, gain=ingain,
                                       readnoise=inreadnoise, deltatime=grouptime)

    # Setup the desired pixel values
    model.data[0, :, 1, 1] = [10019.966, 10057.298, 10078.248, 10096.01,
                              20241.627, 20248.752, 20268.047, 20284.895, 20298.705, 20314.25]
    model.data[0, :, 1, 2] = [10016.457, 10053.907, 10063.568, 10076.166,
                              11655.773, 11654.063, 11681.795, 11693.763, 11712.788, 11736.994]
    model.data[0, :, 1, 3] = [10013.259, 10050.348, 10070.398, 10097.658, 10766.534, 10787.84,
                              10802.418, 10818.872, 10832.695, 10861.175]
    model.data[0, :, 1, 4] = [10016.422, 10053.959, 10070.934, 10090.381, 10104.014, 10127.665,
                              10143.687, 10172.227, 10178.138, 10199.59]
    model.data[0, :, 2, 1] = [10021.067, 10042.973, 10059.062, 10069.323, 18732.406, 18749.602,
                              18771.908, 18794.695, 18803.223, 18819.523]
    model.data[0, :, 2, 2] = [10019.651, 10043.371, 10056.423, 10085.121, 40584.703, 40606.08,
                              40619.51, 40629.574, 40641.9, 40660.145]
    model.data[0, :, 2, 3] = [10021.223, 10042.112, 10052.958, 10067.142, 28188.316, 28202.922,
                              28225.557, 28243.79, 28253.883, 28273.586]
    model.data[0, :, 2, 4] = [10022.608, 10037.174, 10069.476, 10081.729, 11173.748, 11177.344,
                              11201.127, 11219.607, 11229.468, 11243.174]
    model.data[0, :, 2, 5] = [10011.095, 10047.422, 10061.066, 10079.375, 10106.405, 10116.071,
                              10129.348, 10136.305, 10161.373, 10181.479]
    model.data[0, :, 3, 1] = [10011.877, 10052.809, 10075.108, 10085.111, 10397.106, 10409.291,
                              10430.475, 10445.3, 10462.004, 10484.906]
    model.data[0, :, 3, 2] = [10012.124, 10059.202, 10078.984, 10092.74, 11939.488,
                              11958.45, 11977.5625, 11991.776, 12025.897, 12027.326]
    model.data[0, :, 3, 3] = [10013.282, 10046.887, 10062.308, 10085.447, 28308.426, 28318.957,
                              28335.55, 28353.832, 28371.746, 28388.848]
    model.data[0, :, 3, 4] = [10016.784, 10048.249, 10060.097, 10074.606, 21506.082, 21522.027,
                              21542.309, 21558.34, 21576.365, 21595.58]
    model.data[0, :, 3, 5] = [10014.916, 10052.995, 10063.7705, 10092.866, 10538.075,
                              10558.318, 10570.754, 10597.343, 10608.488, 10628.104]
    model.data[0, :, 4, 1] = [10017.438, 10038.94, 10057.657, 10069.987, 10090.22, 10114.296,
                              10133.543, 10148.657, 10158.109, 10172.842]
    model.data[0, :, 4, 2] = [10011.129, 10037.982, 10054.445, 10079.703, 10097.964, 10110.593,
                              10135.701, 10149.448, 10171.771, 10185.874]
    model.data[0, :, 4, 3] = [10021.109, 10043.658, 10063.909, 10072.364, 10766.232, 10774.402,
                              10790.677, 10809.337, 10833.65, 10849.55]
    model.data[0, :, 4, 4] = [10023.877, 10035.997, 10052.321, 10077.937, 10529.645, 10541.947,
                              10571.127, 10577.249, 10599.716, 10609.544]

    # Flag first and last frame as DO_NOT_USE
    model.groupdq[0, 0, :, :] = 1
    model.groupdq[0, -1, :, :] = 1

    # run jump detection
    out_model = run_detect_jumps(model, gain, rnoise, rejection_thresh=200.0, three_grp_thresh=200,
                                 four_grp_thresh=200,
                                 max_cores='none', max_jump_to_flag_neighbors=200,
                                 min_jump_to_flag_neighbors=10, flag_4_neighbors=True)

    # Verify that the correct groups have been flagged. The entries for pixels
    # 2,2 and 3,3 are the ones that had previously been flagged in group 2 instead
    # of group 5.
    assert_array_equal(out_model.groupdq[0, :, 1, 1], [1, 0, 0, 0, 4, 0, 0, 0, 0, 1])
    assert_array_equal(out_model.groupdq[0, :, 1, 2], [1, 0, 0, 0, 4, 0, 0, 0, 0, 1])
    assert_array_equal(out_model.groupdq[0, :, 1, 3], [1, 0, 0, 0, 0, 0, 0, 0, 0, 1])
    assert_array_equal(out_model.groupdq[0, :, 1, 4], [1, 0, 0, 0, 0, 0, 0, 0, 0, 1])
    assert_array_equal(out_model.groupdq[0, :, 2, 1], [1, 0, 0, 0, 4, 0, 0, 0, 0, 1])
    assert_array_equal(out_model.groupdq[0, :, 2, 2], [1, 0, 0, 0, 4, 0, 0, 0, 0, 1])
    assert_array_equal(out_model.groupdq[0, :, 2, 3], [1, 0, 0, 0, 4, 0, 0, 0, 0, 1])
    assert_array_equal(out_model.groupdq[0, :, 2, 4], [1, 0, 0, 0, 0, 0, 0, 0, 0, 1])
    assert_array_equal(out_model.groupdq[0, :, 3, 1], [1, 0, 0, 0, 0, 0, 0, 0, 0, 1])
    assert_array_equal(out_model.groupdq[0, :, 3, 2], [1, 0, 0, 0, 4, 0, 0, 0, 0, 1])
    assert_array_equal(out_model.groupdq[0, :, 3, 3], [1, 0, 0, 0, 4, 0, 0, 0, 0, 1])
    assert_array_equal(out_model.groupdq[0, :, 3, 4], [1, 0, 0, 0, 4, 0, 0, 0, 0, 1])
    assert_array_equal(out_model.groupdq[0, :, 4, 1], [1, 0, 0, 0, 0, 0, 0, 0, 0, 1])
    assert_array_equal(out_model.groupdq[0, :, 4, 2], [1, 0, 0, 0, 0, 0, 0, 0, 0, 1])
    assert_array_equal(out_model.groupdq[0, :, 4, 3], [1, 0, 0, 0, 0, 0, 0, 0, 0, 1])
    assert_array_equal(out_model.groupdq[0, :, 4, 4], [1, 0, 0, 0, 0, 0, 0, 0, 0, 1])


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

    model, rnoise, gain = setup_inputs(ngroups=ngroups, nints=1, nrows=nrows,
                                       ncols=ncols, gain=ingain,
                                       readnoise=inreadnoise, deltatime=grouptime)

    # Setup the needed input pixel and DQ values
    model.data[0, :, 1, 1] = [639854.75, 4872.451, -17861.791, 14022.15, 22320.176,
                              1116.3828, 1936.9746]
    model.groupdq[0, :, 1, 1] = [0, 0, 0, 0, 0, 2, 2]
    model.data[0, :, 2, 2] = [8.25666812e+05, -1.10471914e+05, 1.95755371e+02, 1.83118457e+03,
                              1.72250879e+03, 1.81733496e+03, 1.65188281e+03]
    # 2 non-sat groups means only 1 non-sat diff, so no jumps should be flagged
    model.groupdq[0, :, 2, 2] = [0, 0, 2, 2, 2, 2, 2]
    model.data[0, :, 3, 3] = [1228767., 46392.234, -3245.6553, 7762.413,
                              37190.76, 266611.62, 5072.4434]
    model.groupdq[0, :, 3, 3] = [0, 0, 0, 0, 0, 0, 2]

    # run jump detection
    out_model = run_detect_jumps(model, gain, rnoise, rejection_thresh=200.0, three_grp_thresh=200,
                                 four_grp_thresh=200,
                                 max_cores='none', max_jump_to_flag_neighbors=200,
                                 min_jump_to_flag_neighbors=10, flag_4_neighbors=True)

    # Check the results. There should not be any pixels with DQ values of 6, which
    # is saturated (2) plus jump (4). All the DQ's should be either just 2 or just 4.
    assert_array_equal(out_model.groupdq[0, :, 1, 1], [0, 4, 4, 4, 0, 2, 2])
    # assert that no groups are flagged when theres only 1 non-sat. grp
    assert_array_equal(out_model.groupdq[0, :, 2, 2], [0, 0, 2, 2, 2, 2, 2])
    assert_array_equal(out_model.groupdq[0, :, 3, 3], [0, 4, 4, 0, 0, 4, 2])


def add_crs(model, crs_frac):
    """"
    Randomly add a cosmic ray of magnitude CR_MAG to a fraction (crs_frac)
    of the groups in the model's SCI array
    """

    num_ints = model.data.shape[0]
    num_groups = model.data.shape[1]
    num_rows = model.data.shape[2]
    num_cols = model.data.shape[3]

    tot_cr = 0  # counter
    CR_MAG = 1000.  # consider making a variable ?

    random.seed(0)  # to generate same CRs

    # Add to the model's data, in all but the 0th group
    for ii_int in range(num_ints):  # loop over integrations
        for ii_col in range(num_cols):
            for ii_row in range(num_rows):
                for ii_group in range(1, num_groups):
                    cr_rand = random.random()
                    if cr_rand < crs_frac:
                        tot_cr += 1
                        model.data[ii_int, ii_group:, ii_row, ii_col] += CR_MAG

    return model


def test_flagging_of_CRs_across_slice_boundaries(setup_inputs):
    """"
    A multiprocessing test that has two CRs on the boundary between two slices.
    This makes sure that we are correctly flagging neighbors in different  slices.
    """
    grouptime = 3.0
    ingain = 200
    inreadnoise = 7.0
    ngroups = 10

    model, rnoise, gain = setup_inputs(ngroups=ngroups, nints=2, gain=ingain,
                                       readnoise=inreadnoise, deltatime=grouptime)
    nrows = model.data.shape[3]
    num_cores = multiprocessing.cpu_count()
    max_cores = 'half'
    numslices = num_cores // 2
    if numslices > 1:
        yincrement = int(nrows / numslices)
        # two segments perfect fit, second segment has twice the slope
        # add a CR on the last row of the first slice
        model.data[0, 0, yincrement - 1, 5] = 15.0
        model.data[0, 1, yincrement - 1, 5] = 20.0
        model.data[0, 2, yincrement - 1, 5] = 25.0
        model.data[0, 3, yincrement - 1, 5] = 30.0
        model.data[0, 4, yincrement - 1, 5] = 35.0
        model.data[0, 5, yincrement - 1, 5] = 140.0
        model.data[0, 6, yincrement - 1, 5] = 150.0
        model.data[0, 7, yincrement - 1, 5] = 160.0
        model.data[0, 8, yincrement - 1, 5] = 170.0
        model.data[0, 9, yincrement - 1, 5] = 180.0
        # add a CR on the first row of the second slice
        model.data[1, 0, yincrement, 25] = 15.0
        model.data[1, 1, yincrement, 25] = 20.0
        model.data[1, 2, yincrement, 25] = 25.0
        model.data[1, 3, yincrement, 25] = 30.0
        model.data[1, 4, yincrement, 25] = 35.0
        model.data[1, 5, yincrement, 25] = 40.0
        model.data[1, 6, yincrement, 25] = 50.0
        model.data[1, 7, yincrement, 25] = 160.0
        model.data[1, 8, yincrement, 25] = 170.0
        model.data[1, 9, yincrement, 25] = 180.0
        out_model = run_detect_jumps(model, gain, rnoise, 4.0, 5.0, 6.0, max_cores, 200, 4, True)

        # check that the neighbors of the CR on the last row were flagged
        assert out_model.groupdq[0, 5, yincrement - 1, 5] == JUMP_DET
        assert out_model.groupdq[0, 5, yincrement - 1, 6] == JUMP_DET
        assert out_model.groupdq[0, 5, yincrement - 1, 4] == JUMP_DET
        assert out_model.groupdq[0, 5, yincrement, 5] == JUMP_DET
        assert out_model.groupdq[0, 5, yincrement - 2, 5] == JUMP_DET
        # check that the neighbors of the CR on the first row were flagged
        assert out_model.groupdq[1, 7, yincrement, 25] == JUMP_DET
        assert out_model.groupdq[1, 7, yincrement, 26] == JUMP_DET
        assert out_model.groupdq[1, 7, yincrement, 24] == JUMP_DET
        assert out_model.groupdq[1, 7, yincrement + 1, 25] == JUMP_DET
        assert out_model.groupdq[1, 7, yincrement - 1, 25] == JUMP_DET


def test_twoints_onecr_10_groups_neighbors_flagged_multi(setup_inputs):
    """"
    A multiprocessing test that has two CRs on the boundary between two slices
    in different integrations. This makes sure that we are correctly flagging
    neighbors in different slices and that we are parsing the integrations correctly.
    """
    grouptime = 3.0
    ingain = 200
    inreadnoise = 7.0
    ngroups = 10
    model, rnoise, gain = setup_inputs(ngroups=ngroups, nints=2, gain=ingain,
                                       nrows=40, ncols=10,
                                       readnoise=inreadnoise, deltatime=grouptime)

    # two segments perfect fit, second segment has twice the slope
    model.data[0, 0, 5, 5] = 15.0
    model.data[0, 1, 5, 5] = 20.0
    model.data[0, 2, 5, 5] = 25.0
    model.data[0, 3, 5, 5] = 30.0
    model.data[0, 4, 5, 5] = 35.0
    model.data[0, 5, 5, 5] = 140.0
    model.data[0, 6, 5, 5] = 150.0
    model.data[0, 7, 5, 5] = 160.0
    model.data[0, 8, 5, 5] = 170.0
    model.data[0, 9, 5, 5] = 180.0
    model.data[1, 0, 15, 5] = 15.0
    model.data[1, 1, 15, 5] = 20.0
    model.data[1, 2, 15, 5] = 25.0
    model.data[1, 3, 15, 5] = 30.0
    model.data[1, 4, 15, 5] = 35.0
    model.data[1, 5, 15, 5] = 40.0
    model.data[1, 6, 15, 5] = 45.0
    model.data[1, 7, 15, 5] = 160.0
    model.data[1, 8, 15, 5] = 170.0
    model.data[1, 9, 15, 5] = 180.0
    out_model = run_detect_jumps(model, gain, rnoise, 4.0, 5.0, 6.0, 'half', 200, 4, True)

    assert np.max(out_model.groupdq[0, 5, 5, 5]) == JUMP_DET
    assert out_model.groupdq[0, 5, 5, 6] == JUMP_DET
    assert out_model.groupdq[0, 5, 5, 4] == JUMP_DET
    assert out_model.groupdq[0, 5, 6, 5] == JUMP_DET
    assert out_model.groupdq[0, 5, 4, 5] == JUMP_DET
    assert out_model.groupdq[1, 7, 15, 5] == JUMP_DET
    assert out_model.groupdq[1, 7, 15, 6] == JUMP_DET
    assert out_model.groupdq[1, 7, 15, 4] == JUMP_DET
    assert out_model.groupdq[1, 7, 16, 5] == JUMP_DET
    assert out_model.groupdq[1, 7, 14, 5] == JUMP_DET


@pytest.mark.skip(reason="Test is only used to test performance issue. No need to run every time.")
def test_every_pixel_CR_neighbors_flagged(setup_inputs):
    """"
    A multiprocessing test that has a jump in every pixel. This is used
    to test the performance gain from multiprocessing.
    """
    grouptime = 3.0
    ingain = 200
    inreadnoise = 7.0
    ngroups = 10
    model, rnoise, gain = setup_inputs(ngroups=ngroups, gain=ingain, nrows=100, ncols=100,
                                       readnoise=inreadnoise, deltatime=grouptime)

    # two segments perfect fit, second segment has twice the slope
    model.data[0, 0, :, :] = 15.0
    model.data[0, 1, :, :] = 20.0
    model.data[0, 2, :, :] = 25.0
    model.data[0, 3, :, :] = 30.0
    model.data[0, 4, :, :] = 35.0
    model.data[0, 5, :, :] = 140.0
    model.data[0, 6, :, :] = 150.0
    model.data[0, 7, :, :] = 160.0
    model.data[0, 8, :, :] = 170.0
    model.data[0, 9, :, :] = 180.0
    out_model = run_detect_jumps(model, gain, rnoise, 4.0, 5.0, 6.0, 'half', 200, 4, True)

    assert np.max(out_model.groupdq[0, 5, 5, 5]) == JUMP_DET
    assert out_model.groupdq[0, 5, 5, 6] == JUMP_DET
    assert out_model.groupdq[0, 5, 5, 4] == JUMP_DET
    assert out_model.groupdq[0, 5, 6, 5] == JUMP_DET
    assert out_model.groupdq[0, 5, 4, 5] == JUMP_DET


def test_crs_on_edge_with_neighbor_flagging(setup_inputs):
    """"
    A test to make sure that the neighbors of CRs on the edges of the
    array are flagged correctly.
    """
    grouptime = 3.0
    ingain = 200
    inreadnoise = 7.0
    ngroups = 10
    model, rnoise, gain = setup_inputs(ngroups=ngroups, nrows=20, ncols=20,
                                       gain=ingain, readnoise=inreadnoise,
                                       deltatime=grouptime)

    # two segments perfect fit, second segment has twice the slope
    # CR on 1st row
    model.data[0, 0, 0, 15] = 15.0
    model.data[0, 1, 0, 15] = 20.0
    model.data[0, 2, 0, 15] = 25.0
    model.data[0, 3, 0, 15] = 30.0
    model.data[0, 4, 0, 15] = 35.0
    model.data[0, 5, 0, 15] = 140.0
    model.data[0, 6, 0, 15] = 150.0
    model.data[0, 7, 0, 15] = 160.0
    model.data[0, 8, 0, 15] = 170.0
    model.data[0, 9, 0, 15] = 180.0
    # CR on last row
    model.data[0, 0, -1, 5] = 15.0
    model.data[0, 1, -1, 5] = 20.0
    model.data[0, 2, -1, 5] = 25.0
    model.data[0, 3, -1, 5] = 30.0
    model.data[0, 4, -1, 5] = 35.0
    model.data[0, 5, -1, 5] = 140.0
    model.data[0, 6, -1, 5] = 150.0
    model.data[0, 7, -1, 5] = 160.0
    model.data[0, 8, -1, 5] = 170.0
    model.data[0, 9, -1, 5] = 180.0
    # CR on 1st column
    model.data[0, 0, 5, 0] = 15.0
    model.data[0, 1, 5, 0] = 20.0
    model.data[0, 2, 5, 0] = 25.0
    model.data[0, 3, 5, 0] = 30.0
    model.data[0, 4, 5, 0] = 35.0
    model.data[0, 5, 5, 0] = 140.0
    model.data[0, 6, 5, 0] = 150.0
    model.data[0, 7, 5, 0] = 160.0
    model.data[0, 8, 5, 0] = 170.0
    model.data[0, 9, 5, 0] = 180.0
    # CR on last column
    model.data[0, 0, 15, -1] = 15.0
    model.data[0, 1, 15, -1] = 20.0
    model.data[0, 2, 15, -1] = 25.0
    model.data[0, 3, 15, -1] = 30.0
    model.data[0, 4, 15, -1] = 35.0
    model.data[0, 5, 15, -1] = 140.0
    model.data[0, 6, 15, -1] = 150.0
    model.data[0, 7, 15, -1] = 160.0
    model.data[0, 8, 15, -1] = 170.0
    model.data[0, 9, 15, -1] = 180.0

    out_model = run_detect_jumps(model, gain, rnoise, 4.0, 5.0, 6.0, 'none', 200, 10, True)

    # flag CR and three neighbors of first row CR
    assert out_model.groupdq[0, 5, 0, 15] == JUMP_DET
    assert out_model.groupdq[0, 5, 1, 15] == JUMP_DET
    assert out_model.groupdq[0, 5, 0, 14] == JUMP_DET
    assert out_model.groupdq[0, 5, 0, 16] == JUMP_DET
    assert out_model.groupdq[0, 5, -1, 15] == 0  # The one not to flag
    # flag CR and three neighbors of last row CR
    assert out_model.groupdq[0, 5, -1, 5] == JUMP_DET
    assert out_model.groupdq[0, 5, -2, 5] == JUMP_DET
    assert out_model.groupdq[0, 5, -1, 4] == JUMP_DET
    assert out_model.groupdq[0, 5, -1, 6] == JUMP_DET
    # flag CR and three neighbors of first column CR
    assert out_model.groupdq[0, 5, 5, 0] == JUMP_DET
    assert out_model.groupdq[0, 5, 6, 0] == JUMP_DET
    assert out_model.groupdq[0, 5, 4, 0] == JUMP_DET
    assert out_model.groupdq[0, 5, 5, 1] == JUMP_DET
    assert out_model.groupdq[0, 5, 5, -1] == 0  # The one not to flag
    # flag CR and three neighbors of last column CR
    assert out_model.groupdq[0, 5, 15, -1] == JUMP_DET
    assert out_model.groupdq[0, 5, 15, -2] == JUMP_DET
    assert out_model.groupdq[0, 5, 16, -1] == JUMP_DET
    assert out_model.groupdq[0, 5, 14, -1] == JUMP_DET


def test_onecr_10_groups(setup_inputs):
    """"
    A test to make sure that neighbors are not flagged when they are not requested to be flagged.
    """
    grouptime = 3.0
    ingain = 200
    inreadnoise = 7.0
    ngroups = 10
    model, rnoise, gain = setup_inputs(ngroups=ngroups, gain=ingain, nrows=20, ncols=20,
                                       readnoise=inreadnoise, deltatime=grouptime)
    # two segments perfect fit, second segment has twice the slope
    model.data[0, 0, 5, 5] = 15.0
    model.data[0, 1, 5, 5] = 20.0
    model.data[0, 2, 5, 5] = 25.0
    model.data[0, 3, 5, 5] = 30.0
    model.data[0, 4, 5, 5] = 35.0
    model.data[0, 5, 5, 5] = 140.0
    model.data[0, 6, 5, 5] = 150.0
    model.data[0, 7, 5, 5] = 160.0
    model.data[0, 8, 5, 5] = 170.0
    model.data[0, 9, 5, 5] = 180.0
    out_model = run_detect_jumps(model, gain, rnoise, 4.0, 5.0, 6.0, 'none', 200, 10, False)

    assert out_model.groupdq[0, 5, 5, 5] == JUMP_DET
    assert out_model.groupdq[0, 5, 4, 5] == GOOD
    assert out_model.groupdq[0, 5, 6, 5] == GOOD
    assert out_model.groupdq[0, 5, 5, 6] == GOOD
    assert out_model.groupdq[0, 5, 5, 4] == GOOD


def test_onecr_10_groups_fullarray(setup_inputs):
    """"
    A test that has a cosmic ray in the 5th group for all pixels except column 10. In column
    10 the jump is in the 7th group.
    """
    grouptime = 3.0
    ingain = 5
    inreadnoise = 7.0
    ngroups = 10
    model, rnoise, gain = setup_inputs(ngroups=ngroups, gain=ingain, nrows=20, ncols=20,
                                       readnoise=inreadnoise, deltatime=grouptime)

    model.data[0, 0, 5, :] = 15.0
    model.data[0, 1, 5, :] = 20.0
    model.data[0, 2, 5, :] = 25.0
    model.data[0, 3, 5, :] = 30.0
    model.data[0, 4, 5, :] = 35.0
    model.data[0, 5, 5, :] = 140.0
    model.data[0, 6, 5, :] = 150.0
    model.data[0, 7, 5, :] = 160.0
    model.data[0, 8, 5, :] = 170.0
    model.data[0, 9, 5, :] = 180.0
    # move the CR to group 7 for row 10 and make difference be 300
    model.data[0, 3, 5, 10] = 100
    model.data[0, 4, 5, 10] = 130
    model.data[0, 5, 5, 10] = 160
    model.data[0, 6, 5, 10] = 190
    model.data[0, 7, 5, 10] = 400
    model.data[0, 8, 5, 10] = 410
    model.data[0, 9, 5, 10] = 420
    out_model = run_detect_jumps(model, gain, rnoise, 4.0, 5.0, 6.0, 'none', 200, 10, False)

    # The jump is in group 5 for columns 0-9
    assert_array_equal(out_model.groupdq[0, 5, 5, 0:10], JUMP_DET)
    # The jump is in group 7 for column 10
    assert out_model.groupdq[0, 7, 5, 10] == JUMP_DET
    # The jump is in group 5 for columns 11+
    assert_array_equal(out_model.groupdq[0, 5, 5, 11:], JUMP_DET)


def test_onecr_50_groups(setup_inputs):
    """"
    A test with a fifty group integration. There are two jumps in pixel 5,5. One in group 5 and
    one in group 30.
    """
    grouptime = 3.0
    ingain = 5
    inreadnoise = 7.0
    ngroups = 50
    model, rnoise, gain = setup_inputs(ngroups=ngroups, gain=ingain, nrows=10, ncols=10,
                                       readnoise=inreadnoise, deltatime=grouptime)

    model.data[0, 0, 5, 5] = 15.0
    model.data[0, 1, 5, 5] = 20.0
    model.data[0, 2, 5, 5] = 25.0
    model.data[0, 3, 5, 5] = 30.0
    model.data[0, 4, 5, 5] = 35.0
    model.data[0, 5, 5, 5] = 140.0
    model.data[0, 6, 5, 5] = 150.0
    model.data[0, 7, 5, 5] = 160.0
    model.data[0, 8, 5, 5] = 170.0
    model.data[0, 9, 5, 5] = 180.0
    model.data[0, 10:30, 5, 5] = np.arange(190, 290, 5)
    model.data[0, 30:50, 5, 5] = np.arange(500, 600, 5)
    out_model = run_detect_jumps(model, gain, rnoise, 4.0, 5.0, 6.0, 'none', 200, 10, False)

    # CR in group 5
    assert out_model.groupdq[0, 5, 5, 5] == JUMP_DET
    # CR in group 30
    assert out_model.groupdq[0, 30, 5, 5] == JUMP_DET
    # groups in between are not flagged
    assert_array_equal(out_model.groupdq[0, 6:30, 5, 5], GOOD)


def test_onecr_50_groups_afterjump(setup_inputs):
    """"
    A test with a fifty group integration. There are two jumps in pixel 5,5. One in group 5 and
    one in group 30.  Test includes after jump flagging.
    """
    grouptime = 3.0
    ingain = 5
    inreadnoise = 7.0
    ngroups = 50
    model, rnoise, gain = setup_inputs(ngroups=ngroups, gain=ingain, nrows=10, ncols=10,
                                       readnoise=inreadnoise, deltatime=grouptime)

    model.data[0, 0, 5, 5] = 15.0
    model.data[0, 1, 5, 5] = 20.0
    model.data[0, 2, 5, 5] = 25.0
    model.data[0, 3, 5, 5] = 30.0
    model.data[0, 4, 5, 5] = 35.0
    model.data[0, 5, 5, 5] = 140.0
    model.data[0, 6, 5, 5] = 150.0
    model.data[0, 7, 5, 5] = 160.0
    model.data[0, 8, 5, 5] = 170.0
    model.data[0, 9, 5, 5] = 180.0
    model.data[0, 10:30, 5, 5] = np.arange(190, 290, 5)
    model.data[0, 30:50, 5, 5] = np.arange(500, 600, 5)
    out_model = run_detect_jumps(model, gain, rnoise,
                                 4.0, 5.0, 6.0, 'none', 200, 10, False,
                                 after_jump_flag_dn1=20.,
                                 after_jump_flag_time1=grouptime * 2,
                                 after_jump_flag_dn2=150.,
                                 after_jump_flag_time2=grouptime * 3)

    # CR in group 5 + 2 groups
    for k in range(5, 8):
        assert out_model.groupdq[0, k, 5, 5] == JUMP_DET, f"first cr group {k}"
    # CR in group 30 + 3 groups
    for k in range(30, 34):
        assert out_model.groupdq[0, k, 5, 5] == JUMP_DET, f"second cr group {k}"
    # groups in between are not flagged
    assert_array_equal(out_model.groupdq[0, 8:30, 5, 5], GOOD, err_msg="between crs")


def test_single_CR_neighbor_flag(setup_inputs):
    """"
    A single CR in a 10 group exposure. Tests that:
    - if neighbor-flagging is set, the 4 neighboring pixels *ARE* flagged, and
    - if neighbor-flagging is *NOT* set, the 4 neighboring pixels are *NOT* flagged
    """
    grouptime = 3.0
    ingain = 5
    inreadnoise = 7.0
    ngroups = 10

    model, rnoise, gain = setup_inputs(ngroups=ngroups, nrows=5, ncols=6,
                                       gain=ingain, readnoise=inreadnoise,
                                       deltatime=grouptime)

    # two segments perfect fit, second segment has twice the slope
    model.data[0, 0, 3, 3] = 15.0
    model.data[0, 1, 3, 3] = 20.0
    model.data[0, 2, 3, 3] = 25.0
    model.data[0, 3, 3, 3] = 30.0
    model.data[0, 4, 3, 3] = 35.0
    model.data[0, 5, 3, 3] = 140.0
    model.data[0, 6, 3, 3] = 150.0
    model.data[0, 7, 3, 3] = 160.0
    model.data[0, 8, 3, 3] = 170.0
    model.data[0, 9, 3, 3] = 180.0

    # Flag neighbors
    out_model = run_detect_jumps(model, gain, rnoise, 4.0, 5.0, 6.0, 'none', 200, 4, True)

    assert np.max(out_model.groupdq[0, 5, 3, 3]) == JUMP_DET
    assert out_model.groupdq[0, 5, 3, 4] == JUMP_DET
    assert out_model.groupdq[0, 5, 3, 2] == JUMP_DET
    assert out_model.groupdq[0, 5, 2, 3] == JUMP_DET
    assert out_model.groupdq[0, 5, 4, 3] == JUMP_DET

    # Do not flag neighbors
    out_model = run_detect_jumps(model, gain, rnoise, 4.0, 5.0, 6.0, 'none', 200, 4, False)

    assert np.max(out_model.groupdq[0, 5, 3, 3]) == JUMP_DET
    assert out_model.groupdq[0, 5, 3, 4] == GOOD
    assert out_model.groupdq[0, 5, 3, 2] == GOOD
    assert out_model.groupdq[0, 5, 2, 3] == GOOD
    assert out_model.groupdq[0, 5, 4, 3] == GOOD


def test_proc(setup_inputs):
    """"
    A single CR in a 10 group exposure. Verify that the pixels flagged using
    multiprocessing are identical to the pixels flagged when no
    multiprocessing is done.
    """
    grouptime = 3.0
    ingain = 5
    inreadnoise = 7.0
    ngroups = 10

    model, rnoise, gain = setup_inputs(ngroups=ngroups, nrows=25, ncols=6,
                                       nints=2, gain=ingain,
                                       readnoise=inreadnoise, deltatime=grouptime)

    model.data[0, 0, 2, 3] = 15.0
    model.data[0, 1, 2, 3] = 21.0
    model.data[0, 2, 2, 3] = 25.0
    model.data[0, 3, 2, 3] = 30.2
    model.data[0, 4, 2, 3] = 35.0
    model.data[0, 5, 2, 3] = 140.0
    model.data[0, 6, 2, 3] = 151.0
    model.data[0, 7, 2, 3] = 160.0
    model.data[0, 8, 2, 3] = 170.0
    model.data[0, 9, 2, 3] = 180.0

    out_model_a = run_detect_jumps(model, gain, rnoise, 4.0, 5.0, 6.0, 'none', 200, 4, True)
    out_model_b = run_detect_jumps(model, gain, rnoise, 4.0, 5.0, 6.0, 'half', 200, 4, True)
    assert_array_equal(out_model_a.groupdq, out_model_b.groupdq)

    out_model_c = run_detect_jumps(model, gain, rnoise, 4.0, 5.0, 6.0, 'all', 200, 4, True)
    assert_array_equal(out_model_a.groupdq, out_model_c.groupdq)


def test_adjacent_CRs(setup_inputs):
    """
    Three CRs in a 10 group exposure; the CRs have overlapping neighboring
    pixels. This test makes sure that the correct pixels are flagged.
    """
    grouptime = 3.0
    ingain = 5
    inreadnoise = 7.0
    ngroups = 10
    model, rnoise, gain = setup_inputs(ngroups=ngroups, nrows=15, ncols=6,
                                       gain=ingain, readnoise=inreadnoise,
                                       deltatime=grouptime)

    # Populate arrays for 1st CR, centered at (x=2, y=3)
    x = 2
    y = 3
    model.data[0, 0, y, x] = 15.0
    model.data[0, 1, y, x] = 20.0
    model.data[0, 2, y, x] = 26.0
    model.data[0, 3, y, x] = 30.0
    model.data[0, 4, y, x] = 35.0
    model.data[0, 5, y, x] = 140.0
    model.data[0, 6, y, x] = 150.0
    model.data[0, 7, y, x] = 161.0
    model.data[0, 8, y, x] = 170.0
    model.data[0, 9, y, x] = 180.0

    # Populate arrays for 2nd CR, centered at (x=2, y=2)
    x = 2
    y = 2
    model.data[0, 0, y, x] = 20.0
    model.data[0, 1, y, x] = 30.0
    model.data[0, 2, y, x] = 41.0
    model.data[0, 3, y, x] = 51.0
    model.data[0, 4, y, x] = 62.0
    model.data[0, 5, y, x] = 170.0
    model.data[0, 6, y, x] = 200.0
    model.data[0, 7, y, x] = 231.0
    model.data[0, 8, y, x] = 260.0
    model.data[0, 9, y, x] = 290.0

    # Populate arrays for 3rd CR, centered at (x=3, y=2)
    x = 3
    y = 2
    model.data[0, 0, y, x] = 120.0
    model.data[0, 1, y, x] = 140.0
    model.data[0, 2, y, x] = 161.0
    model.data[0, 3, y, x] = 181.0
    model.data[0, 4, y, x] = 202.0
    model.data[0, 5, y, x] = 70.0
    model.data[0, 6, y, x] = 100.0
    model.data[0, 7, y, x] = 131.0
    model.data[0, 8, y, x] = 160.0
    model.data[0, 9, y, x] = 190.0

    out_model = run_detect_jumps(model, gain, rnoise, 4.0, 5.0, 6.0, 'half', 200, 4, True)

    # 1st CR (centered at x=2, y=3)
    assert out_model.groupdq[0, 5, 2, 2] == JUMP_DET
    assert out_model.groupdq[0, 5, 3, 1] == JUMP_DET
    assert out_model.groupdq[0, 5, 3, 2] == JUMP_DET
    assert out_model.groupdq[0, 5, 3, 3] == JUMP_DET
    assert out_model.groupdq[0, 5, 4, 2] == JUMP_DET

    # 2nd CR (centered at x=2, y=2)
    assert out_model.groupdq[0, 5, 1, 2] == JUMP_DET
    assert out_model.groupdq[0, 5, 2, 1] == JUMP_DET
    assert out_model.groupdq[0, 5, 2, 3] == JUMP_DET

    # 3rd CR (centered at x=3, y=2)
    assert out_model.groupdq[0, 5, 1, 3] == JUMP_DET
    assert out_model.groupdq[0, 5, 2, 4] == JUMP_DET


# Need test for multi-ints near zero with positive and negative slopes


def test_cr_neighbor_sat_flagging(setup_inputs):
    """"
    For pixels impacted by cosmic rays, neighboring pixels that are not
    already flagged as saturated will also be flagged as JUMP_DET. If a
    neigboring pixel has been flagged as saturated, it should not also be
    flagged as JUMP_DET.  This test checks if both of those types of
    neighboring pixels are correctly flagged, and tests more cases than
    covered by the tests test_nirspec_saturated_pix() and
    test_crs_on_edge_with_neighbor_flagging().
    """
    SAT_SCI = 1E4  # dummy saturation level

    grouptime = 3.0
    ingain = 200
    inreadnoise = 7.0
    ngroups = 8
    model, rnoise, gain = setup_inputs(ngroups=ngroups, gain=ingain, nrows=6,
                                       ncols=7, readnoise=inreadnoise,
                                       deltatime=grouptime)

    # Construct ramps with cosmic rays having some neighboring SAT pixels:
    # CR (A): group 2, in corner
    model.data[0, 0, 0, 0] = 15.0
    model.data[0, 1, 0, 0] = 20.0
    model.data[0, 2, 0, 0] = 120.0
    model.data[0, 3, 0, 0] = 125.0
    model.data[0, 4, 0, 0] = 130.0
    model.data[0, 5, 0, 0] = 135.0
    model.data[0, 6, 0, 0] = 140.0
    model.data[0, 7, 0, 0] = 145.0
    # ... and set 1 of the 2 neighbors to SAT, and adjust its ramp
    model.groupdq[0, 2:, 0, 1] = SATURATED
    model.data[0, 2:, 0, 1] = SAT_SCI

    # CR (B): group 2, not on edge
    model.data[0, 0, 1, 4] = 30.0
    model.data[0, 1, 1, 4] = 40.0
    model.data[0, 2, 1, 4] = 160.0
    model.data[0, 3, 1, 4] = 165.0
    model.data[0, 4, 1, 4] = 170.0
    model.data[0, 5, 1, 4] = 175.0
    model.data[0, 6, 1, 4] = 180.0
    model.data[0, 7, 1, 4] = 185.0
    # ... and set 2 of the 4 neighbors to SAT, and adjust their ramps
    model.groupdq[0, 2:, 0, 4] = SATURATED
    model.groupdq[0, 2:, 1, 3] = SATURATED
    model.data[0, 2:, 0, 4] = SAT_SCI
    model.data[0, 2:, 1, 3] = SAT_SCI

    # CR (C): group 2, on edge
    model.data[0, 0, 4, 0] = 20.0
    model.data[0, 1, 4, 0] = 35.0
    model.data[0, 2, 4, 0] = 260.0
    model.data[0, 3, 4, 0] = 275.0
    model.data[0, 4, 4, 0] = 295.0
    model.data[0, 5, 4, 0] = 310.0
    model.data[0, 6, 4, 0] = 325.0
    model.data[0, 7, 4, 0] = 340.0
    # ... and set 1 of the 3 neighbors to SAT, and adjust its ramp
    model.groupdq[0, 2:, 3, 0] = SATURATED
    model.data[0, 2:, 3, 0] = SAT_SCI

    # CR (D): group 2, not on edge
    model.data[0, 0, 4, 5] = 60.0
    model.data[0, 1, 4, 5] = 75.0
    model.data[0, 2, 4, 5] = 160.0
    model.data[0, 3, 4, 5] = 175.0
    model.data[0, 4, 4, 5] = 195.0
    model.data[0, 5, 4, 5] = 210.0
    model.data[0, 6, 4, 5] = 225.0
    model.data[0, 7, 4, 5] = 240.0
    # ... and set 2 of the 4 neighbors to SAT, and adjust their ramps
    model.groupdq[0, 2:, 4, 4] = SATURATED
    model.groupdq[0, 2:, 4, 6] = SATURATED
    model.data[0, 2:, 4, 4] = SAT_SCI
    model.data[0, 2:, 4, 6] = SAT_SCI

    # CR (E): group 4, not on edge
    model.data[0, 0, 2, 5] = 20.0
    model.data[0, 1, 2, 5] = 25.0
    model.data[0, 2, 2, 5] = 31.0
    model.data[0, 3, 2, 5] = 36.0
    model.data[0, 4, 2, 5] = 190.0
    model.data[0, 5, 2, 5] = 196.0
    model.data[0, 6, 2, 5] = 201.0
    model.data[0, 7, 2, 5] = 207.0
    # ... and set 2 of the 4 neighbors to SAT, and adjust their ramps
    model.groupdq[0, 4:, 1, 5] = SATURATED
    model.groupdq[0, 4:, 3, 5] = SATURATED
    model.data[0, 4:, 1, 5] = SAT_SCI
    model.data[0, 4:, 3, 5] = SAT_SCI

    # CR (F): group 4, on edge
    model.data[0, 0, 5, 2] = 50.0
    model.data[0, 1, 5, 2] = 55.0
    model.data[0, 2, 5, 2] = 61.0
    model.data[0, 3, 5, 2] = 66.0
    model.data[0, 4, 5, 2] = 290.0
    model.data[0, 5, 5, 2] = 296.0
    model.data[0, 6, 5, 2] = 301.0
    model.data[0, 7, 5, 2] = 307.0
    # ... and set 1 of the 3 neighbors to SAT, and adjust its ramp
    model.groupdq[0, 4:, 5, 1] = SATURATED
    model.data[0, 4:, 5, 1] = SAT_SCI

    out_model = run_detect_jumps(model, gain, rnoise, 4.0, 5.0, 6.0, 'none', 200, 4, True)
    dq_out = out_model.groupdq

    # Do assertions for each cosmic ray and its neighbors
    # CR (A):
    # Check the CR-impacted pixel and non-SAT neighbor
    assert dq_out[0, 2, 0, 0] == JUMP_DET
    assert dq_out[0, 2, 1, 0] == JUMP_DET
    # Check the SAT neighbor
    assert dq_out[0, 2, 0, 1] == SATURATED

    # CR (B):
    # Check the CR-impacted pixel and 2 non-SAT neighbors
    assert dq_out[0, 2, 1, 4] == JUMP_DET
    assert dq_out[0, 2, 1, 5] == JUMP_DET
    assert dq_out[0, 2, 2, 4] == JUMP_DET
    # Check the 2 SAT neighbors
    assert dq_out[0, 2, 0, 4] == SATURATED
    assert dq_out[0, 2, 1, 3] == SATURATED

    # CR (C):
    # Check the CR-impacted pixel and 2 non-SAT neighbors
    assert dq_out[0, 2, 4, 0] == JUMP_DET
    assert dq_out[0, 2, 4, 1] == JUMP_DET
    assert dq_out[0, 2, 5, 0] == JUMP_DET
    # Check the SAT neighbor
    assert dq_out[0, 2, 3, 0] == SATURATED

    # CR (D):
    # Check the CR-impacted pixel and 2 non-SAT neighbors
    assert dq_out[0, 2, 4, 5] == JUMP_DET
    assert dq_out[0, 2, 3, 5] == JUMP_DET
    assert dq_out[0, 2, 5, 5] == JUMP_DET
    # Check the 2 SAT neighbors
    assert dq_out[0, 2, 4, 4] == SATURATED
    assert dq_out[0, 2, 4, 6] == SATURATED

    # CR (E):
    # Check the CR-impacted pixel and 2 non-SAT neighbors
    assert dq_out[0, 4, 2, 5] == JUMP_DET
    assert dq_out[0, 4, 2, 4] == JUMP_DET
    assert dq_out[0, 4, 2, 6] == JUMP_DET
    # Check the 2 SAT neighbors
    assert dq_out[0, 4, 1, 5] == SATURATED
    assert dq_out[0, 4, 3, 5] == SATURATED

    # CR (F):
    # Check the CR-impacted pixel and 2 non-SAT neighbors
    assert dq_out[0, 4, 5, 2] == JUMP_DET
    assert dq_out[0, 4, 4, 2] == JUMP_DET
    assert dq_out[0, 4, 5, 3] == JUMP_DET
    # Check the 1 SAT neighbor
    assert dq_out[0, 4, 5, 1] == SATURATED


@pytest.fixture
def setup_inputs():
    def _setup(ngroups=10, readnoise=10, nints=1, nrows=102, ncols=103,
               nframes=1, grouptime=1.0, gain=1, deltatime=1, subarray=False):

        if subarray:
            shape = (nints, ngroups, 20, 20)
        else:
            shape = (nints, ngroups, nrows, ncols)

        model = RampModel(shape)
        model.meta.instrument.name = 'MIRI'
        model.meta.instrument.detector = 'MIRIMAGE'
        model.meta.instrument.filter = 'F480M'
        model.meta.observation.date = '2015-10-13'
        model.meta.exposure.type = 'MIR_IMAGE'
        model.meta.exposure.group_time = deltatime
        model.meta.subarray.name = 'FULL'
        model.meta.subarray.xstart = 1
        model.meta.subarray.ystart = 1
        model.meta.subarray.xsize = shape[3]
        model.meta.subarray.ysize = shape[2]
        model.meta.exposure.frame_time = deltatime
        model.meta.exposure.ngroups = ngroups
        model.meta.exposure.group_time = deltatime
        model.meta.exposure.nframes = nframes
        model.meta.exposure.groupgap = 0

        gain_model = GainModel((nrows, ncols))
        gain_model.data += gain
        gain_model.meta.instrument.name = 'MIRI'
        gain_model.meta.subarray.xstart = 1
        gain_model.meta.subarray.ystart = 1
        gain_model.meta.subarray.xsize = ncols
        gain_model.meta.subarray.ysize = nrows

        read_noise_model = ReadnoiseModel((nrows, ncols))
        read_noise_model.data += readnoise
        read_noise_model.meta.instrument.name = 'MIRI'
        read_noise_model.meta.subarray.xstart = 1
        read_noise_model.meta.subarray.ystart = 1
        read_noise_model.meta.subarray.xsize = ncols
        read_noise_model.meta.subarray.ysize = nrows

        return model, read_noise_model, gain_model

    return _setup
