import pytest
import numpy as np

from jwst.jump.twopoint_difference import find_crs
from jwst.jump.twopoint_difference import get_clipped_median_vector
from jwst.jump.twopoint_difference import get_clipped_median_array
from jwst.datamodels import dqflags


def test_nocrs_noflux(setup_cube):
    ngroups = 5
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups)
    out_gdq, row_below_gdq, row_above_gdq = find_crs(data, gdq, read_noise, rej_threshold, rej_threshold, rej_threshold,
                                                     nframes, False, 200, 10)
    assert 0 == np.max(out_gdq)  # no CR found


def test_5grps_cr3_noflux(setup_cube):
    """"
       A test that has a jump in the third group of pixel 100, 100
       """
    ngroups = 5
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups)

    data[0, 0:2, 100, 100] = 10.0
    data[0, 2:5, 100, 100] = 1000
    out_gdq, row_below_gdq, row_above_gdq = find_crs(data, gdq, read_noise, rej_threshold, rej_threshold, rej_threshold,
                                                     nframes, False, 200, 10)
    assert 4 == np.max(out_gdq)  # a CR was found
    assert 2 == np.argmax(out_gdq[0, :, 100, 100])  # find the CR in the expected group


def test_5grps_cr2_noflux(setup_cube):
    ngroups = 5
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups)

    data[0, 0, 100, 100] = 10.0
    data[0, 1:6, 100, 100] = 1000
    out_gdq, row_below_gdq, row_above_gdq = find_crs(data, gdq, read_noise, rej_threshold, rej_threshold, rej_threshold,
                                                     nframes, False, 200, 10)
    assert 4 == np.max(out_gdq)  # a CR was found
    assert 1 == np.argmax(out_gdq[0, :, 100, 100])  # find the CR in the expected group


def test_6grps_negative_differences_zeromedian(setup_cube):
    ngroups = 6
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups)

    data[0, 0, 100, 100] = 100
    data[0, 1, 100, 100] = 90
    data[0, 2, 100, 100] = 95
    data[0, 3, 100, 100] = 105
    data[0, 4, 100, 100] = 100
    data[0, 5, 100, 100] = 100
    out_gdq, row_below_gdq, row_above_gdq = find_crs(data, gdq, read_noise, rej_threshold, rej_threshold, rej_threshold,
                                                     nframes, False, 200, 10)
    assert 0 == np.max(out_gdq)  # no CR was found


def test_5grps_cr2_negjumpflux(setup_cube):
    ngroups = 5
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups)

    data[0, 0, 100, 100] = 1000.0
    data[0, 1:6, 100, 100] = 10
    out_gdq, row_below_gdq, row_above_gdq = find_crs(data, gdq, read_noise, rej_threshold, rej_threshold, rej_threshold,
                                                     nframes, False, 200, 10)
    assert 4 == np.max(out_gdq)  # a CR was found
    assert 1 == np.argmax(out_gdq[0, :, 100, 100])  # find the CR in the expected group


def test_3grps_cr2_noflux(setup_cube):
    ngroups = 3
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups)
    data[0, 0, 100, 100] = 10.0
    data[0, 1:4, 100, 100] = 1000
    out_gdq, row_below_gdq, row_above_gdq = find_crs(data, gdq, read_noise, rej_threshold, rej_threshold, rej_threshold,
                                                     nframes, False, 200, 10)
    assert 4 == np.max(out_gdq)  # a CR was found
    assert np.array_equal([0, 4, 0], out_gdq[0, :, 100, 100])


@pytest.mark.xfail
def test_4grps_cr2_noflux(setup_cube):
    # This test should fail because with 2 CRs and only 4 groups we cannot detect the jump
    ngroups = 4
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups)
    data[0, 0, 100, 100] = 10.0
    data[0, 1:4, 100, 100] = 1000
    out_gdq, row_below_gdq, row_above_gdq = find_crs(data, gdq, read_noise, rej_threshold, rej_threshold, rej_threshold,
                                                     nframes, False, 200, 10)
    assert 4 == np.max(out_gdq)  # a CR was found
    assert 1 == np.argmax(out_gdq[0, :, 100, 100])  # find the CR in the expected group


def test_5grps_cr2_nframe2(setup_cube):
    ngroups = 5
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups)
    nframes = 2
    data[0, 0, 100, 100] = 10.0
    data[0, 1, 100, 100] = 500
    data[0, 2, 100, 100] = 1002
    data[0, 3, 100, 100] = 1001
    data[0, 4, 100, 100] = 1005
    out_gdq, row_below_gdq, row_above_gdq = find_crs(data, gdq, read_noise, rej_threshold, rej_threshold, rej_threshold,
                                                     nframes, False, 200, 10)
    assert 4 == np.max(out_gdq)  # a CR was found
    assert np.array_equal([0, 4, 4, 0, 0], out_gdq[0, :, 100, 100])


@pytest.mark.xfail
def test_4grps_twocrs_2nd_4th(setup_cube):
    # This test should fail because two jumps with four groups we cannot find the second jump.
    ngroups = 4
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups)
    nframes = 1
    data[0, 0, 100, 100] = 10.0
    data[0, 1, 100, 100] = 60
    data[0, 2, 100, 100] = 60
    data[0, 3, 100, 100] = 115
    out_gdq, row_below_gdq, row_above_gdq = find_crs(data, gdq, read_noise, rej_threshold, rej_threshold, rej_threshold,
                                                     nframes, False, 200, 10)
    assert np.max(out_gdq) == 4  # a CR was found


def test_5grps_twocrs_2nd_5th(setup_cube):
    ngroups = 5
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups)
    nframes = 1
    data[0, 0, 100, 100] = 10.0
    data[0, 1, 100, 100] = 60
    data[0, 2, 100, 100] = 60
    data[0, 3, 100, 100] = 60
    data[0, 4, 100, 100] = 115
    out_gdq, row_below_gdq, row_above_gdq = find_crs(data, gdq, read_noise, rej_threshold, rej_threshold, rej_threshold,
                                                     nframes, False, 200, 10)
    assert 4 == np.max(out_gdq)  # a CR was found
    assert np.array_equal([0, 4, 0, 0, 4],out_gdq[0, :, 100, 100])


def test_5grps_twocrs_2nd_5thbig(setup_cube):
    ngroups = 5
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups)
    nframes = 1
    data[0, 0, 100, 100] = 10.0
    data[0, 1, 100, 100] = 60
    data[0, 2, 100, 100] = 60
    data[0, 3, 100, 100] = 60
    data[0, 4, 100, 100] = 2115
    out_gdq, row_below_gdq, row_above_gdq = find_crs(data, gdq, read_noise, rej_threshold, rej_threshold, rej_threshold,
                                                     nframes, False, 200, 10)
    assert 4 == np.max(out_gdq)  # a CR was found
    assert np.array_equal([0, 4, 0, 0, 4], out_gdq[0, :, 100, 100])


def test_10grps_twocrs_2nd_8th_big(setup_cube):
    ngroups = 10
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups)
    nframes = 1
    data[0, 0, 100, 100] = 10.0
    data[0, 1, 100, 100] = 60
    data[0, 2, 100, 100] = 60
    data[0, 3, 100, 100] = 60
    data[0, 4, 100, 100] = 60
    data[0, 5, 100, 100] = 60
    data[0, 6, 100, 100] = 60
    data[0, 7, 100, 100] = 2115
    data[0, 8, 100, 100] = 2115
    data[0, 9, 100, 100] = 2115
    out_gdq, row_below_gdq, row_above_gdq = find_crs(data, gdq, read_noise, rej_threshold, rej_threshold, rej_threshold,
                                                     nframes, False, 200, 10)
    assert 4 == np.max(out_gdq)  # a CR was found
    assert np.array_equal([0, 4, 0, 0, 0, 0, 0, 4, 0, 0], out_gdq[0, :, 100, 100])


def test_10grps_twocrs_10percenthit(setup_cube):
    ngroups = 10
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups)
    nframes = 2
    data[0:200, 0, 100, 100] = 10.0
    data[0:200, 1, 100, 100] = 60
    data[0:200, 2, 100, 100] = 60
    data[0:200, 3, 100, 100] = 60
    data[0:200, 4, 100, 100] = 60
    data[0:200, 5, 100, 100] = 60
    data[0:200, 6, 100, 100] = 60
    data[0:200, 7, 100, 100] = 2115
    data[0:200, 8, 100, 100] = 2115
    data[0:200, 9, 100, 100] = 2115
    out_gdq, row_below_gdq, row_above_gdq = find_crs(data, gdq, read_noise, rej_threshold, rej_threshold, rej_threshold,
                                                     nframes, False, 200, 10)
    assert 4 == np.max(out_gdq)  # a CR was found
    assert np.array_equal([0, 4, 0, 0, 0, 0, 0, 4, 0, 0], out_gdq[0, :, 100, 100])


def test_5grps_twocrs_2nd_5thbig_nframes2(setup_cube):
    ngroups = 5
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups, readnoise=10 * np.sqrt(2))
    nframes = 2
    data[0, 0, 100, 100] = 10.0
    data[0, 1, 100, 100] = 60
    data[0, 2, 100, 100] = 60
    data[0, 3, 100, 100] = 60
    data[0, 4, 100, 100] = 2115
    out_gdq, row_below_gdq, row_above_gdq = find_crs(data, gdq, read_noise, rej_threshold, rej_threshold, rej_threshold,
                                                     nframes, False, 200, 10)
    assert 4 == np.max(out_gdq)  # a CR was found
    assert np.array_equal([0, 4, 0, 0, 4], out_gdq[0, :, 100, 100])


def test_6grps_twocrs_2nd_5th(setup_cube):
    ngroups = 6
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups)
    nframes = 1
    data[0, 0, 100, 100] = 10.0
    data[0, 1, 100, 100] = 60
    data[0, 2, 100, 100] = 60
    data[0, 3, 100, 100] = 60
    data[0, 4, 100, 100] = 115
    data[0, 5, 100, 100] = 115
    out_gdq, row_below_gdq, row_above_gdq = find_crs(data, gdq, read_noise, rej_threshold, rej_threshold, rej_threshold,
                                                     nframes, False, 200, 10)
    assert 4 == np.max(out_gdq)  # a CR was found
    assert np.array_equal([0, 4, 0, 0, 4, 0], out_gdq[0, :, 100, 100])


def test_6grps_twocrs_2nd_5th_nframes2(setup_cube):
    ngroups = 6
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups, readnoise=10 * np.sqrt(2))
    nframes = 2
    data[0, 0, 100, 100] = 10.0
    data[0, 1, 100, 100] = 60
    data[0, 2, 100, 100] = 60
    data[0, 3, 100, 100] = 60
    data[0, 4, 100, 100] = 115
    data[0, 5, 100, 100] = 115
    out_gdq, row_below_gdq, row_above_gdq = find_crs(data, gdq, read_noise, rej_threshold, rej_threshold, rej_threshold,
                                                     nframes, False, 200, 10)
    assert 4 == np.max(out_gdq)  # a CR was found
    assert np.array_equal([0, 4, 0, 0, 4, 0], out_gdq[0, :, 100, 100])


def test_6grps_twocrs_twopixels_nframes2(setup_cube):
    ngroups = 6
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups, readnoise=10 * np.sqrt(2))
    nframes = 2
    data[0, 0, 100, 100] = 10.0
    data[0, 1, 100, 100] = 60
    data[0, 2, 100, 100] = 60
    data[0, 3, 100, 100] = 60
    data[0, 4, 100, 100] = 115
    data[0, 5, 100, 100] = 115
    data[0, 0, 200, 100] = 10.0
    data[0, 1, 200, 100] = 10.0
    data[0, 2, 200, 100] = 60
    data[0, 3, 200, 100] = 60
    data[0, 4, 200, 100] = 115
    data[0, 5, 200, 100] = 115
    out_gdq, row_below_gdq, row_above_gdq = find_crs(data, gdq, read_noise, rej_threshold, rej_threshold, rej_threshold,
                                                     nframes, False, 200, 10)
    assert 4 == np.max(out_gdq)  # a CR was found
    assert np.array_equal([0, 4, 0, 0, 4, 0], out_gdq[0, :, 100, 100])
    assert np.array_equal([0, 0, 4, 0, 4, 0], out_gdq[0, :, 200, 100])


def test_5grps_cr2_negslope(setup_cube):
    ngroups = 5
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups)
    nframes = 1
    data[0, 0, 100, 100] = 100.0
    data[0, 1, 100, 100] = 0
    data[0, 2, 100, 100] = -200
    data[0, 3, 100, 100] = -260
    data[0, 4, 100, 100] = -360
    out_gdq, row_below_gdq, row_above_gdq = find_crs(data, gdq, read_noise, rej_threshold, rej_threshold, rej_threshold,
                                                     nframes, False, 200, 10)
    assert 4 == np.max(out_gdq)  # a CR was found
    assert np.array_equal([0, 0, 4, 0, 0], out_gdq[0, :, 100, 100])


def test_6grps_1cr(setup_cube):
    ngroups = 6
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups, readnoise=10)
    nframes = 1
    data[0, 0, 100, 100] = 0
    data[0, 1, 100, 100] = 10
    data[0, 2, 100, 100] = 21
    data[0, 3, 100, 100] = 33
    data[0, 4, 100, 100] = 46
    data[0, 5, 100, 100] = 1146
    out_gdq, row_below_gdq, row_above_gdq = find_crs(data, gdq, read_noise, rej_threshold, rej_threshold, rej_threshold,
                                                     nframes, False, 200, 10)
    assert 4 == out_gdq[0, 5, 100, 100]


def test_7grps_1cr(setup_cube):
    ngroups = 7
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups, readnoise=10)
    nframes = 1
    data[0, 0, 100, 100] = 0
    data[0, 1, 100, 100] = 10
    data[0, 2, 100, 100] = 21
    data[0, 3, 100, 100] = 33
    data[0, 4, 100, 100] = 46
    data[0, 5, 100, 100] = 60
    data[0, 6, 100, 100] = 1160
    out_gdq, row_below_gdq, row_above_gdq = find_crs(data, gdq, read_noise, rej_threshold, rej_threshold, rej_threshold,
                                                     nframes, False, 200, 10)
    assert 4 == out_gdq[0, 6, 100, 100]


def test_8grps_1cr(setup_cube):
    ngroups = 8
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups, readnoise=10)
    nframes = 1
    data[0, 0, 100, 100] = 0
    data[0, 1, 100, 100] = 10
    data[0, 2, 100, 100] = 21
    data[0, 3, 100, 100] = 33
    data[0, 4, 100, 100] = 46
    data[0, 5, 100, 100] = 60
    data[0, 6, 100, 100] = 1160
    data[0, 7, 100, 100] = 1175
    out_gdq, row_below_gdq, row_above_gdq = find_crs(data, gdq, read_noise, rej_threshold, rej_threshold, rej_threshold,
                                                     nframes, False, 200, 10)
    assert 4 == out_gdq[0, 6, 100, 100]


def test_9grps_1cr_1sat(setup_cube):
    ngroups = 9
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups, readnoise=10)
    nframes = 1
    data[0, 0, 100, 100] = 0
    data[0, 1, 100, 100] = 10
    data[0, 2, 100, 100] = 21
    data[0, 3, 100, 100] = 33
    data[0, 4, 100, 100] = 46
    data[0, 5, 100, 100] = 60
    data[0, 6, 100, 100] = 1160
    data[0, 7, 100, 100] = 1175
    data[0, 8, 100, 100] = 6175
    gdq[0, 8, 100, 100] = dqflags.group['SATURATED']
    out_gdq, row_below_gdq, row_above_gdq = find_crs(data, gdq, read_noise, rej_threshold, rej_threshold, rej_threshold,
                                                     nframes, False, 200, 10)
    assert 4 == out_gdq[0, 6, 100, 100]


def test_10grps_1cr_2sat(setup_cube):
    ngroups = 10
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups, readnoise=10)
    nframes = 1
    data[0, 0, 100, 100] = 0
    data[0, 1, 100, 100] = 10
    data[0, 2, 100, 100] = 21
    data[0, 3, 100, 100] = 33
    data[0, 4, 100, 100] = 46
    data[0, 5, 100, 100] = 60
    data[0, 6, 100, 100] = 1160
    data[0, 7, 100, 100] = 1175
    data[0, 8, 100, 100] = 6175
    data[0, 9, 100, 100] = 6175
    gdq[0, 8, 100, 100] = dqflags.group['SATURATED']
    gdq[0, 9, 100, 100] = dqflags.group['SATURATED']
    out_gdq, row_below_gdq, row_above_gdq = find_crs(data, gdq, read_noise, rej_threshold, rej_threshold, rej_threshold,
                                                     nframes, False, 200, 10)
    assert 4 == out_gdq[0, 6, 100, 100]


def test_11grps_1cr_3sat(setup_cube):
    ngroups = 11
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups, readnoise=10)
    nframes = 1
    data[0, 0, 100, 100] = 0
    data[0, 1, 100, 100] = 20
    data[0, 2, 100, 100] = 39
    data[0, 3, 100, 100] = 57
    data[0, 4, 100, 100] = 74
    data[0, 5, 100, 100] = 90
    data[0, 6, 100, 100] = 1160
    data[0, 7, 100, 100] = 1175
    data[0, 8, 100, 100] = 6175
    data[0, 9, 100, 100] = 6175
    data[0, 10, 100, 100] = 6175
    gdq[0, 8, 100, 100] = dqflags.group['SATURATED']
    gdq[0, 9, 100, 100] = dqflags.group['SATURATED']
    gdq[0, 10, 100, 100] = dqflags.group['SATURATED']
    out_gdq, row_below_gdq, row_above_gdq = find_crs(data, gdq, read_noise, rej_threshold, rej_threshold, rej_threshold,
                                                     nframes, False, 200, 10)
    assert 4 == out_gdq[0, 6, 100, 100]


def test_11grps_0cr_3donotuse(setup_cube):
    ngroups = 11
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups, readnoise=10)
    nframes = 1
    data[0, 0, 100, 100] = 0
    data[0, 1, 100, 100] = 18
    data[0, 2, 100, 100] = 39
    data[0, 3, 100, 100] = 57
    data[0, 4, 100, 100] = 74
    data[0, 5, 100, 100] = 90
    data[0, 6, 100, 100] = 115
    data[0, 7, 100, 100] = 131
    data[0, 8, 100, 100] = 150
    data[0, 9, 100, 100] = 6175
    data[0, 10, 100, 100] = 6175
    gdq[0, 0, 100, 100] = dqflags.group['DO_NOT_USE']
    gdq[0, 9, 100, 100] = dqflags.group['DO_NOT_USE']
    gdq[0, 10, 100, 100] = dqflags.group['DO_NOT_USE']
    out_gdq, row_below_gdq, row_above_gdq = find_crs(data, gdq, read_noise, rej_threshold, rej_threshold, rej_threshold,
                                                     nframes, False, 200, 10)
    assert np.array_equal([0, 0, 0, 0, 0, 0, 0, 0], out_gdq[0, 1:-2, 100, 100])


def test_10grps_cr2_gt3sigma(setup_cube):
    ngroups = 10
    crmag = 16
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups, readnoise=5)
    nframes = 1
    data[0, 0, 100, 100] = 0
    data[0, 1:11, 100, 100] = crmag
    out_gdq, row_below_gdq, row_above_gdq = find_crs(data, gdq, read_noise, rej_threshold, rej_threshold, rej_threshold,
                                                     nframes, False, 200, 10)
    assert 4 == np.max(out_gdq)  # a CR was found
    assert np.array_equal([0, 4, 0, 0, 0, 0, 0, 0, 0, 0], out_gdq[0, :, 100, 100])


def test_10grps_cr2_3sigma_nocr(setup_cube):
    ngroups = 10
    crmag = 15
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups, readnoise=5)
    nframes = 1
    data[0, 0, 100, 100] = 0
    data[0, 1:11, 100, 100] = crmag
    out_gdq, row_below_gdq, row_above_gdq = find_crs(data, gdq, read_noise, rej_threshold, rej_threshold, rej_threshold,
                                                     nframes, False, 200, 10)
    assert 0 == np.max(out_gdq)  # a CR was found
    assert np.array_equal([0, 0, 0, 0, 0, 0, 0, 0, 0, 0], out_gdq[0, :, 100, 100])


def test_10grps_cr2_gt3sigma_2frames(setup_cube):
    ngroups = 10
    crmag = 16
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups, readnoise=5 * np.sqrt(2))
    nframes = 2
    data[0, 0, 100, 100] = 0
    data[0, 1:11, 100, 100] = crmag
    out_gdq, row_below_gdq, row_above_gdq = find_crs(data, gdq, read_noise, rej_threshold, rej_threshold, rej_threshold,
                                                     nframes, False, 200, 10)
    assert 4 == np.max(out_gdq)  # a CR was found
    assert np.array_equal([0, 4, 0, 0, 0, 0, 0, 0, 0, 0], out_gdq[0, :, 100, 100])


def test_10grps_cr2_gt3sigma_2frames_offdiag(setup_cube):
    ngroups = 10
    crmag = 16
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups, readnoise=5 * np.sqrt(2))
    nframes = 2
    data[0, 0, 100, 110] = 0
    data[0, 1:11, 100, 110] = crmag
    out_gdq, row_below_gdq, row_above_gdq = find_crs(data, gdq, read_noise, rej_threshold, rej_threshold, rej_threshold,
                                                     nframes, False, 200, 10)
    assert 4 == np.max(out_gdq)  # a CR was found
    assert np.array_equal([0, 4, 0, 0, 0, 0, 0, 0, 0, 0], out_gdq[0, :, 100, 110])


def test_10grps_cr2_3sigma_2frames_nocr(setup_cube):
    ngroups = 10
    crmag = 15
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups, readnoise=5 * np.sqrt(2))
    nframes = 2
    data[0, 0, 100, 100] = 0
    data[0, 1:11, 100, 100] = crmag
    out_gdq, row_below_gdq, row_above_gdq = find_crs(data, gdq, read_noise, rej_threshold, rej_threshold, rej_threshold,
                                                     nframes, False, 200, 10)
    assert 0 == np.max(out_gdq)  # a CR was found
    assert np.array_equal([0, 0, 0, 0, 0, 0, 0, 0, 0, 0], out_gdq[0, :, 100, 100])


def test_10grps_nocr_2pixels_sigma0(setup_cube):
    ngroups = 10
    crmag = 15
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups, readnoise=5 * np.sqrt(2))
    nframes = 1
    data[0, 0, 100, 100] = crmag
    data[0, 1:11, 100, 100] = crmag
    read_noise[50, 50] = 0.0
    read_noise[60, 60] = 0.0
    out_gdq, row_below_gdq, row_above_gdq = find_crs(data, gdq, read_noise, rej_threshold, rej_threshold, rej_threshold,
                                                     nframes, False, 200, 10)
    assert 0 == np.max(out_gdq)  # no CR was found


def test_5grps_satat4_crat3(setup_cube):
    ngroups = 5
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups, readnoise=5 * np.sqrt(2))
    nframes = 1
    data[0, 0, 100, 100] = 10000
    data[0, 1, 100, 100] = 30000
    data[0, 2, 100, 100] = 60000
    data[0, 3, 100, 100] = 61000
    data[0, 4, 100, 100] = 61000
    gdq[0, 3, 100, 100] = dqflags.group['SATURATED']
    gdq[0, 4, 100, 100] = dqflags.group['SATURATED']
    out_gdq, row_below_gdq, row_above_gdq = find_crs(data, gdq, read_noise, rej_threshold, rej_threshold, rej_threshold,
                                                     nframes, False, 200, 10)
    # assert 4 == np.max(out_gdq)  # no CR was found
    assert np.array_equal(
        [0, 0, dqflags.group['JUMP_DET'], dqflags.group['SATURATED'], dqflags.group['SATURATED']],
        out_gdq[0, :, 100, 100]
    )


def test_6grps_satat6_crat1(setup_cube):
    ngroups = 6
    # crmag = 1000
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups, readnoise=5 * np.sqrt(2))
    nframes = 1
    data[0, 0, 100, 100] = 10000
    data[0, 1, 100, 100] = 35000  # CR
    data[0, 2, 100, 100] = 40005
    data[0, 3, 100, 100] = 45029
    data[0, 4, 100, 100] = 50014
    data[0, 5, 100, 101] = 61000
    data[0, 0, 100, 101] = 10000
    data[0, 1, 100, 101] = 15001
    data[0, 2, 100, 101] = 20003
    data[0, 3, 100, 101] = 25006
    data[0, 4, 100, 101] = 30010
    data[0, 5, 100, 101] = 35015
    gdq[0, 5, 100, 100] = dqflags.group['SATURATED']
    out_gdq, row_below_gdq, row_above_gdq = find_crs(data, gdq, read_noise, rej_threshold, rej_threshold, rej_threshold,
                                                     nframes, False, 200, 10)
    assert np.array_equal([0, dqflags.group['JUMP_DET'], 0, 0, 0, dqflags.group['SATURATED']], out_gdq[0, :, 100, 100])


def test_6grps_satat6_crat1_flagadjpixelsoff(setup_cube):
    # Check that neighbor pixels are not flagged when flag_4_neighbors is set to false
    ngroups = 6
    # crmag = 1000
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups, readnoise=5 * np.sqrt(2))
    nframes = 1
    data[0, 0, 100, 100] = 10000
    data[0, 1, 100, 100] = 35000  # CR
    data[0, 2, 100, 100] = 40005
    data[0, 3, 100, 100] = 45029
    data[0, 4, 100, 100] = 50014
    data[0, 5, 100, 101] = 61000
    data[0, 0, 100, 101] = 10000
    data[0, 1, 100, 101] = 15001
    data[0, 2, 100, 101] = 20003
    data[0, 3, 100, 101] = 25006
    data[0, 4, 100, 101] = 30010
    data[0, 5, 100, 101] = 35015
    gdq[0, 5, 100, 100] = dqflags.group['SATURATED']
    out_gdq, row_below_gdq, row_above_gdq = find_crs(data, gdq, read_noise, rej_threshold, rej_threshold, rej_threshold,
                                                     nframes, False, 200, 10)
    assert np.array_equal([0, dqflags.group['JUMP_DET'], 0, 0, 0, dqflags.group['SATURATED']], out_gdq[0, :, 100, 100])
    assert np.array_equal([0, 0, 0, 0, 0, 0], out_gdq[0, :, 99, 100])


def test_10grps_satat8_crsat3and6(setup_cube):
    ngroups = 10
    # crmag = 1000
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups, readnoise=5 * np.sqrt(2))
    nframes = 1
    data[0, 0, 100, 100] = 0
    data[0, 1, 100, 100] = 5000
    data[0, 2, 100, 100] = 15000  # CR
    data[0, 3, 100, 100] = 20000
    data[0, 4, 100, 100] = 25000
    data[0, 5, 100, 100] = 40000  # CR
    data[0, 6, 100, 100] = 45000
    data[0, 7:11, 100, 100] = 61000
    gdq[0, 7:11, 100, 100] = dqflags.group['SATURATED']
    out_gdq, row_below_gdq, row_above_gdq = find_crs(data, gdq, read_noise, rej_threshold, rej_threshold, rej_threshold,
                                                     nframes, False, 200, 10)

    # assert 4 == np.max(out_gdq)  # no CR was found
    assert np.array_equal(
        [0, 0, dqflags.group['JUMP_DET'], 0, 0, dqflags.group['JUMP_DET'], 0,
            dqflags.group['SATURATED'], dqflags.group['SATURATED'], dqflags.group['SATURATED']],
        out_gdq[0, :, 100, 100])


def test_median_with_saturation(setup_cube):
    ngroups = 10
    # crmag = 1000
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups, readnoise=5 * np.sqrt(2))
    nframes = 1
    data[0, 0, 100, 100] = 0
    data[0, 1, 100, 100] = 4500
    data[0, 2, 100, 100] = 9100
    data[0, 3, 100, 100] = 13800
    data[0, 4, 100, 100] = 18600
    data[0, 5, 100, 100] = 40000  # CR
    data[0, 6, 100, 100] = 44850
    data[0, 7, 100, 100] = 49900
    data[0, 8:10, 100, 100] = 60000
    gdq[0, 7:10, 100, 100] = dqflags.group['SATURATED']
    out_gdq, row_below_gdq, row_above_gdq = find_crs(data, gdq, read_noise, rej_threshold, rej_threshold, rej_threshold,
                                                     nframes, False, 200, 10)
    assert np.array_equal([0, 0, 0, 0, 0, 4, 0, 2, 2, 2], out_gdq[0, :, 100, 100])


def test_median_with_saturation_even_num_sat_frames(setup_cube):
    ngroups = 10
    # crmag = 1000
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups, readnoise=5 * np.sqrt(2))
    nframes = 1
    data[0, 0, 100, 100] = 0
    data[0, 1, 100, 100] = 4500
    data[0, 2, 100, 100] = 9100
    data[0, 3, 100, 100] = 13800
    data[0, 4, 100, 100] = 18600
    data[0, 5, 100, 100] = 40000  # CR
    data[0, 6, 100, 100] = 44850
    data[0, 7, 100, 100] = 49900
    data[0, 8:10, 100, 100] = 60000
    gdq[0, 6:10, 100, 100] = dqflags.group['SATURATED']
    out_gdq, row_below_gdq, row_above_gdq = find_crs(data, gdq, read_noise, rej_threshold, rej_threshold, rej_threshold,
                                                     nframes, False, 200, 10)
    assert np.array_equal([0, 0, 0, 0, 0, 4, 2, 2, 2, 2], out_gdq[0, :, 100, 100])


def test_median_with_saturation_odd_number_final_difference(setup_cube):
    ngroups = 9
    # crmag = 1000
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups, readnoise=5 * np.sqrt(2))
    nframes = 1
    data[0, 0, 100, 100] = 0
    data[0, 1, 100, 100] = 4500
    data[0, 2, 100, 100] = 9100
    data[0, 3, 100, 100] = 13800
    data[0, 4, 100, 100] = 18600
    data[0, 5, 100, 100] = 40000  # CR
    data[0, 6, 100, 100] = 44850
    data[0, 7, 100, 100] = 49900
    data[0, 8:9, 100, 100] = 60000
    gdq[0, 6:9, 100, 100] = dqflags.group['SATURATED']

    out_gdq, row_below_gdq, row_above_gdq = find_crs(data, gdq, read_noise, rej_threshold, rej_threshold, rej_threshold,
                                                     nframes, False, 200, 10)
    assert np.array_equal([0, 0, 0, 0, 0, 4, 2, 2, 2], out_gdq[0, :, 100, 100])


def test_first_last_group(setup_cube):
    ngroups = 7
    nframes = 1
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups, readnoise=25.0)

    #  set up the data so that if the first and last group are used in jump
    #  detection it would cause a jump to be detected between group 1-2
    #  and group 6-7. Add a jump between 3 and 4 just to make sure jump detection is working
    #  set group 1 to be 10,000
    data[0, 0, 100, 100] = 10000.0
    #  set groups 1, 2 - to be around 30,000
    data[0, 1, 100, 100] = 30000.0
    data[0, 2, 100, 100] = 30020.0
    #  set up a jump to make sure it is detected
    data[0, 3, 100, 100] = 40000.0
    data[0, 4, 100, 100] = 40020.0
    data[0, 5, 100, 100] = 40040.0
    #  set group 6 to be 50,000
    data[0, 6, 100, 100] = 50000.0

    gdq[0, 0, 100, 100] = dqflags.group['DO_NOT_USE']
    gdq[0, 6, 100, 100] = dqflags.group['DO_NOT_USE']
    outgdq, row_below_gdq, row_above_gdq = find_crs(data, gdq, read_noise, rej_threshold, rej_threshold, rej_threshold,
                                                    nframes, False, 200, 10)

    assert outgdq[0, 0, 100, 100] == dqflags.group['DO_NOT_USE']
    assert outgdq[0, 6, 100, 100] == dqflags.group['DO_NOT_USE']
    assert outgdq[0, 3, 100, 100] == dqflags.group['JUMP_DET']


def test_4grps_1cr(setup_cube):
    """"
    A test of a four group integration that has a jump in the third group of pixel 100, 100
    """
    ngroups = 4
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups, readnoise=5 * np.sqrt(2))
    nframes = 1
    data[0, 0, 100, 100] = 0
    data[0, 1, 100, 100] = 4500
    data[0, 2, 100, 100] = 40000
    data[0, 3, 100, 100] = 44500
    out_gdq, row_below_gdq, row_above_gdq = find_crs(data, gdq, read_noise, rej_threshold,
                                                     rej_threshold, rej_threshold, nframes, False, 200, 10)
    assert np.array_equal([0, 0, 4, 0], out_gdq[0, :, 100, 100])


def test_3grps_1cr(setup_cube):
    """"
        A test of a three group integration that has a jump in the third group of pixel 1, 1
        """
    ngroups = 3
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups, readnoise=5 * np.sqrt(2), nrows=2, ncols=2)
    nframes = 1
    data[0, 0, 1, 1] = 0
    data[0, 1, 1, 1] = 40000
    data[0, 2, 1, 1] = 44500
    out_gdq, row_below_gdq, row_above_gdq = find_crs(data, gdq, read_noise, rej_threshold,
                                                     rej_threshold, rej_threshold, nframes, False, 200, 10)
    assert np.array_equal([0, 4, 0], out_gdq[0, :, 1, 1])


def test_4grps_2crs(setup_cube):
    """"
        A test of a four group integration that has a jump in the 2nd and fourth groups of pixel 1, 1
        """
    ngroups = 4
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups, readnoise=5 * np.sqrt(2), nrows=2, ncols=2)
    nframes = 1
    data[0, 0, 1, 1] = 0
    data[0, 1, 1, 1] = 1000
    data[0, 2, 1, 1] = 1100
    data[0, 3, 1, 1] = 10000
    out_gdq, row_below_gdq, row_above_gdq = find_crs(data, gdq, read_noise, rej_threshold,
                                                     rej_threshold, rej_threshold, nframes, False, 200, 10)
    assert np.array_equal([0, 4, 0, 4], out_gdq[0, :, 1, 1])


def test_5grps_2crs(setup_cube):
    ngroups = 5
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups, readnoise=5 * np.sqrt(2), nrows=2, ncols=2)
    nframes = 1
    data[0, 0, 1, 1] = 0
    data[0, 1, 1, 1] = 1000
    data[0, 2, 1, 1] = 1100
    data[0, 3, 1, 1] = 10000
    data[0, 4, 1, 1] = 10102
    out_gdq, row_below_gdq, row_above_gdq = find_crs(data, gdq, read_noise, rej_threshold,
                                                     rej_threshold, rej_threshold, nframes, False, 200, 10)
    assert np.array_equal([0, 4, 0, 4, 0], out_gdq[0, :, 1, 1])


def test_6grps_2crs(setup_cube):
    ngroups = 6
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups, readnoise=5 * np.sqrt(2), nrows=2, ncols=2)
    nframes = 1
    data[0, 0, 1, 1] = 0
    data[0, 1, 1, 1] = 1000
    data[0, 2, 1, 1] = 1100
    data[0, 3, 1, 1] = 10000
    data[0, 4, 1, 1] = 10102
    data[0, 5, 1, 1] = 10208
    out_gdq, row_below_gdq, row_above_gdq = find_crs(data, gdq, read_noise, rej_threshold,
                                                     rej_threshold, rej_threshold, nframes, False, 200, 10)
    assert np.array_equal([0, 4, 0, 4, 0, 0], out_gdq[0, :, 1, 1])


def test_6grps_different_valid_grps_each_pixel(setup_cube):
    ngroups = 6
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups, readnoise=5 * np.sqrt(2), nrows=3, ncols=2)
    nframes = 1
    # pixel 0, 0 no crs
    data[0, 0, 0, 0] = 0
    data[0, 1, 0, 0] = 10
    data[0, 2, 0, 0] = 20
    data[0, 3, 0, 0] = 25
    data[0, 4, 0, 0] = 35
    data[0, 5, 0, 0] = 42
    # pixel 0, 1 one crs grp 3
    data[0, 0, 0, 1] = 0
    data[0, 1, 0, 1] = 10
    data[0, 2, 0, 1] = 20
    data[0, 3, 0, 1] = 2500
    data[0, 4, 0, 1] = 2510
    data[0, 5, 0, 1] = 2522
    # pixel 1, 0 one cr, two last grps saturated
    data[0, 0, 1, 0] = 0
    data[0, 1, 1, 0] = 10
    data[0, 2, 1, 0] = 20
    data[0, 3, 1, 0] = 2500
    data[0, 4, 1, 0] = 2510
    data[0, 5, 1, 0] = 1522
    gdq[0, 4, 1, 0] = dqflags.group['SATURATED']
    gdq[0, 5, 1, 0] = dqflags.group['SATURATED']
    # pixel 1, 1 two crs
    data[0, 0, 1, 1] = 0
    data[0, 1, 1, 1] = 1000
    data[0, 2, 1, 1] = 1100
    data[0, 3, 1, 1] = 10000
    data[0, 4, 1, 1] = 10102
    data[0, 5, 1, 1] = 10208
    # pixel 2, 0 one cr, three last groups saturated
    data[0, 0, 2, 0] = 0
    data[0, 1, 2, 0] = 1000
    data[0, 2, 2, 0] = 10100
    data[0, 3, 2, 0] = 10000
    data[0, 4, 2, 0] = 10102
    data[0, 5, 2, 0] = 10208
    gdq[0, 3, 2, 0] = dqflags.group['SATURATED']
    gdq[0, 4, 2, 0] = dqflags.group['SATURATED']
    gdq[0, 5, 2, 0] = dqflags.group['SATURATED']
    # pixel 2, 1 two crs, last group saturated
    data[0, 0, 2, 1] = 0
    data[0, 1, 2, 1] = 10000
    data[0, 2, 2, 1] = 20100
    data[0, 3, 2, 1] = 20200
    data[0, 4, 2, 1] = 20302
    data[0, 5, 2, 1] = 10208
    gdq[0, 5, 2, 1] = dqflags.group['SATURATED']
    out_gdq, row_below_gdq, row_above_gdq = find_crs(data, gdq, read_noise, rej_threshold,
                                                     rej_threshold, rej_threshold, nframes, False, 200, 10)
    assert np.array_equal([0, 0, 0, 0, 0, 0], out_gdq[0, :, 0, 0])
    assert np.array_equal([0, 0, 0, 4, 0, 0], out_gdq[0, :, 0, 1])
    assert np.array_equal([0, 0, 0, 4, 2, 2], out_gdq[0, :, 1, 0])
    assert np.array_equal([0, 4, 0, 4, 0, 0], out_gdq[0, :, 1, 1])
    assert np.array_equal([0, 0, 4, 2, 2, 2], out_gdq[0, :, 2, 0])
    assert np.array_equal([0, 4, 4, 0, 0, 2], out_gdq[0, :, 2, 1])


def test_3diff_median_vector():
    indiffs = np.asarray([7, 5, 10])
    indices = np.asarray([1, 0, 2])
    median = get_clipped_median_vector(num_differences=3, diffs_to_ignore=0, input_vector=indiffs, sorted_index=indices)
    assert median == 7


def test_3diff_1sat_median_vector():
    indiffs = np.asarray([7, 5, 10])
    indices = np.asarray([1, 0, 2])
    median = get_clipped_median_vector(num_differences=3, diffs_to_ignore=1, input_vector=indiffs, sorted_index=indices)
    assert median == 5


def test_4diff_median_vector():
    indiffs = np.asarray([7, 5, 10, 12])
    indices = np.asarray([1, 0, 2, 13])
    median = get_clipped_median_vector(num_differences=4, diffs_to_ignore=0, input_vector=indiffs, sorted_index=indices)
    assert median == 7


def test_5diff_median_vector():
    indiffs = np.asarray([7, 5, 10, 12, 13])
    indices = np.asarray([1, 0, 2, 3, 4])
    median = get_clipped_median_vector(num_differences=5, diffs_to_ignore=0, input_vector=indiffs, sorted_index=indices)
    assert median == 8.5


def test_5diff_1sat_median_vector():
    indiffs = np.asarray([7, 5, 10, 12, 100000])
    indices = np.asarray([1, 0, 2, 3, 4])
    median = get_clipped_median_vector(num_differences=5, diffs_to_ignore=1, input_vector=indiffs, sorted_index=indices)
    assert median == 7


def test_5diff_2sat_median_vector():
    indiffs = np.asarray([7, 5, 10, 10000, 100000])
    indices = np.asarray([1, 0, 2, 3, 4])
    median = get_clipped_median_vector(num_differences=5, diffs_to_ignore=2, input_vector=indiffs, sorted_index=indices)
    assert median == 7


def test_5diff_3sat_median_vector():
    indiffs = np.asarray([7, 5, 10, 10000, 100000])
    indices = np.asarray([1, 0, 2, 3, 4])
    median = get_clipped_median_vector(num_differences=5, diffs_to_ignore=3, input_vector=indiffs, sorted_index=indices)
    assert median == 5


def test_4diff_median_array():
    diffs_to_skip = np.zeros(shape=(2, 2), dtype=np.int32)
    indiffs = np.zeros(shape=(2, 2, 4),dtype=np.float32)
    indices = np.zeros(shape=(2, 2, 4),dtype=np.int32)
    indiffs[0, 0] = [11, 6, 9, 13]
    indiffs[0, 1] = [9, 4, 7, 3]
    indiffs[1, 0] = [10, 5, 3, 12]
    indiffs[1, 1] = [1, 3, 8, 12]
    indices[0, 0] = [1, 2, 0, 3]
    indices[0, 1] = [3, 1, 2, 0]
    indices[1, 0] = [2, 1, 0, 3]
    indices[1, 1] = [0, 1, 2, 3]
    med_array = get_clipped_median_array(num_differences=4, diffs_to_ignore=diffs_to_skip, input_array=indiffs,
                                         sorted_index=indices)
    assert med_array[0, 0] == 9
    assert med_array[0, 1] == 4
    assert med_array[1, 0] == 5
    assert med_array[1, 1] == 3


def test_4diff_median_2pixsat_array():
    diffs_to_skip = np.zeros(shape=(2, 2), dtype=np.int32)
    indiffs = np.zeros(shape=(2, 2, 4),dtype=np.float32)
    indices = np.zeros(shape=(2, 2, 4),dtype=np.int32)
    indiffs[0, 0] = [11, 6, 9, 13]
    indiffs[0, 1] = [9, 4, 7, 3]
    indiffs[1, 0] = [10, 5, 3, 12]
    indiffs[1, 1] = [1, 3, 8, 12]
    indices[0, 0] = [1, 2, 0, 3]
    indices[0, 1] = [3, 1, 2, 0]
    indices[1, 0] = [2, 1, 0, 3]
    indices[1, 1] = [0, 1, 2, 3]
    diffs_to_skip[1, 0] = 1
    diffs_to_skip[1, 1] = 1
    med_array = get_clipped_median_array(num_differences=4, diffs_to_ignore=diffs_to_skip, input_array=indiffs,
                                         sorted_index=indices)
    assert med_array[0, 0] == 9
    assert med_array[0, 1] == 4
    assert med_array[1, 0] == 5
    assert med_array[1, 1] == 3


def test_4diff_median_2pixsat_2pixverysat_array():
    diffs_to_skip = np.zeros(shape=(2, 2), dtype=np.int32)
    indiffs = np.zeros(shape=(2, 2, 4),dtype=np.float32)
    indices = np.zeros(shape=(2, 2, 4),dtype=np.int32)
    indiffs[0, 0] = [11, 6, 9, 13]
    indiffs[0, 1] = [9, 4, 7, 3]
    indiffs[1, 0] = [10, 5, 3, 12]
    indiffs[1, 1] = [1, 3, 8, 12]
    indices[0, 0] = [1, 2, 0, 3]
    indices[0, 1] = [3, 1, 2, 0]
    indices[1, 0] = [2, 1, 0, 3]
    indices[1, 1] = [0, 1, 2, 3]
    diffs_to_skip[1, 0] = 1
    diffs_to_skip[1, 1] = 1
    diffs_to_skip[0, 0] = 2
    diffs_to_skip[0, 1] = 2
    med_array = get_clipped_median_array(num_differences=4, diffs_to_ignore=diffs_to_skip, input_array=indiffs,
                                         sorted_index=indices)
    assert med_array[0, 0] == 6  # min value
    assert med_array[0, 1] == 3  # min value
    assert med_array[1, 0] == 5  # median
    assert med_array[1, 1] == 3  # median


@pytest.fixture(scope='function')
def setup_cube():

    def _cube(ngroups, readnoise=10, nrows=204, ncols=204):
        nints = 1
        rej_threshold = 3
        nframes = 1
        data = np.zeros(shape=(nints, ngroups, nrows, ncols), dtype=np.float32)
        read_noise = np.full((nrows, ncols), readnoise, dtype=np.float32)
        gdq = np.zeros(shape=(nints, ngroups, nrows, ncols), dtype=np.uint32)

        return data, gdq, nframes, read_noise, rej_threshold

    return _cube
