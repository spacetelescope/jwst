import pytest
import numpy as np

from jwst.jump.twopoint_difference import find_crs
from jwst.datamodels import dqflags


def test_nocrs_noflux(setup_cube):
    ngroups = 5
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups)
    median_diff, out_gdq = find_crs(data, gdq, read_noise, rej_threshold, nframes)
    assert(0 == np.max(out_gdq)) # no CR found


def test_5grps_cr3_noflux(setup_cube):
    ngroups = 5
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups)

    data[0, 0:2, 100, 100] = 10.0
    data[0, 2:5, 100, 100] = 1000
    median_diff, out_gdq =find_crs(data, gdq, read_noise, rej_threshold, nframes)
    assert(4 == np.max(out_gdq)) #a CR was found
    assert(2 == np.argmax(out_gdq[0,:,100,100])) #find the CR in the expected group


def test_5grps_cr2_noflux(setup_cube):
    ngroups = 5
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups)

    data[0, 0, 100, 100] = 10.0
    data[0, 1:6, 100, 100] = 1000
    median_diff, out_gdq =find_crs(data, gdq, read_noise, rej_threshold, nframes)
    assert(4 == np.max(out_gdq)) #a CR was found
    assert(1 == np.argmax(out_gdq[0,:,100,100])) #find the CR in the expected group


def test_6grps_negative_differences_zeromedian(setup_cube):
    ngroups = 6
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups)

    data[0, 0, 100, 100] = 100
    data[0, 1, 100, 100] = 90
    data[0, 2, 100, 100] = 95
    data[0, 3, 100, 100] = 105
    data[0, 4, 100, 100] = 100
    data[0, 5, 100, 100] = 100
    median_diff, out_gdq = find_crs(data, gdq, read_noise, rej_threshold, nframes)
    assert(0 == np.max(out_gdq)) #no CR was found
    assert(0 == median_diff[0,100,100]) #Median difference is zero


def test_5grps_cr2_negjumpflux(setup_cube):
    ngroups = 5
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups)

    data[0, 0, 100, 100] = 1000.0
    data[0, 1:6, 100, 100] = 10
    median_diff, out_gdq =find_crs(data, gdq, read_noise, rej_threshold, nframes)
    assert(4 == np.max(out_gdq)) #a CR was found
    assert(1 == np.argmax(out_gdq[0,:,100,100])) #find the CR in the expected group


def test_3grps_cr2_noflux(setup_cube):
    ngroups = 3
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups)
    data[0, 0, 100, 100] = 10.0
    data[0, 1:4, 100, 100] = 1000
    median_diff, out_gdq = find_crs(data, gdq, read_noise, rej_threshold, nframes)
    assert(4 == np.max(out_gdq)) # a CR was found
    #    assert(1,np.argmax(out_gdq[0,:,100,100])) #find the CR in the expected group
    assert(np.array_equal([0, 4, 0], out_gdq[0, :, 100, 100]))


def test_4grps_cr2_noflux(setup_cube):
    ngroups = 4
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups)
    data[0, 0, 100, 100] = 10.0
    data[0, 1:4, 100, 100] = 1000
    median_diff, out_gdq = find_crs(data, gdq, read_noise, rej_threshold, nframes)
    assert(4 == np.max(out_gdq)) #a CR was found
    assert(1 == np.argmax(out_gdq[0,:,100,100])) #find the CR in the expected group


def test_5grps_cr2_nframe2(setup_cube):
    ngroups = 5
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups)
    nframes=2
    data[0, 0, 100, 100] = 10.0
    data[0, 1, 100, 100] = 500
    data[0, 2, 100, 100] = 1002
    data[0, 3, 100, 100] = 1001
    data[0, 4, 100, 100] = 1005
    median_diff, out_gdq = find_crs(data, gdq, read_noise, rej_threshold, nframes)
    assert(4 == np.max(out_gdq)) #a CR was found
    assert(np.array_equal([0,4,4,0,0], out_gdq[0, :, 100, 100]) )


@pytest.mark.xfail
def test_4grps_twocrs_2nd_4th(setup_cube):
    ngroups = 4
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups)
    nframes=1
    data[0, 0, 100, 100] = 10.0
    data[0, 1, 100, 100] = 60
    data[0, 2, 100, 100] = 60
    data[0, 3, 100, 100] = 115
    median_diff, out_gdq = find_crs(data, gdq, read_noise, rej_threshold, nframes)
    assert(4 == np.max(out_gdq)) #a CR was found
    assert(np.array_equal([0,4,0,4] , out_gdq[0, :, 100, 100]) )


def test_5grps_twocrs_2nd_5th(setup_cube):
    ngroups = 5
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups)
    nframes=1
    data[0, 0, 100, 100] = 10.0
    data[0, 1, 100, 100] = 60
    data[0, 2, 100, 100] = 60
    data[0, 3, 100, 100] = 60
    data[0, 4, 100, 100] = 115
    median_diff, out_gdq = find_crs(data, gdq, read_noise, rej_threshold, nframes)
    assert(4 == np.max(out_gdq)) #a CR was found
    assert(np.array_equal([0,4,0,0,4] ,out_gdq[0, :, 100, 100]) )


def test_5grps_twocrs_2nd_5thbig(setup_cube):
    ngroups = 5
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups)
    nframes=1
    data[0, 0, 100, 100] = 10.0
    data[0, 1, 100, 100] = 60
    data[0, 2, 100, 100] = 60
    data[0, 3, 100, 100] = 60
    data[0, 4, 100, 100] = 2115
    median_diff, out_gdq = find_crs(data, gdq, read_noise, rej_threshold, nframes)
    assert(4 == np.max(out_gdq)) #a CR was found
    assert(np.array_equal([0,4,0,0,4] , out_gdq[0, :, 100, 100]) )


def test_10grps_twocrs_2nd_8th_big(setup_cube):
    ngroups = 10
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups)
    nframes=1
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
    median_diff, out_gdq = find_crs(data, gdq, read_noise, rej_threshold, nframes)
    assert(4 == np.max(out_gdq)) #a CR was found
    assert(np.array_equal([0,4,0,0,0,0,0,4,0,0] , out_gdq[0, :, 100, 100]) )


def test_10grps_twocrs_10percenthit(setup_cube):
    ngroups = 10
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups)
    nframes=2
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
    median_diff, out_gdq = find_crs(data, gdq, read_noise, rej_threshold, nframes)
    assert(4 == np.max(out_gdq)) #a CR was found
    assert(np.array_equal([0,4,0,0,0,0,0,4,0,0] , out_gdq[0, :, 100, 100]) )


def test_5grps_twocrs_2nd_5thbig_nframes2(setup_cube):
    ngroups = 5
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups,readnoise=10*np.sqrt(2))
    nframes=2
    data[0, 0, 100, 100] = 10.0
    data[0, 1, 100, 100] = 60
    data[0, 2, 100, 100] = 60
    data[0, 3, 100, 100] = 60
    data[0, 4, 100, 100] = 2115
    median_diff, out_gdq = find_crs(data, gdq, read_noise, rej_threshold, nframes)
    assert(4 == np.max(out_gdq)) #a CR was found
    assert(np.array_equal([0,4,0,0,4] , out_gdq[0, :, 100, 100]) )


def test_6grps_twocrs_2nd_5th(setup_cube):
    ngroups = 6
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups)
    nframes=1
    data[0, 0, 100, 100] = 10.0
    data[0, 1, 100, 100] = 60
    data[0, 2, 100, 100] = 60
    data[0, 3, 100, 100] = 60
    data[0, 4, 100, 100] = 115
    data[0, 5, 100, 100] = 115
    median_diff, out_gdq = find_crs(data, gdq, read_noise, rej_threshold, nframes)
    assert(4 == np.max(out_gdq)) #a CR was found
    assert(np.array_equal([0,4,0,0,4,0] , out_gdq[0, :, 100, 100]) )


def test_6grps_twocrs_2nd_5th_nframes2(setup_cube):
    ngroups = 6
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups,readnoise=10*np.sqrt(2))
    nframes=2
    data[0, 0, 100, 100] = 10.0
    data[0, 1, 100, 100] = 60
    data[0, 2, 100, 100] = 60
    data[0, 3, 100, 100] = 60
    data[0, 4, 100, 100] = 115
    data[0, 5, 100, 100] = 115
    median_diff, out_gdq = find_crs(data, gdq, read_noise, rej_threshold, nframes)
    assert(4 == np.max(out_gdq)) #a CR was found
    assert(np.array_equal([0,4,0,0,4,0] , out_gdq[0, :, 100, 100]) )


def test_6grps_twocrs_twopixels_nframes2(setup_cube):
    ngroups = 6
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups,readnoise=10*np.sqrt(2))
    nframes=2
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
    median_diff, out_gdq = find_crs(data, gdq, read_noise, rej_threshold, nframes)
    assert(4 == np.max(out_gdq)) #a CR was found
    assert(np.array_equal([0,4,0,0,4,0] , out_gdq[0, :, 100, 100]) )
    assert(np.array_equal([0, 0, 4, 0, 4, 0] , out_gdq[0, :, 200, 100]))


def test_5grps_cr2_negslope(setup_cube):
    ngroups = 5
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups)
    nframes = 1
    data[0, 0, 100, 100] = 100.0
    data[0, 1, 100, 100] = 0
    data[0, 2, 100, 100] = -200
    data[0, 3, 100, 100] = -260
    data[0, 4, 100, 100] = -360
    median_diff, out_gdq = find_crs(data, gdq, read_noise, rej_threshold, nframes)
    assert(4 == np.max(out_gdq))  # a CR was found
    assert(np.array_equal([0, 0, 4, 0, 0] , out_gdq[0, :, 100, 100]))


def test_6grps_1cr(setup_cube):
    ngroups = 6
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups,readnoise=10)
    nframes = 1
    data[0, 0, 100, 100] = 0
    data[0, 1, 100, 100] = 10
    data[0, 2, 100, 100] = 21
    data[0, 3, 100, 100] = 33
    data[0, 4, 100, 100] = 46
    data[0, 5, 100, 100] = 1146
    median_diff, out_gdq = find_crs(data, gdq, read_noise, rej_threshold, nframes)
    assert (4 == out_gdq[0, 5, 100, 100])
    assert(11 == median_diff[0, 100, 100])


def test_7grps_1cr(setup_cube):
    ngroups = 7
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups,readnoise=10)
    nframes = 1
    data[0, 0, 100, 100] = 0
    data[0, 1, 100, 100] = 10
    data[0, 2, 100, 100] = 21
    data[0, 3, 100, 100] = 33
    data[0, 4, 100, 100] = 46
    data[0, 5, 100, 100] = 60
    data[0, 6, 100, 100] = 1160
    median_diff, out_gdq = find_crs(data, gdq, read_noise, rej_threshold, nframes)
    assert(4 == out_gdq[0, 6,100,100])
    assert(11.5 == median_diff[0, 100, 100])

def test_5grps_nocr(setup_cube):
    ngroups = 6
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups,readnoise=10)
    nframes = 1
    data[0, 0, 100, 100] = 0
    data[0, 1, 100, 100] = 10
    data[0, 2, 100, 100] = 21
    data[0, 3, 100, 100] = 33
    data[0, 4, 100, 100] = 46
    median_diff, out_gdq = find_crs(data, gdq, read_noise, rej_threshold, nframes)
    assert(11 == median_diff[0, 100, 100])


def test_6grps_nocr(setup_cube):
    ngroups = 6
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups,readnoise=10)
    nframes = 1
    data[0, 0, 100, 100] = 0
    data[0, 1, 100, 100] = 10
    data[0, 2, 100, 100] = 21
    data[0, 3, 100, 100] = 33
    data[0, 4, 100, 100] = 46
    data[0, 5, 100, 100] = 60
    median_diff, out_gdq = find_crs(data, gdq, read_noise, rej_threshold, nframes)
    assert(11.5 == median_diff[0, 100, 100])


def test_10grps_cr2_gt3sigma(setup_cube):
    ngroups = 10
    crmag=16
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups,readnoise=5)
    nframes = 1
    data[0, 0, 100, 100] = 0
    data[0, 1:11, 100, 100] = crmag
    median_diff, out_gdq = find_crs(data, gdq, read_noise, rej_threshold, nframes)
    assert(4 == np.max(out_gdq))  # a CR was found
    assert(np.array_equal([0, 4, 0, 0, 0,0,0,0,0,0] , out_gdq[0, :, 100, 100]))


def test_10grps_cr2_3sigma_nocr(setup_cube):
    ngroups = 10
    crmag=15
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups,readnoise=5)
    nframes = 1
    data[0, 0, 100, 100] = 0
    data[0, 1:11, 100, 100] = crmag
    median_diff, out_gdq = find_crs(data, gdq, read_noise, rej_threshold, nframes)
    assert(0 == np.max(out_gdq))  # a CR was found
    assert(np.array_equal([0, 0, 0, 0, 0,0,0,0,0,0] , out_gdq[0, :, 100, 100]))


def test_10grps_cr2_gt3sigma_2frames(setup_cube):
    ngroups = 10
    crmag=16
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups,readnoise=5*np.sqrt(2))
    nframes = 2
    data[0, 0, 100, 100] = 0
    data[0, 1:11, 100, 100] = crmag
    median_diff, out_gdq = find_crs(data, gdq, read_noise, rej_threshold, nframes)
    assert(4 == np.max(out_gdq))  # a CR was found
    assert(np.array_equal([0, 4, 0, 0, 0,0,0,0,0,0] , out_gdq[0, :, 100, 100]))

def test_10grps_cr2_gt3sigma_2frames_offdiag(setup_cube):
    ngroups = 10
    crmag=16
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups,readnoise=5*np.sqrt(2))
    nframes = 2
    data[0, 0, 100, 110] = 0
    data[0, 1:11, 100, 110] = crmag
    median_diff, out_gdq = find_crs(data, gdq, read_noise, rej_threshold, nframes)
    assert(4 == np.max(out_gdq))  # a CR was found
    assert(np.array_equal([0, 4, 0, 0, 0,0,0,0,0,0] , out_gdq[0, :, 100, 110]))

def test_10grps_cr2_3sigma_2frames_nocr(setup_cube):
    ngroups = 10
    crmag = 15
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups,readnoise=5*np.sqrt(2))
    nframes = 2
    data[0, 0, 100, 100] = 0
    data[0, 1:11, 100, 100] = crmag
    median_diff, out_gdq = find_crs(data, gdq, read_noise, rej_threshold, nframes)
    assert(0 == np.max(out_gdq))  # a CR was found
    assert(np.array_equal([0, 0, 0, 0, 0, 0, 0, 0, 0, 0] , out_gdq[0, :, 100, 100]))


def test_10grps_nocr_2pixels_sigma0(setup_cube):
    ngroups = 10
    crmag = 15
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups, readnoise=5 * np.sqrt(2))
    nframes=1
    data[0, 0, 100, 100] = crmag
    data[0, 1:11, 100, 100] = crmag
    read_noise[50, 50] = 0.0
    read_noise[60, 60] = 0.0
    median_diff, out_gdq = find_crs(data, gdq, read_noise, rej_threshold, nframes)
    assert(0 == np.max(out_gdq))  # no CR was found


def test_5grps_satat4_crat3(setup_cube):
    ngroups = 5
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups, readnoise=5 * np.sqrt(2))
    nframes=1
    data[0, 0, 100, 100] = 10000
    data[0, 1, 100, 100] = 30000
    data[0, 2, 100, 100] = 60000
    data[0, 3, 100, 100] = 61000
    data[0, 4, 100, 100] = 61000
    gdq[0, 3, 100, 100] = dqflags.group['SATURATED']
    gdq[0, 4, 100, 100] = dqflags.group['SATURATED']
    median_diff, out_gdq = find_crs(data, gdq, read_noise, rej_threshold, nframes)
    # assert(4 == np.max(out_gdq))  # no CR was found
    assert (np.array_equal([0, 0, dqflags.group['JUMP_DET'],dqflags.group['SATURATED'], dqflags.group['SATURATED']], out_gdq[0, :, 100, 100]))


def test_6grps_satat6_crat1(setup_cube):
    ngroups = 6
    #crmag = 1000
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups, readnoise=5 * np.sqrt(2))
    nframes=1
    data[0, 0, 100, 100] = 10000
    data[0, 1, 100, 100] = 35000 #CR
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
    median_diff, out_gdq = find_crs(data, gdq, read_noise, rej_threshold, nframes)
    # assert(4 == np.max(out_gdq))  # no CR was found
    assert (np.array_equal([0, dqflags.group['JUMP_DET'], 0,0,0, dqflags.group['SATURATED']], out_gdq[0, :, 100, 100]))


@pytest.mark.xfail
def test_6grps_satat6_crat1_flagadjpixels(setup_cube):
    ngroups = 6
    #crmag = 1000
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups, readnoise=5 * np.sqrt(2))
    nframes=1
    data[0, 0, 100, 100] = 10000
    data[0, 1, 100, 100] = 35000 #CR
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
    median_diff, out_gdq = find_crs(data, gdq, read_noise, rej_threshold, nframes)
   # assert(4 == np.max(out_gdq))  # no CR was found
    assert (np.array_equal([0, dqflags.group['JUMP_DET'], 0,0,0, dqflags.group['SATURATED']], out_gdq[0, :, 100, 100]))
    assert (np.array_equal([0, dqflags.group['JUMP_DET'], 0, 0, 0, dqflags.group['SATURATED']], out_gdq[0, :, 99, 100]))
    assert (np.array_equal([0, dqflags.group['JUMP_DET'], 0, 0, 0, dqflags.group['SATURATED']], out_gdq[0, :, 101, 100]))
    assert (np.array_equal([0, dqflags.group['JUMP_DET'], 0, 0, 0, dqflags.group['SATURATED']], out_gdq[0, :, 100, 99]))
    assert (np.array_equal([0, dqflags.group['JUMP_DET'], 0, 0, 0, dqflags.group['SATURATED']], out_gdq[0, :, 100, 101]))
    
    
def test_10grps_satat8_crsat3and6(setup_cube):
    ngroups = 10
    #crmag = 1000
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups, readnoise=5 * np.sqrt(2))
    nframes=1
    data[0, 0, 100, 100] = 0
    data[0, 1, 100, 100] = 5000
    data[0, 2, 100, 100] = 15000 #CR
    data[0, 3, 100, 100] = 20000
    data[0, 4, 100, 100] = 25000
    data[0, 5, 100, 100] = 40000 #CR
    data[0, 6, 100, 100] = 45000
    data[0, 7:11, 100, 100] = 61000
    gdq[0, 7:11, 100, 100] = dqflags.group['SATURATED']
    median_diff, out_gdq = find_crs(data, gdq, read_noise, rej_threshold, nframes)
   # assert(4 == np.max(out_gdq))  # no CR was found
    assert (np.array_equal([0, 0, dqflags.group['JUMP_DET'], 0, 0, dqflags.group['JUMP_DET'],
                            0,dqflags.group['SATURATED'],dqflags.group['SATURATED'],dqflags.group['SATURATED']], out_gdq[0, :, 100, 100]))


@pytest.fixture(scope='function')
def setup_cube():

    def _cube(ngroups, readnoise=10):
        nints = 1
        nrows = 204
        ncols = 204
        rej_threshold = 3
        nframes = 1
        data = np.zeros(shape=(nints, ngroups, nrows, ncols), dtype=np.float32)
        read_noise = np.full((nrows, ncols), readnoise, dtype=np.float32)
        gdq = np.zeros(shape=(nints, ngroups, nrows, ncols), dtype=np.int32)

        return data, gdq, nframes, read_noise, rej_threshold

    return _cube
