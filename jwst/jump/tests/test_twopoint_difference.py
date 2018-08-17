from jwst.jump.twopoint_difference import find_CRs
import numpy as np
import pytest
from jwst.datamodels import dqflags

#@pytest.mark.skip(reason="testing skipping")
def test_noCRs_NoFlux():
    ngroups=5
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups)
    assert(0 == np.max(find_CRs(data,gdq,read_noise,rej_threshold,nframes))) # no CR found

#@pytest.mark.skip
def test_5grps_cr3_NoFlux():
    ngroups=5
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups)

    data[0, 0:2, 100, 100] = 10.0
    data[0, 2:5, 100, 100] = 1000
    medianDiff = find_CRs(data, gdq, read_noise, rej_threshold, nframes)
    assert(4 == np.max(gdq)) #a CR was found
    assert(2 == np.argmax(gdq[0,:,100,100])) #find the CR in the expected group

#@pytest.mark.skip("testing skipping")
def test_5grps_cr2_NoFlux():
    ngroups=5
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups)

    data[0, 0, 100, 100] = 10.0
    data[0, 1:6, 100, 100] = 1000
    medianDiff = find_CRs(data, gdq, read_noise, rej_threshold, nframes)
    assert(4 == np.max(gdq)) #a CR was found
    assert(1 == np.argmax(gdq[0,:,100,100])) #find the CR in the expected group

##@pytest.mark.skip("testing skipping")
def test_6grps_negative_differences_zeromedian():
    ngroups=6
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups)

    data[0, 0, 100, 100] = 100
    data[0, 1, 100, 100] = 90
    data[0, 2, 100, 100] = 95
    data[0, 3, 100, 100] = 105
    data[0, 4, 100, 100] = 100
    data[0, 5, 100, 100] = 100
    medianDiff = find_CRs(data, gdq, read_noise, rej_threshold, nframes)
    print('shape '+str(medianDiff.shape))
    print('median differences of 100,100 '+str(medianDiff[0,100,100]))
    assert(0 == np.max(gdq)) #no CR was found
    assert(0 == medianDiff[0,100,100]) #Median difference is zero


#@pytest.mark.skip("testing skipping")
def test_5grps_cr2_NegJumpFlux():
    ngroups=5
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups)

    data[0, 0, 100, 100] = 1000.0
    data[0, 1:6, 100, 100] = 10
    medianDiff = find_CRs(data, gdq, read_noise, rej_threshold, nframes)
    assert(4 == np.max(gdq)) #a CR was found
    assert(1 == np.argmax(gdq[0,:,100,100])) #find the CR in the expected group

#@pytest.mark.skip("testing skipping")
def test_3grps_cr2_NoFlux():
    ngroups=3
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups)
    data[0, 0, 100, 100] = 10.0
    data[0, 1:4, 100, 100] = 1000
    print("test data "+repr(data[0,:,100,100]))
    medianDiff = find_CRs(data, gdq, read_noise, rej_threshold, nframes)
    print(repr(gdq[0, :, 100, 100]))
    print("calculated median diff of pixel is ", medianDiff[0, 100, 100])
    assert(4 == np.max(gdq)) #a CR was found
    #    assert(1,np.argmax(gdq[0,:,100,100])) #find the CR in the expected group
    assert(np.array_equal([0, 4, 0], gdq[0, :, 100, 100]))


##@pytest.mark.xfail("testing skipping")
def test_4grps_cr2_NoFlux():
    ngroups=4
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups)
    data[0, 0, 100, 100] = 10.0
    data[0, 1:4, 100, 100] = 1000
    medianDiff = find_CRs(data, gdq, read_noise, rej_threshold, nframes)
    assert(4 == np.max(gdq)) #a CR was found
    assert(1 == np.argmax(gdq[0,:,100,100])) #find the CR in the expected group

#@pytest.mark.skip("testing skipping")
def test_5grps_cr2_nframe2():
    ngroups=5
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups)
    nframes=2
    data[0, 0, 100, 100] = 10.0
    data[0, 1, 100, 100] = 500
    data[0, 2, 100, 100] = 1002
    data[0, 3, 100, 100] = 1001
    data[0, 4, 100, 100] = 1005
    medianDiff = find_CRs(data, gdq, read_noise, rej_threshold, nframes)
    print("calculated median diff of pixel is ", medianDiff[0, 100, 100])
    print(repr(gdq[0, :, 100, 100]))
    assert(4 == np.max(gdq)) #a CR was found
    assert(np.array_equal([0,4,4,0,0], gdq[0, :, 100, 100]) )

@pytest.mark.xfail
def test_4grps_twocrs_2nd_4th():
    ngroups=4
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups)
    nframes=1
    data[0, 0, 100, 100] = 10.0
    data[0, 1, 100, 100] = 60
    data[0, 2, 100, 100] = 60
    data[0, 3, 100, 100] = 115
    medianDiff = find_CRs(data, gdq, read_noise, rej_threshold, nframes)
    print("calculated median diff of pixel is ", medianDiff[0, 100, 100])
    print(repr(gdq[0, :, 100, 100]))
    assert(4 == np.max(gdq)) #a CR was found
    assert(np.array_equal([0,4,0,4] , gdq[0, :, 100, 100]) )

#@pytest.mark.skip("testing skipping")
def test_5grps_twocrs_2nd_5th():
    ngroups=5
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups)
    nframes=1
    data[0, 0, 100, 100] = 10.0
    data[0, 1, 100, 100] = 60
    data[0, 2, 100, 100] = 60
    data[0, 3, 100, 100] = 60
    data[0, 4, 100, 100] = 115
    medianDiff = find_CRs(data, gdq, read_noise, rej_threshold, nframes)
    print("calculated median diff of pixel is ", medianDiff[0, 100, 100])
    assert(4 == np.max(gdq)) #a CR was found
    assert(np.array_equal([0,4,0,0,4] ,gdq[0, :, 100, 100]) )

#@pytest.mark.skip("testing skipping")
def test_5grps_twocrs_2nd_5thbig():
    ngroups=5
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups)
    nframes=1
    data[0, 0, 100, 100] = 10.0
    data[0, 1, 100, 100] = 60
    data[0, 2, 100, 100] = 60
    data[0, 3, 100, 100] = 60
    data[0, 4, 100, 100] = 2115
    medianDiff = find_CRs(data, gdq, read_noise, rej_threshold, nframes)
    print("calculated median diff of pixel is ", medianDiff[0, 100, 100])
    assert(4 == np.max(gdq)) #a CR was found
    assert(np.array_equal([0,4,0,0,4] , gdq[0, :, 100, 100]) )

#@pytest.mark.skip("testing skipping")
def test_10grps_twocrs_2nd_8th_big():
    ngroups=10
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
    medianDiff = find_CRs(data, gdq, read_noise, rej_threshold, nframes)
    print("calculated median diff of pixel is ", medianDiff[0, 100, 100])
    assert(4 == np.max(gdq)) #a CR was found
    assert(np.array_equal([0,4,0,0,0,0,0,4,0,0] , gdq[0, :, 100, 100]) )

#@pytest.mark.skip("testing skipping")
def test_10grps_twocrs_10percenthit():
    ngroups=10
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
    medianDiff = find_CRs(data, gdq, read_noise, rej_threshold, nframes)
    print("calculated median diff of pixel is ", medianDiff[0, 100, 100])
    assert(4 == np.max(gdq)) #a CR was found
    assert(np.array_equal([0,4,0,0,0,0,0,4,0,0] , gdq[0, :, 100, 100]) )

#@pytest.mark.skip("testing skipping")
def test_5grps_twocrs_2nd_5thbig_nframes2():
    ngroups=5
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups,readnoise=10*np.sqrt(2))
    nframes=2
    data[0, 0, 100, 100] = 10.0
    data[0, 1, 100, 100] = 60
    data[0, 2, 100, 100] = 60
    data[0, 3, 100, 100] = 60
    data[0, 4, 100, 100] = 2115
    medianDiff = find_CRs(data, gdq, read_noise, rej_threshold, nframes)
    print("calculated median diff of pixel is ", medianDiff[0, 100, 100])
    assert(4 == np.max(gdq)) #a CR was found
    assert(np.array_equal([0,4,0,0,4] , gdq[0, :, 100, 100]) )


#@pytest.mark.skip("testing skipping")
def test_6grps_twocrs_2nd_5th():
    ngroups=6
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups)
    nframes=1
    data[0, 0, 100, 100] = 10.0
    data[0, 1, 100, 100] = 60
    data[0, 2, 100, 100] = 60
    data[0, 3, 100, 100] = 60
    data[0, 4, 100, 100] = 115
    data[0, 5, 100, 100] = 115
    medianDiff = find_CRs(data, gdq, read_noise, rej_threshold, nframes)
    print("calculated median diff of pixel is ",medianDiff[0,100,100])
    assert(4 == np.max(gdq)) #a CR was found
    assert(np.array_equal([0,4,0,0,4,0] , gdq[0, :, 100, 100]) )

#@pytest.mark.skip("testing skipping")
def test_6grps_twocrs_2nd_5th_nframes2():
    ngroups=6
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups,readnoise=10*np.sqrt(2))
    nframes=2
    data[0, 0, 100, 100] = 10.0
    data[0, 1, 100, 100] = 60
    data[0, 2, 100, 100] = 60
    data[0, 3, 100, 100] = 60
    data[0, 4, 100, 100] = 115
    data[0, 5, 100, 100] = 115
    medianDiff = find_CRs(data, gdq, read_noise, rej_threshold, nframes)
    assert(4 == np.max(gdq)) #a CR was found
    assert(np.array_equal([0,4,0,0,4,0] , gdq[0, :, 100, 100]) )

#@pytest.mark.skip("testing skipping")
def test_6grps_twocrs_twopixels_nframes2():
    ngroups=6
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
    medianDiff = find_CRs(data, gdq, read_noise, rej_threshold, nframes)
    assert(4 == np.max(gdq)) #a CR was found
    print("100 100 dq",repr(gdq[0,:,100,100]))
    print("200 100 dq",repr(gdq[0,:,200,100]))
    assert(np.array_equal([0,4,0,0,4,0] , gdq[0, :, 100, 100]) )
    assert(np.array_equal([0, 0, 4, 0, 4, 0] , gdq[0, :, 200, 100]))

#@pytest.mark.skip("testing skipping")
def test_5grps_cr2_negslope():
    ngroups = 5
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups)
    nframes = 1
    data[0, 0, 100, 100] = 100.0
    data[0, 1, 100, 100] = 0
    data[0, 2, 100, 100] = -200
    data[0, 3, 100, 100] = -260
    data[0, 4, 100, 100] = -360
    medianDiff = find_CRs(data, gdq, read_noise, rej_threshold, nframes)
    print("calculated median diff of pixel is ", medianDiff[0, 100, 100])
    assert(4 == np.max(gdq))  # a CR was found
    assert(np.array_equal([0, 0, 4, 0, 0] , gdq[0, :, 100, 100]))
#@pytest.mark.skip
def test_6grps_1CR():
    ngroups = 6
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups,readnoise=10)
    nframes = 1
    data[0, 0, 100, 100] = 0
    data[0, 1, 100, 100] = 10
    data[0, 2, 100, 100] = 21
    data[0, 3, 100, 100] = 33
    data[0, 4, 100, 100] = 46
    data[0, 5, 100, 100] = 1146
    medianDiff = find_CRs(data, gdq, read_noise, rej_threshold, nframes)
    print("calculated median diff of pixel is ", medianDiff[0, 100, 100])
    assert(11 == medianDiff[0, 100, 100])
#@pytest.mark.skip
def test_7grps_1CR():
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
    medianDiff = find_CRs(data, gdq, read_noise, rej_threshold, nframes)
    print("calculated median diff of pixel is ", medianDiff[0, 100, 100])
    assert(11.5 == medianDiff[0, 100, 100])
#@pytest.mark.skip
def test_5grps_noCR():
    ngroups = 6
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups,readnoise=10)
    nframes = 1
    data[0, 0, 100, 100] = 0
    data[0, 1, 100, 100] = 10
    data[0, 2, 100, 100] = 21
    data[0, 3, 100, 100] = 33
    data[0, 4, 100, 100] = 46
    medianDiff = find_CRs(data, gdq, read_noise, rej_threshold, nframes)
    print("calculated median diff of pixel is ", medianDiff[0, 100, 100])
    assert(11 == medianDiff[0, 100, 100])
#@pytest.mark.skip
def test_6grps_noCR():
    ngroups = 6
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups,readnoise=10)
    nframes = 1
    data[0, 0, 100, 100] = 0
    data[0, 1, 100, 100] = 10
    data[0, 2, 100, 100] = 21
    data[0, 3, 100, 100] = 33
    data[0, 4, 100, 100] = 46
    data[0, 5, 100, 100] = 60
    medianDiff = find_CRs(data, gdq, read_noise, rej_threshold, nframes)
    print("calculated median diff of pixel is ", medianDiff[0, 100, 100])
    assert(11.5 == medianDiff[0, 100, 100])


#@pytest.mark.skip("testing skipping")
def test_10grps_cr2_gt3sigma():
    ngroups = 10
    crmag=16
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups,readnoise=5)
    nframes = 1
    data[0, 0, 100, 100] = 0
    data[0, 1:11, 100, 100] = crmag
    medianDiff = find_CRs(data, gdq, read_noise, rej_threshold, nframes)
    assert(4 == np.max(gdq))  # a CR was found
    assert(np.array_equal([0, 4, 0, 0, 0,0,0,0,0,0] , gdq[0, :, 100, 100]))

#@pytest.mark.skip(reason="testing skipping")
def test_10grps_cr2_3sigma_noCR():
    ngroups = 10
    crmag=15
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups,readnoise=5)
    nframes = 1
    data[0, 0, 100, 100] = 0
    data[0, 1:11, 100, 100] = crmag
    medianDiff = find_CRs(data, gdq, read_noise, rej_threshold, nframes)
    assert(0 == np.max(gdq))  # a CR was found
    assert(np.array_equal([0, 0, 0, 0, 0,0,0,0,0,0] , gdq[0, :, 100, 100]))

#@pytest.mark.skip(reason="testing skipping")
def test_10grps_cr2_gt3sigma_2frames():
    ngroups = 10
    crmag=16
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups,readnoise=5*np.sqrt(2))
    nframes = 2
    data[0, 0, 100, 100] = 0
    data[0, 1:11, 100, 100] = crmag
    medianDiff = find_CRs(data, gdq, read_noise, rej_threshold, nframes)
    assert(4 == np.max(gdq))  # a CR was found
    assert(np.array_equal([0, 4, 0, 0, 0,0,0,0,0,0] , gdq[0, :, 100, 100]))

#@pytest.mark.skip("testing skipping")
def test_10grps_cr2_3sigma_2frames_noCR():
    ngroups = 10
    crmag=15
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups,readnoise=5*np.sqrt(2))
    nframes = 2
    data[0, 0, 100, 100] = 0
    data[0, 1:11, 100, 100] = crmag
    medianDiff = find_CRs(data, gdq, read_noise, rej_threshold, nframes)
    assert(0 == np.max(gdq))  # a CR was found
    assert(np.array_equal([0, 0, 0, 0, 0, 0, 0, 0, 0, 0] , gdq[0, :, 100, 100]))

#@pytest.mark.skip(reason="testing skipping")
def test_10grps_noCR_2pixels_sigma0():
    ngroups = 10
    crmag = 15
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups, readnoise=5 * np.sqrt(2))
    nframes=1
    data[0, 0, 100, 100] = crmag
    data[0, 1:11, 100, 100] = crmag
    read_noise[500, 500] = 0.0
    read_noise[600, 600] = 0.0
    medianDiff = find_CRs(data, gdq, read_noise, rej_threshold, nframes)
    assert(0 == np.max(gdq))  # no CR was found
#@pytest.mark.skip
def test_5grps_Satat4_CRat3():
    ngroups = 5
    #crmag = 1000
    data, gdq, nframes, read_noise, rej_threshold = setup_cube(ngroups, readnoise=5 * np.sqrt(2))
    nframes=1
    data[0, 0, 100, 100] = 10000
    data[0, 1, 100, 100] = 30000
    data[0, 2, 100, 100] = 60000
    data[0, 3, 100, 100] = 61000
    data[0, 4, 100, 100] = 61000
    gdq[0, 3, 100, 100] = dqflags.group['SATURATED']
    gdq[0, 4, 100, 100] = dqflags.group['SATURATED']
    medianDiff = find_CRs(data, gdq, read_noise, rej_threshold, nframes)
   # assert(4 == np.max(gdq))  # no CR was found
    assert (np.array_equal([0, 0, dqflags.group['JUMP_DET'],dqflags.group['SATURATED'], dqflags.group['SATURATED']], gdq[0, :, 100, 100]))

#@pytest.mark.skip
def test_6grps_Satat6_CRat1():
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
    medianDiff = find_CRs(data, gdq, read_noise, rej_threshold, nframes)
   # assert(4 == np.max(gdq))  # no CR was found
    assert (np.array_equal([0, dqflags.group['JUMP_DET'], 0,0,0, dqflags.group['SATURATED']], gdq[0, :, 100, 100]))


@pytest.mark.xfail
def test_6grps_Satat6_CRat1_flagadjpixels():
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
    medianDiff = find_CRs(data, gdq, read_noise, rej_threshold, nframes)
   # assert(4 == np.max(gdq))  # no CR was found
    assert (np.array_equal([0, dqflags.group['JUMP_DET'], 0,0,0, dqflags.group['SATURATED']], gdq[0, :, 100, 100]))
    assert (np.array_equal([0, dqflags.group['JUMP_DET'], 0, 0, 0, dqflags.group['SATURATED']], gdq[0, :, 99, 100]))
    assert (np.array_equal([0, dqflags.group['JUMP_DET'], 0, 0, 0, dqflags.group['SATURATED']], gdq[0, :, 101, 100]))
    assert (np.array_equal([0, dqflags.group['JUMP_DET'], 0, 0, 0, dqflags.group['SATURATED']], gdq[0, :, 100, 99]))
    assert (np.array_equal([0, dqflags.group['JUMP_DET'], 0, 0, 0, dqflags.group['SATURATED']], gdq[0, :, 100, 101]))
    
#@pytest.mark.skip
def test_10grps_Satat8_CRsat3and6():
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
    medianDiff = find_CRs(data, gdq, read_noise, rej_threshold, nframes)
   # assert(4 == np.max(gdq))  # no CR was found
    assert (np.array_equal([0, 0, dqflags.group['JUMP_DET'], 0, 0, dqflags.group['JUMP_DET'],
                            0,dqflags.group['SATURATED'],dqflags.group['SATURATED'],dqflags.group['SATURATED']], gdq[0, :, 100, 100]))


def setup_cube(ngroups,readnoise=10):
    nints = 1
    nrows = 2048
    ncols = 2048
    rej_threshold = 3
    nframes = 1
    data = np.zeros(shape=(nints, ngroups, nrows, ncols), dtype=np.float32)
    read_noise = np.zeros(shape=(nrows, ncols), dtype=np.float32)
    read_noise[:, :] = readnoise
    #primary_gdq = np.zeros(shape=(nints, ngroups, nrows, ncols), dtype=np.int32)
    gdq = np.zeros(shape=(nints, ngroups, nrows, ncols), dtype=np.int32)
    return data, gdq, nframes, read_noise, rej_threshold

