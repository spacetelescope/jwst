import pytest
import numpy as np

from jwst.jump.jump import detect_jumps
from jwst.datamodels import dqflags
from jwst.datamodels import MIRIRampModel
from jwst.datamodels import GainModel, ReadnoiseModel


def test_nocrs_noflux():
    # all pixel values are zero. So slope should be zero
    model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=5)
    out_model = detect_jumps(model1, gain, rnModel, 4.0, False, 4.0)
    assert (0 == np.max(out_model.groupdq))



def test_oneCR_10_groups():
    grouptime = 3.0
    deltaDN = 5
    ingain = 200  # use large gain to show that Poisson noise doesn't affect the recombination
    inreadnoise = np.float64(7)
    ngroups = 10
    model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=ngroups,
                                                          gain=ingain, readnoise=inreadnoise, deltatime=grouptime)
    # two segments perfect fit, second segment has twice the slope
    model1.data[0, 0, 500, 500] = 15.0
    model1.data[0, 1, 500, 500] = 20.0
    model1.data[0, 2, 500, 500] = 25.0
    model1.data[0, 3, 500, 500] = 30.0
    model1.data[0, 4, 500, 500] = 35.0
    model1.data[0, 5, 500, 500] = 140.0
    model1.data[0, 6, 500, 500] = 150.0
    model1.data[0, 7, 500, 500] = 160.0
    model1.data[0, 8, 500, 500] = 170.0
    model1.data[0, 9, 500, 500] = 180.0
    out_model = detect_jumps(model1, gain, rnModel, 4.0, False, 4.0)
    assert (4 == np.max(out_model.groupdq[0, 5, 500, 500]))

def test_oneCR_10_groups_fullarray():
    grouptime = 3.0
    deltaDN = 5
    ingain = 5  # use large gain to show that Poisson noise doesn't affect the recombination
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
    # move the CR to group 3 for row 102 and make difference be 30
    model1.data[0, 3, :, 102] = 100
    model1.data[0, 4, :, 102] = 130
    model1.data[0, 5, :, 102] = 160
    model1.data[0, 6, :, 102] = 190
    model1.data[0, 7, :, 102] = 220
    model1.data[0, 8, :, 102] = 250
    model1.data[0, 9, :, 102] = 280
    out_model = detect_jumps(model1, gain, rnModel, 4.0, False, 4.0)
    assert (4 == np.max(out_model.groupdq[0, 5, :, :]))

def test_oneCR_100_groups_fullarray():
    grouptime = 3.0
    deltaDN = 5
    ingain = 5  # use large gain to show that Poisson noise doesn't affect the recombination
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
    model1.data[0, 10:99, :, :] =190.0
    model1.data[0, 30:99, :, :] = 490.0
    out_model = detect_jumps(model1, gain, rnModel, 4.0, False, 4.0)
    assert (4 == np.max(out_model.groupdq[0, 5, :, :]))
    outdqcr = out_model.groupdq[0, 5, :, :]
 #   np.testing.assert_equal(4,out_model.groupdq[0, 5, :, :])
    badpixel = np.where(outdqcr != 4)
    np.testing.assert_allclose(4, outdqcr)



# Need test for multi-ints near zero with positive and negative slopes

def setup_inputs(ngroups=10, readnoise=10, nints=1,
                 nrows=1032, ncols=1024, nframes=1, grouptime=1.0, gain=1, deltatime=1):
    print('readnoise', readnoise)
    print('gain', gain)
    times = np.array(list(range(ngroups)), dtype=np.float64) * deltatime
    gain = np.ones(shape=(nrows, ncols), dtype=np.float64) * gain
    err = np.ones(shape=(nints, ngroups, nrows, ncols), dtype=np.float64)
    data = np.zeros(shape=(nints, ngroups, nrows, ncols), dtype=np.float64)
    pixdq = np.zeros(shape=(nrows, ncols), dtype=np.float64)
    read_noise = np.full((nrows, ncols), readnoise, dtype=np.float64)
    gdq = np.zeros(shape=(nints, ngroups, nrows, ncols), dtype=np.int32)
    model1 = MIRIRampModel(data=data, err=err, pixeldq=pixdq, groupdq=gdq, times=times)
    model1.meta.instrument.name = 'MIRI'
    model1.meta.instrument.detector = 'MIRIMAGE'
    model1.meta.instrument.filter = 'F480M'
    model1.meta.observation.date = '2015-10-13'
    model1.meta.exposure.type = 'MIR_IMAGE'
    model1.meta.exposure.group_time = deltatime
    model1.meta.subarray.name = 'FULL'
    model1.meta.subarray.xstart = 1
    model1.meta.subarray.ystart = 1
    model1.meta.subarray.xsize = 1024
    model1.meta.subarray.ysize = 1032
    model1.meta.exposure.frame_time = deltatime
    model1.meta.exposure.ngroups = ngroups
    model1.meta.exposure.group_time = deltatime
    model1.meta.exposure.nframes = 1
    model1.meta.exposure.groupgap = 0
    gain = GainModel(data=gain)
    gain.meta.instrument.name = 'MIRI'
    gain.meta.subarray.xstart = 1
    gain.meta.subarray.ystart = 1
    gain.meta.subarray.xsize = 1024
    gain.meta.subarray.ysize = 1032
    rnModel = ReadnoiseModel(data=read_noise)
    rnModel.meta.instrument.name = 'MIRI'
    rnModel.meta.subarray.xstart = 1
    rnModel.meta.subarray.ystart = 1
    rnModel.meta.subarray.xsize = 1024
    rnModel.meta.subarray.ysize = 1032
    return model1, gdq, rnModel, pixdq, err, gain


def setup_subarray_inputs(ngroups=10, readnoise=10, nints=1,
                          nrows=1032, ncols=1024, subxstart=1, subystart=1,
                          subxsize=1024, subysize=1032, nframes=1,
                          grouptime=1.0, gain=1, deltatime=1):
    times = np.array(list(range(ngroups)), dtype=np.float64) * deltatime
    gain = np.ones(shape=(nrows, ncols), dtype=np.float64) * gain
    err = np.ones(shape=(nints, ngroups, nrows, ncols), dtype=np.float64)
    data = np.zeros(shape=(nints, ngroups, subysize, subxsize), dtype=np.float64)
    pixdq = np.zeros(shape=(subysize, subxsize), dtype=np.float64)
    read_noise = np.full((nrows, ncols), readnoise, dtype=np.float64)
    gdq = np.zeros(shape=(nints, ngroups, subysize, subxsize), dtype=np.int32)
    model1 = MIRIRampModel(data=data, err=err, pixeldq=pixdq, groupdq=gdq, times=times)
    model1.meta.instrument.name = 'MIRI'
    model1.meta.instrument.detector = 'MIRIMAGE'
    model1.meta.instrument.filter = 'F480M'
    model1.meta.observation.date = '2015-10-13'
    model1.meta.exposure.type = 'MIR_IMAGE'
    model1.meta.exposure.group_time = deltatime
    model1.meta.subarray.name = 'FULL'
    model1.meta.subarray.xstart = subxstart
    model1.meta.subarray.ystart = subystart
    model1.meta.subarray.xsize = subxsize
    model1.meta.subarray.ysize = subysize
    model1.meta.exposure.frame_time = deltatime
    model1.meta.exposure.ngroups = ngroups
    model1.meta.exposure.group_time = deltatime
    model1.meta.exposure.nframes = 1
    model1.meta.exposure.groupgap = 0
    gain = GainModel(data=gain)
    gain.meta.instrument.name = 'MIRI'
    gain.meta.subarray.xstart = 1
    gain.meta.subarray.ystart = 1
    gain.meta.subarray.xsize = 1024
    gain.meta.subarray.ysize = 1032
    rnModel = ReadnoiseModel(data=read_noise)
    rnModel.meta.instrument.name = 'MIRI'
    rnModel.meta.subarray.xstart = 1
    rnModel.meta.subarray.ystart = 1
    rnModel.meta.subarray.xsize = 1024
    rnModel.meta.subarray.ysize = 1032
    return model1, gdq, rnModel, pixdq, err, gain
