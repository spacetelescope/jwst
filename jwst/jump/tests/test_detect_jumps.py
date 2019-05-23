import numpy as np
import pytest

from jwst.datamodels import GainModel, ReadnoiseModel
from jwst.datamodels import MIRIRampModel
from jwst.jump.jump import detect_jumps


def test_nocrs_noflux(setup_inputs):
    # all pixel values are zero. So slope should be zero
    model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=5)
    out_model = detect_jumps(model1, gain, rnModel, 4.0, False, 4.0)
    assert (0 == np.max(out_model.groupdq))


def test_onecr_10_groups(setup_inputs):
    grouptime = 3.0
    ingain = 200  # use large gain to show that Poisson noise doesn't affect the recombination
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
    out_model = detect_jumps(model1, gain, rnModel, 4.0, False, 4.0)
    assert (4 == np.max(out_model.groupdq[0, 5, 5, 5]))

def test_onecr_10_groups_fullarray(setup_inputs):
    grouptime = 3.0
    ingain = 5  # use large gain to show that Poisson noise doesn't affect the recombination
    inreadnoise = np.float64(7)
    ngroups = 10
    model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=ngroups,
                                                          gain=ingain, readnoise=inreadnoise, deltatime=grouptime)
    # two segments perfect fit, second segment has twice the slope
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
    # move the CR to group 3 for row 10 and make difference be 30
    model1.data[0, 3, 5, 10] = 100
    model1.data[0, 4, 5, 10] = 130
    model1.data[0, 5, 5, 10] = 160
    model1.data[0, 6, 5, 10] = 190
    model1.data[0, 7, 5, 10] = 220
    model1.data[0, 8, 5, 10] = 250
    model1.data[0, 9, 5, 10] = 280
    out_model = detect_jumps(model1, gain, rnModel, 4.0, False, 4.0)
    assert (4 == np.max(out_model.groupdq[0, 5, :, :]))


def test_onecr_50_groups(setup_inputs):
    grouptime = 3.0
    ingain = 5  # use large gain to show that Poisson noise doesn't affect the recombination
    inreadnoise = np.float64(7)
    ngroups = 50
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
    model1.data[0, 10:29, 5, 5] = 190.0
    model1.data[0, 30:49, 5, 5] = 490.0
    out_model = detect_jumps(model1, gain, rnModel, 4.0, False, 4.0)
    assert (4 == np.max(out_model.groupdq[0, 5, 5, 5]))
    outdqcr = out_model.groupdq[0, 5, 5, 5]
    np.testing.assert_allclose(4, outdqcr)

# Need test for multi-ints near zero with positive and negative slopes

@pytest.fixture
def setup_inputs():
    def _setup(ngroups=10, readnoise=10, nints=1,
                 nrows=1032, ncols=1024, nframes=1, grouptime=1.0, gain=1, deltatime=1):
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
        model1.meta.subarray.xsize = 20
        model1.meta.subarray.ysize = 20
        model1.meta.exposure.frame_time = deltatime
        model1.meta.exposure.ngroups = ngroups
        model1.meta.exposure.group_time = deltatime
        model1.meta.exposure.nframes = 1
        model1.meta.exposure.groupgap = 0
        gain = GainModel(data=gain)
        gain.meta.instrument.name = 'MIRI'
        gain.meta.subarray.xstart = 1
        gain.meta.subarray.ystart = 1
        gain.meta.subarray.xsize = 20
        gain.meta.subarray.ysize = 20
        rnModel = ReadnoiseModel(data=read_noise)
        rnModel.meta.instrument.name = 'MIRI'
        rnModel.meta.subarray.xstart = 1
        rnModel.meta.subarray.ystart = 1
        rnModel.meta.subarray.xsize = 20
        rnModel.meta.subarray.ysize = 20
        return model1, gdq, rnModel, pixdq, err, gain

    return _setup
