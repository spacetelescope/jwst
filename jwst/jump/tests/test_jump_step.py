from itertools import cycle

import pytest
import numpy as np

from jwst.datamodels import RampModel
from jwst.datamodels import GainModel, ReadnoiseModel
from jwst.jump import JumpStep

MAXIMUM_CORES = ['none', 'quarter','half','all']

@pytest.fixture(scope="module")
def generate_miri_reffiles(tmpdir_factory):
    gainfile = str(tmpdir_factory.mktemp("data").join("gain.fits"))
    readnoisefile = str(tmpdir_factory.mktemp("data").join('readnoise.fits'))

    ingain = 6
    xsize = 103
    ysize = 102
    gain = np.ones(shape=(ysize, xsize), dtype=np.float64) * ingain
    gain_model = GainModel(data=gain)
    gain_model.meta.instrument.name = "MIRI"
    gain_model.meta.subarray.name = "FULL"
    gain_model.meta.subarray.xstart = 1
    gain_model.meta.subarray.ystart = 1
    gain_model.meta.subarray.xsize = xsize
    gain_model.meta.subarray.ysize = ysize
    gain_model.save(gainfile)

    inreadnoise = 5
    rnoise = np.ones(shape=(ysize, xsize), dtype=np.float64) * inreadnoise
    readnoise_model = ReadnoiseModel(data=rnoise)
    readnoise_model.meta.instrument.name = "MIRI"
    readnoise_model.meta.subarray.xstart = 1
    readnoise_model.meta.subarray.ystart = 1
    readnoise_model.meta.subarray.xsize = xsize
    readnoise_model.meta.subarray.ysize = ysize
    readnoise_model.save(readnoisefile)

    return gainfile, readnoisefile


@pytest.fixture(scope="module")
def generate_nircam_reffiles(tmpdir_factory):
    gainfile = str(tmpdir_factory.mktemp("ndata").join("gain.fits"))
    readnoisefile = str(tmpdir_factory.mktemp("ndata").join('readnoise.fits'))

    ingain = 6
    xsize = 20
    ysize = 20
    gain = np.ones(shape=(ysize, xsize), dtype=np.float64) * ingain
    gain_model = GainModel(data=gain)
    gain_model.meta.instrument.name = "NIRCAM"
    gain_model.meta.subarray.name = "FULL"
    gain_model.meta.subarray.xstart = 1
    gain_model.meta.subarray.ystart = 1
    gain_model.meta.subarray.xsize = xsize
    gain_model.meta.subarray.ysize = ysize
    gain_model.save(gainfile)

    inreadnoise = 5
    rnoise = np.ones(shape=(ysize, xsize), dtype=np.float64) * inreadnoise
    readnoise_model = ReadnoiseModel(data=rnoise)
    readnoise_model.meta.instrument.name = "NIRCAM"
    readnoise_model.meta.subarray.xstart = 1
    readnoise_model.meta.subarray.ystart = 1
    readnoise_model.meta.subarray.xsize = xsize
    readnoise_model.meta.subarray.ysize = ysize
    readnoise_model.save(readnoisefile)

    return gainfile, readnoisefile


@pytest.fixture
def setup_inputs():

    def _setup(ngroups=10, readnoise=10, nints=1, nrows=1024, ncols=1032,
                                nframes=1, grouptime=1.0, gain=1, deltatime=1):
        times = np.array(list(range(ngroups)), dtype=np.float64) * deltatime
        gain = np.ones(shape=(nrows, ncols), dtype=np.float64) * gain
        err = np.ones(shape=(nints, ngroups, nrows, ncols), dtype=np.float64)
        data = np.zeros(shape=(nints, ngroups, nrows, ncols), dtype=np.float64)
        pixdq = np.zeros(shape=(nrows, ncols), dtype=np.uint32)
        read_noise = np.full((nrows, ncols), readnoise, dtype=np.float64)
        gdq = np.zeros(shape=(nints, ngroups, nrows, ncols), dtype=np.uint32)

        rampmodel = RampModel(data=data, err=err, pixeldq=pixdq, groupdq=gdq, times=times)
        rampmodel.meta.instrument.name = 'MIRI'
        rampmodel.meta.instrument.detector = 'MIRIMAGE'
        rampmodel.meta.instrument.filter = 'F480M'
        rampmodel.meta.observation.date = '2015-10-13'
        rampmodel.meta.exposure.type = 'MIR_IMAGE'
        rampmodel.meta.exposure.group_time = deltatime
        rampmodel.meta.subarray.name = 'FULL'
        rampmodel.meta.subarray.xstart = 1
        rampmodel.meta.subarray.ystart = 1
        rampmodel.meta.subarray.xsize = ncols
        rampmodel.meta.subarray.ysize = nrows
        rampmodel.meta.exposure.frame_time = deltatime
        rampmodel.meta.exposure.ngroups = ngroups
        rampmodel.meta.exposure.group_time = deltatime
        rampmodel.meta.exposure.nframes = 1
        rampmodel.meta.exposure.groupgap = 0

        gain = GainModel(data=gain)
        gain.meta.instrument.name = 'MIRI'
        gain.meta.subarray.xstart = 1
        gain.meta.subarray.ystart = 1
        gain.meta.subarray.xsize = ncols
        gain.meta.subarray.ysize = nrows

        rnmodel = ReadnoiseModel(data=read_noise)
        rnmodel.meta.instrument.name = 'MIRI'
        rnmodel.meta.subarray.xstart = 1
        rnmodel.meta.subarray.ystart = 1
        rnmodel.meta.subarray.xsize = ncols
        rnmodel.meta.subarray.ysize = nrows

        return rampmodel, gdq, rnmodel, pixdq, err, gain

    return _setup


@pytest.mark.parametrize("max_cores", MAXIMUM_CORES)
def test_one_CR(generate_miri_reffiles, max_cores, setup_inputs):
    override_gain, override_readnoise = generate_miri_reffiles
    print("max_cores = ",max_cores)
    grouptime = 3.0
    deltaDN = 5
    ingain = 6
    inreadnoise = np.float64(7)
    ngroups = 100
    CR_fraction = 3
    xsize = 103
    ysize = 102
    model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=ngroups,
        nrows=ysize, ncols=xsize,
        gain=ingain, readnoise=inreadnoise, deltatime=grouptime)
    for i in range(ngroups):
        model1.data[0, i, :, :] = deltaDN * i
    first_CR_group_locs = [x for x in range(1,89) if x % 5 == 0]
    CR_locs = [x for x in range(xsize*ysize) if x % CR_fraction == 0]
    CR_x_locs = [x % ysize for x in CR_locs]
    CR_y_locs = [np.int(x / xsize) for x in CR_locs]
    CR_pool = cycle(first_CR_group_locs)
    for i in range(len(CR_x_locs)):
        CR_group = next(CR_pool)
        model1.data[0,CR_group:,CR_y_locs[i], CR_x_locs[i]] = \
            model1.data[0,CR_group:, CR_y_locs[i], CR_x_locs[i]] + 500

    print("number of CRs {}".format(len(CR_x_locs)))

    out_model = JumpStep.call(model1, override_gain=override_gain,
        override_readnoise = override_readnoise, maximum_cores = max_cores)
    CR_pool = cycle(first_CR_group_locs)
    for i in range(len(CR_x_locs)):
        CR_group = next(CR_pool)
        assert (4 == np.max(out_model.groupdq[0, CR_group, CR_y_locs[i], CR_x_locs[i]]))


@pytest.mark.parametrize("max_cores", MAXIMUM_CORES)
def test_nircam(generate_nircam_reffiles, setup_inputs, max_cores):
    override_gain, override_readnoise = generate_nircam_reffiles

    grouptime = 3.0
    deltaDN = 5
    ingain = 6
    inreadnoise = np.float64(7)
    ngroups = 100
    CR_fraction = 5
    nrows = 20
    ncols = 20
    model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=ngroups,
        nrows=nrows, ncols=ncols, gain=ingain, readnoise=inreadnoise, deltatime=grouptime)
    for i in range(ngroups):
        model1.data[0, i, :, :] = deltaDN * i
    first_CR_group_locs = [x for x in range(1,89) if x % 5 == 0]
    CR_locs = [x for x in range(nrows*ncols) if x % CR_fraction == 0]
    CR_x_locs = [x % ncols for x in CR_locs]
    CR_y_locs = [np.int(x / nrows) for x in CR_locs]
    CR_pool = cycle(first_CR_group_locs)
    for i in range(len(CR_x_locs)):
        CR_group = next(CR_pool)
        model1.data[0,CR_group:,CR_y_locs[i], CR_x_locs[i]] = \
            model1.data[0,CR_group:, CR_y_locs[i], CR_x_locs[i]] + 500

    print("number of CRs {}".format(len(CR_x_locs)))

    out_model = JumpStep.call(model1, override_gain=override_gain,
        override_readnoise=override_readnoise, maximum_cores=max_cores)
    CR_pool = cycle(first_CR_group_locs)
    for i in range(len(CR_x_locs)):
        CR_group = next(CR_pool)
        assert (4 == np.max(out_model.groupdq[0, CR_group, CR_y_locs[i], CR_x_locs[i]]))


@pytest.mark.parametrize("max_cores", MAXIMUM_CORES)
def test_two_CRs(generate_miri_reffiles, max_cores, setup_inputs):
    override_gain, override_readnoise = generate_miri_reffiles
    grouptime = 3.0
    deltaDN = 5
    ingain = 6
    inreadnoise = np.float64(7)
    ngroups = 100
    CR_fraction = 5
    xsize = 103
    ysize = 102
    model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=ngroups,
        nrows=ysize, ncols=xsize,
        gain=ingain, readnoise=inreadnoise, deltatime=grouptime)
    for i in range(ngroups):
        model1.data[0, i, :, :] = deltaDN * i
    first_CR_group_locs = [x for x in range(1,89) if x % 5 == 0]
    CR_locs = [x for x in range(xsize*ysize) if x % CR_fraction == 0]
    CR_x_locs = [x % ysize for x in CR_locs]
    CR_y_locs = [np.int(x / xsize) for x in CR_locs]
    CR_pool = cycle(first_CR_group_locs)
    for i in range(len(CR_x_locs)):
        CR_group = next(CR_pool)
        model1.data[0,CR_group:,CR_y_locs[i], CR_x_locs[i]] = \
            model1.data[0,CR_group:, CR_y_locs[i], CR_x_locs[i]] + 500
        model1.data[0, CR_group+8:, CR_y_locs[i], CR_x_locs[i]] = \
            model1.data[0, CR_group+8:, CR_y_locs[i], CR_x_locs[i]] + 700
    out_model = JumpStep.call(model1, override_gain=override_gain,
        override_readnoise=override_readnoise, maximum_cores=max_cores)
    CR_pool = cycle(first_CR_group_locs)
    for i in range(len(CR_x_locs)):
        CR_group = next(CR_pool)
        assert (4 == np.max(out_model.groupdq[0, CR_group, CR_y_locs[i], CR_x_locs[i]]))
        assert (4 == np.max(out_model.groupdq[0, CR_group+8, CR_y_locs[i], CR_x_locs[i]]))


@pytest.mark.parametrize("max_cores", MAXIMUM_CORES)
def test_two_group_integration(generate_miri_reffiles, max_cores, setup_inputs):
    override_gain, override_readnoise = generate_miri_reffiles
    grouptime = 3.0
    ingain = 6
    inreadnoise = np.float64(7)
    ngroups = 2
    xsize = 103
    ysize = 102
    model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=ngroups,
        nrows=ysize, ncols=xsize,
        gain=ingain, readnoise=inreadnoise, deltatime=grouptime)
    out_model = JumpStep.call(model1, override_gain=override_gain,
                              override_readnoise=override_readnoise, maximum_cores=max_cores)
    assert(out_model.meta.cal_step.jump == 'SKIPPED')


def test_four_group_integration(generate_miri_reffiles, setup_inputs):
    override_gain, override_readnoise = generate_miri_reffiles
    grouptime = 3.0
    ingain = 6
    inreadnoise = np.float64(7)
    ngroups = 4
    xsize = 103
    ysize = 102
    model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=ngroups,
                                                          nrows=ysize, ncols=xsize,
                                                          gain=ingain, readnoise=inreadnoise, deltatime=grouptime)
    out_model = JumpStep.call(model1, override_gain=override_gain,
                              override_readnoise=override_readnoise, maximum_cores='none')
    assert(out_model.meta.cal_step.jump == 'SKIPPED')


def test_five_group_integration(generate_miri_reffiles, setup_inputs):
    override_gain, override_readnoise = generate_miri_reffiles
    grouptime = 3.0
    ingain = 6
    inreadnoise = np.float64(7)
    ngroups = 5
    xsize = 103
    ysize = 102
    model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=ngroups,
                                                          nrows=ysize, ncols=xsize,
                                                          gain=ingain, readnoise=inreadnoise, deltatime=grouptime)
    out_model = JumpStep.call(model1, override_gain=override_gain,
                                  override_readnoise=override_readnoise, maximum_cores='none')
    assert (out_model.meta.cal_step.jump == 'COMPLETE')
