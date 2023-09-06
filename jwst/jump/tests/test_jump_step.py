from itertools import cycle

import numpy as np
import pytest

from stdatamodels.jwst.datamodels import GainModel, ReadnoiseModel, RampModel

from jwst.jump import JumpStep

MAXIMUM_CORES = ['2', 'none', 'quarter', 'half', 'all']


@pytest.fixture(scope="module")
def generate_miri_reffiles(tmpdir_factory):

    def _generate_miri_reffiles(xsize=103, ysize=102, ingain=6):

        gainfile = str(tmpdir_factory.mktemp("data").join("gain.fits"))
        readnoisefile = str(tmpdir_factory.mktemp("data").join('readnoise.fits'))

        ingain = ingain
        xsize = xsize
        ysize = ysize
        gain = np.ones(shape=(ysize, xsize), dtype=np.float64) * ingain
        gain_model = GainModel(data=gain)
        gain_model.meta.instrument.name = "MIRI"
        gain_model.meta.subarray.name = "FULL"
        gain_model.meta.subarray.xstart = 1
        gain_model.meta.subarray.ystart = 1
        gain_model.meta.subarray.xsize = xsize
        gain_model.meta.subarray.ysize = ysize
        gain_model.save(gainfile)
        gain_model.close()

        inreadnoise = 5
        rnoise = np.ones(shape=(ysize, xsize), dtype=np.float64) * inreadnoise
        readnoise_model = ReadnoiseModel(data=rnoise)
        readnoise_model.meta.instrument.name = "MIRI"
        readnoise_model.meta.subarray.xstart = 1
        readnoise_model.meta.subarray.ystart = 1
        readnoise_model.meta.subarray.xsize = xsize
        readnoise_model.meta.subarray.ysize = ysize
        readnoise_model.save(readnoisefile)
        readnoise_model.close()

        return gainfile, readnoisefile

    return _generate_miri_reffiles


@pytest.fixture(scope="module")
def generate_nircam_reffiles(tmpdir_factory):

    def _generate_nircam_reffiles(xsize=20, ysize=20, ingain=6):
        gainfile = str(tmpdir_factory.mktemp("ndata").join("gain.fits"))
        readnoisefile = str(tmpdir_factory.mktemp("ndata").join('readnoise.fits'))

        ingain = ingain
        xsize = xsize
        ysize = ysize
        gain = np.ones(shape=(ysize, xsize), dtype=np.float64) * ingain
        gain_model = GainModel(data=gain)
        gain_model.meta.instrument.name = "NIRCAM"
        gain_model.meta.subarray.name = "FULL"
        gain_model.meta.subarray.xstart = 1
        gain_model.meta.subarray.ystart = 1
        gain_model.meta.subarray.xsize = xsize
        gain_model.meta.subarray.ysize = ysize
        gain_model.save(gainfile)
        gain_model.close()

        inreadnoise = 5
        rnoise = np.ones(shape=(ysize, xsize), dtype=np.float64) * inreadnoise
        readnoise_model = ReadnoiseModel(data=rnoise)
        readnoise_model.meta.instrument.name = "NIRCAM"
        readnoise_model.meta.subarray.xstart = 1
        readnoise_model.meta.subarray.ystart = 1
        readnoise_model.meta.subarray.xsize = xsize
        readnoise_model.meta.subarray.ysize = ysize
        readnoise_model.save(readnoisefile)
        readnoise_model.close()

        return gainfile, readnoisefile

    return _generate_nircam_reffiles


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
        rampmodel.meta.observation.date = '2023-01-13'
        rampmodel.meta.observation.time = '00:00:00'
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


def add_circles_to_data(data, center_coords, radii, fill_val=None):
    """Modifies `data` to add circles at positions specified by `center_coords`
       of sizes specified by `radii`. The magnitude of each circle is 10x its
       radius (bigger snowballs are brighter).
    """

    X, Y = np.ogrid[:data.shape[0], :data.shape[1]]

    for i, _ in enumerate(radii):

        radius = radii[i]
        x_center = center_coords[i][0]
        y_center = center_coords[i][1]
        dist_from_center = np.sqrt((X - x_center)**2 + (Y - y_center)**2)
        circular_mask = ~(dist_from_center >= radius)
        if fill_val is None:
            data[circular_mask] += 10 * radius
        else:
            data[circular_mask] = fill_val


@pytest.mark.parametrize("max_cores", MAXIMUM_CORES)
def test_one_CR(generate_miri_reffiles, max_cores, setup_inputs):
    override_gain, override_readnoise = generate_miri_reffiles()
    print("max_cores = ", max_cores)
    grouptime = 3.0
    deltaDN = 5
    ingain = 6
    inreadnoise = 7
    ngroups = 100
    CR_fraction = 3
    xsize = 103
    ysize = 102
    model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=ngroups,
                                                          nrows=ysize, ncols=xsize,
                                                          gain=ingain, readnoise=inreadnoise,
                                                          deltatime=grouptime)
    for i in range(ngroups):
        model1.data[0, i, :, :] = deltaDN * i
    first_CR_group_locs = [x for x in range(1, 89) if x % 5 == 0]
    CR_locs = [x for x in range(xsize * ysize) if x % CR_fraction == 0]
    CR_x_locs = [x % ysize for x in CR_locs]
    CR_y_locs = [int(x / xsize) for x in CR_locs]
    CR_pool = cycle(first_CR_group_locs)
    for i in range(len(CR_x_locs)):
        CR_group = next(CR_pool)
        model1.data[0, CR_group:, CR_y_locs[i], CR_x_locs[i]] = \
            model1.data[0, CR_group:, CR_y_locs[i], CR_x_locs[i]] + 500

    print("number of CRs {}".format(len(CR_x_locs)))

    out_model = JumpStep.call(model1, override_gain=override_gain,
                              override_readnoise=override_readnoise, maximum_cores=max_cores)
    CR_pool = cycle(first_CR_group_locs)
    for i in range(len(CR_x_locs)):
        CR_group = next(CR_pool)
        assert 4 == np.max(out_model.groupdq[0, CR_group, CR_y_locs[i], CR_x_locs[i]])


@pytest.mark.parametrize("max_cores", MAXIMUM_CORES)
def test_nircam(generate_nircam_reffiles, setup_inputs, max_cores):
    override_gain, override_readnoise = generate_nircam_reffiles()

    grouptime = 3.0
    deltaDN = 5
    ingain = 6
    inreadnoise = 7
    ngroups = 100
    CR_fraction = 5
    nrows = 20
    ncols = 20
    model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=ngroups,
                                                          nrows=nrows, ncols=ncols,
                                                          gain=ingain, readnoise=inreadnoise,
                                                          deltatime=grouptime)
    for i in range(ngroups):
        model1.data[0, i, :, :] = deltaDN * i
    first_CR_group_locs = [x for x in range(1, 89) if x % 5 == 0]
    CR_locs = [x for x in range(nrows * ncols) if x % CR_fraction == 0]
    CR_x_locs = [x % ncols for x in CR_locs]
    CR_y_locs = [int(x / nrows) for x in CR_locs]
    CR_pool = cycle(first_CR_group_locs)
    for i in range(len(CR_x_locs)):
        CR_group = next(CR_pool)
        model1.data[0, CR_group:, CR_y_locs[i], CR_x_locs[i]] = \
            model1.data[0, CR_group:, CR_y_locs[i], CR_x_locs[i]] + 500

    print("number of CRs {}".format(len(CR_x_locs)))

    out_model = JumpStep.call(model1, override_gain=override_gain,
                              override_readnoise=override_readnoise, maximum_cores=max_cores)
    CR_pool = cycle(first_CR_group_locs)
    for i in range(len(CR_x_locs)):
        CR_group = next(CR_pool)
        assert 4 == np.max(out_model.groupdq[0, CR_group, CR_y_locs[i], CR_x_locs[i]])


@pytest.mark.parametrize("max_cores", MAXIMUM_CORES)
def test_two_CRs(generate_miri_reffiles, max_cores, setup_inputs):
    override_gain, override_readnoise = generate_miri_reffiles()
    grouptime = 3.0
    deltaDN = 5
    ingain = 6
    inreadnoise = 7
    ngroups = 100
    CR_fraction = 5
    xsize = 103
    ysize = 102
    model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=ngroups,
                                                          nrows=ysize, ncols=xsize,
                                                          gain=ingain, readnoise=inreadnoise,
                                                          deltatime=grouptime)
    for i in range(ngroups):
        model1.data[0, i, :, :] = deltaDN * i
    first_CR_group_locs = [x for x in range(1, 89) if x % 5 == 0]
    CR_locs = [x for x in range(xsize * ysize) if x % CR_fraction == 0]
    CR_x_locs = [x % ysize for x in CR_locs]
    CR_y_locs = [int(x / xsize) for x in CR_locs]
    CR_pool = cycle(first_CR_group_locs)
    for i in range(len(CR_x_locs)):
        CR_group = next(CR_pool)
        model1.data[0, CR_group:, CR_y_locs[i], CR_x_locs[i]] = \
            model1.data[0, CR_group:, CR_y_locs[i], CR_x_locs[i]] + 500
        model1.data[0, CR_group + 8:, CR_y_locs[i], CR_x_locs[i]] = \
            model1.data[0, CR_group + 8:, CR_y_locs[i], CR_x_locs[i]] + 700
    out_model = JumpStep.call(model1, override_gain=override_gain,
                              override_readnoise=override_readnoise, maximum_cores=max_cores)
    CR_pool = cycle(first_CR_group_locs)
    for i in range(len(CR_x_locs)):
        CR_group = next(CR_pool)
        assert 4 == np.max(out_model.groupdq[0, CR_group, CR_y_locs[i], CR_x_locs[i]])
        assert 4 == np.max(out_model.groupdq[0, CR_group + 8, CR_y_locs[i], CR_x_locs[i]])


@pytest.mark.parametrize("max_cores", MAXIMUM_CORES)
def test_two_group_integration(generate_miri_reffiles, max_cores, setup_inputs):
    override_gain, override_readnoise = generate_miri_reffiles()
    grouptime = 3.0
    ingain = 6
    inreadnoise = 7
    ngroups = 2
    xsize = 103
    ysize = 102
    model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=ngroups,
                                                          nrows=ysize, ncols=xsize,
                                                          gain=ingain, readnoise=inreadnoise,
                                                          deltatime=grouptime)
    out_model = JumpStep.call(model1, override_gain=override_gain,
                              override_readnoise=override_readnoise, maximum_cores=max_cores)
    assert out_model.meta.cal_step.jump == 'SKIPPED'


def test_three_group_integration(generate_miri_reffiles, setup_inputs):
    override_gain, override_readnoise = generate_miri_reffiles()
    grouptime = 3.0
    ingain = 6
    inreadnoise = 7
    ngroups = 3
    xsize = 103
    ysize = 102
    model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=ngroups,
                                                          nrows=ysize, ncols=xsize,
                                                          gain=ingain, readnoise=inreadnoise,
                                                          deltatime=grouptime)
    out_model = JumpStep.call(model1, override_gain=override_gain,
                              override_readnoise=override_readnoise, maximum_cores='none')
    assert out_model.meta.cal_step.jump == 'COMPLETE'


def test_snowball_flagging_nosat(generate_nircam_reffiles, setup_inputs):
    """Test that snowballs are properly flagged when the `sat_required_snowball`,
    which requires there to be a saturated pixel within the cluster of
    jump-flagged pixels for a snowball to be flagged, is set to FALSE.

    This test also tests if the `expand_factor` keyword, which controls how many
    pixels out a snowball is grown, is being used properly"""

    # make datamodel
    override_gain, override_readnoise = generate_nircam_reffiles(xsize=100,
                                                                 ysize=100)
    mod, _, _, _, _, _ = setup_inputs(ngroups=5, nrows=100, ncols=100, gain=6,
                                      readnoise=7, deltatime=3.0)

    # add 'snowballs' to data array in 0th integration, 1st read. when run though
    # jump step, the DQ array in the 1st groupdq group should have clusters of
    # jump-flagged pixels, which should then be detected as snowballs in that group.

    center_coords = [(20, 20), (60, 60)]
    radii = [10, 20]
    add_circles_to_data(mod.data[0, 1], center_coords, radii)

    expand_factor = 2
    jump_result = JumpStep.call(mod, override_gain=override_gain,
                                override_readnoise=override_readnoise,
                                expand_large_events=True, sat_required_snowball=False,
                                expand_factor=expand_factor)

    # both clusters should be detected as jumps then flagged as snowballs,
    # resulting in a circle of x2 radius of the original having a jump flag in
    # that group (and only that group). verify this is the case
    for i, coord in enumerate(center_coords):
        x, y = coord
        rad = radii[i]
        expanded_rad = expand_factor * rad

        initial_area = np.sum((mod.data[0, 1, y - rad: y + rad, x - rad: x + rad]).astype(bool))
        expanded_area = np.sum((jump_result.groupdq[0, 1,
                                                    y - expanded_rad: y + expanded_rad,
                                                    x - expanded_rad: x + expanded_rad]).astype(bool))

        assert (np.floor(expanded_area / initial_area) == (expand_factor**2))


def test_snowball_flagging_sat(generate_nircam_reffiles, setup_inputs):
    """Test that snowballs are properly flagged when the `sat_required_snowball`,
    which requires there to be a saturated pixel within the cluster of
    jump-flagged pixels for a snowball to be flagged, is set to FALSE.

    This test also tests if the `expand_factor` keyword, which controls how many
    pixels out a snowball is grown, is being used properly"""

    # make datamodel
    override_gain, override_readnoise = generate_nircam_reffiles(xsize=100,
                                                                 ysize=100)
    mod, _, _, _, _, _ = setup_inputs(ngroups=5, nrows=100, ncols=100, gain=6,
                                      readnoise=7, deltatime=3.0)

    # add 'snowballs' to data array in 0th integration, 1st read. when run though
    # jump step, the DQ array in the 1st groupdq group should have clusters of
    # jump-flagged pixels, which should then be detected as snowballs in that group.

    center_coords = [(20, 20), (60, 60)]
    radii = [10, 20]
    add_circles_to_data(mod.data[0, 1], center_coords, radii)

    # add circles of saturation flags in the input groupdq in the center of each snowball
    # in every read except the 0th

    for i in range(1, mod.data.shape[1]):
        add_circles_to_data(mod.groupdq[0, i], center_coords, [4, 4], fill_val=2)

    expand_factor = 2
    jump_result = JumpStep.call(mod, override_gain=override_gain,
                                override_readnoise=override_readnoise,
                                expand_large_events=True, sat_required_snowball=True,
                                expand_factor=expand_factor)

    # both clusters should be detected as jumps then flagged as snowballs,
    # resulting in a circle of x2 radius of the original having a jump flag in
    # that group (and only that group). verify this is the case
    for i, coord in enumerate(center_coords):
        x, y = coord
        rad = radii[i]
        expanded_rad = expand_factor * rad

        initial_area = np.sum((mod.data[0, 1, y - rad: y + rad, x - rad: x + rad]).astype(bool))
        expanded_area = np.sum((jump_result.groupdq[0, 1,
                                                    y - expanded_rad: y + expanded_rad,
                                                    x - expanded_rad: x + expanded_rad]).astype(bool))

    assert (np.floor(expanded_area / initial_area) == (expand_factor**2))
