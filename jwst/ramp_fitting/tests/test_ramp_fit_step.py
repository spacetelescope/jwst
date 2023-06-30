import pytest
import numpy as np

from stdatamodels.jwst.datamodels import dqflags, RampModel, GainModel, ReadnoiseModel

from jwst.ramp_fitting.ramp_fit_step import RampFitStep

DELIM = "-" * 80

test_dq_flags = dqflags.pixel

GOOD = test_dq_flags["GOOD"]
DNU = test_dq_flags["DO_NOT_USE"]
JUMP = test_dq_flags["JUMP_DET"]
SAT = test_dq_flags["SATURATED"]


@pytest.fixture(scope="module")
def generate_miri_reffiles():
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

    inreadnoise = 5
    rnoise = np.ones(shape=(ysize, xsize), dtype=np.float64) * inreadnoise
    readnoise_model = ReadnoiseModel(data=rnoise)
    readnoise_model.meta.instrument.name = "MIRI"
    readnoise_model.meta.subarray.xstart = 1
    readnoise_model.meta.subarray.ystart = 1
    readnoise_model.meta.subarray.xsize = xsize
    readnoise_model.meta.subarray.ysize = ysize

    return gain_model, readnoise_model


@pytest.fixture
def setup_inputs():

    def _setup(ngroups=10, readnoise=10, nints=1, nrows=1024, ncols=1032,
               nframes=1, grouptime=1.0, gain=1, deltatime=1):
        gain = np.ones(shape=(nrows, ncols), dtype=np.float64) * gain
        err = np.ones(shape=(nints, ngroups, nrows, ncols), dtype=np.float64)
        pixdq = np.zeros(shape=(nrows, ncols), dtype=np.uint32)
        read_noise = np.full((nrows, ncols), readnoise, dtype=np.float64)
        gdq = np.zeros(shape=(nints, ngroups, nrows, ncols), dtype=np.uint32)
        int_times = np.zeros((nints,))

        rampmodel = RampModel((nints, ngroups, nrows, ncols), int_times=int_times)

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


def setup_subarray_inputs(
        nints=1, ngroups=10, nrows=1032, ncols=1024,
        subxstart=1, subxsize=1024, subystart=1, subysize=1032,
        nframes=1, grouptime=1.0, deltatime=1,
        readnoise=10, gain=1):

    data = np.zeros(shape=(nints, ngroups, subysize, subxsize), dtype=np.float32)
    err = np.ones(shape=(nints, ngroups, nrows, ncols), dtype=np.float32)
    pixdq = np.zeros(shape=(subysize, subxsize), dtype=np.uint32)
    gdq = np.zeros(shape=(nints, ngroups, subysize, subxsize), dtype=np.uint8)
    gain = np.ones(shape=(nrows, ncols), dtype=np.float64) * gain
    read_noise = np.full((nrows, ncols), readnoise, dtype=np.float32)
    times = np.array(list(range(ngroups)), dtype=np.float64) * deltatime

    model1 = RampModel(data=data, err=err, pixeldq=pixdq, groupdq=gdq, times=times)
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


def test_ramp_fit_step(generate_miri_reffiles, setup_inputs):
    """
    Create a simple input to instantiate RampFitStep and execute a call to test
    the step class and class method.
    """
    override_gain, override_readnoise = generate_miri_reffiles
    ingain, inreadnoise = 6, 7
    grouptime = 3.0
    nints, ngroups, nrows, ncols = 1, 5, 2, 2
    model, gdq, rnModel, pixdq, err, gain = setup_inputs(
        ngroups=ngroups, readnoise=inreadnoise, nints=nints, nrows=nrows,
        ncols=ncols, gain=ingain, deltatime=grouptime)

    # Add basic ramps to each pixel
    pix = [(0, 0), (0, 1), (1, 0), (1, 1)]
    base_ramp = np.array([k + 1 for k in range(ngroups)])
    ans_slopes = np.zeros(shape=(2,2))
    for k, p in enumerate(pix):
        ramp = base_ramp * (k + 1)  # A simple linear ramp
        x, y = p
        model.data[0, :, x, y] = ramp
        ans_slopes[x, y] = ramp[0] / grouptime

    # Call ramp fit through the step class
    slopes, cube_model = RampFitStep.call(
        model, override_gain=override_gain, override_readnoise=override_readnoise,
        maximum_cores="none")

    assert slopes is not None
    assert cube_model is not None

    # Test to make sure the ramps are as expected and that the step is complete
    np.testing.assert_allclose(slopes.data, ans_slopes, rtol=1e-5)
    assert slopes.meta.cal_step.ramp_fit == "COMPLETE"


def test_subarray_5groups(tmpdir_factory):
    # all pixel values are zero. So slope should be zero
    gainfile = str(tmpdir_factory.mktemp("data").join("gain.fits"))
    readnoisefile = str(tmpdir_factory.mktemp("data").join('readnoise.fits'))

    model1, gdq, rnModel, pixdq, err, gain = setup_subarray_inputs(
        ngroups=5, subxstart=10, subystart=20, subxsize=5, subysize=15, readnoise=50)
    gain.save(gainfile)
    rnModel.save(readnoisefile)

    model1.meta.exposure.ngroups = 11
    model1.data[0, 0, 12, 1] = 10.0
    model1.data[0, 1, 12, 1] = 15.0
    model1.data[0, 2, 12, 1] = 25.0
    model1.data[0, 3, 12, 1] = 33.0
    model1.data[0, 4, 12, 1] = 60.0

    # Call ramp fit through the step class
    slopes, cube_model = RampFitStep.call(
        model1, override_gain=gainfile, override_readnoise=readnoisefile,
        maximum_cores="none", save_opt=True)

    assert slopes is not None
    assert cube_model is not None

    xvalues = np.arange(5) * 1.0
    yvalues = np.array([10, 15, 25, 33, 60])
    coeff = np.polyfit(xvalues, yvalues, 1)

    np.testing.assert_allclose(slopes.data[12, 1], coeff[0], 1e-6)


def test_int_times1(generate_miri_reffiles, setup_inputs):
    # Test whether int_times table gets copied to output when it should
    override_gain, override_readnoise = generate_miri_reffiles
    ingain, inreadnoise = 6, 7
    grouptime = 3.0
    nints, ngroups, nrows, ncols = 5, 3, 2, 2
    model, gdq, rnModel, pixdq, err, gain = setup_inputs(
        ngroups=ngroups, readnoise=inreadnoise, nints=nints, nrows=nrows,
        ncols=ncols, gain=ingain, deltatime=grouptime)

    # Set TSOVISIT false, despite which the int_times table should come back populated
    model.meta.visit.tsovisit = False

    # Call ramp fit through the step class
    slopes, cube_model = RampFitStep.call(
        model, override_gain=override_gain, override_readnoise=override_readnoise,
        maximum_cores="none")

    assert slopes is not None
    assert cube_model is not None

    assert len(cube_model.int_times) == 5


def test_int_times2(generate_miri_reffiles, setup_inputs):
    # Test whether int_times table gets copied to output when it should
    override_gain, override_readnoise = generate_miri_reffiles
    ingain, inreadnoise = 6, 7
    grouptime = 3.0
    nints, ngroups, nrows, ncols = 5, 3, 2, 2
    model, gdq, rnModel, pixdq, err, gain = setup_inputs(
        ngroups=ngroups, readnoise=inreadnoise, nints=nints, nrows=nrows,
        ncols=ncols, gain=ingain, deltatime=grouptime)

    # Set TSOVISIT true, in which case the int_times table should come back with all content
    model.meta.visit.tsovisit = True

    # Call ramp fit through the step class
    slopes, cube_model = RampFitStep.call(
        model, override_gain=override_gain, override_readnoise=override_readnoise,
        maximum_cores="none")

    assert slopes is not None
    assert cube_model is not None

    assert len(cube_model.int_times) == nints


def one_group_suppressed(nints, suppress, setup_inputs):
    """
    Creates three pixel ramps.
    The first ramp has no good groups.
    The second ramp has one good group.
    The third ramp has all good groups.

    Sets up the models to be used by the tests for the one
    group suppression flag.
    """
    # Define the data.
    ngroups, nrows, ncols = 5, 1, 3
    dims = nints, ngroups, nrows, ncols
    rnoise, gain = 10, 1
    group_time, frame_time = 5.0, 1
    rampmodel, gdq, rnModel, pixdq, err, gmodel = setup_inputs(
        ngroups=ngroups, readnoise=rnoise, nints=nints, nrows=nrows,
        ncols=ncols, gain=gain, deltatime=group_time)

    rampmodel.meta.exposure.frame_time = frame_time

    # Setup the ramp data and DQ.
    arr = np.array([k + 1 for k in range(ngroups)], dtype=float)
    sat = dqflags.pixel["SATURATED"]
    sat_dq = np.array([sat] * ngroups, dtype=rampmodel.groupdq.dtype)
    zdq = np.array([0] * ngroups, dtype=rampmodel.groupdq.dtype)

    rampmodel.data[0, :, 0, 0] = arr
    rampmodel.data[0, :, 0, 1] = arr
    rampmodel.data[0, :, 0, 2] = arr

    rampmodel.groupdq[0, :, 0, 0] = sat_dq  # All groups sat
    rampmodel.groupdq[0, :, 0, 1] = sat_dq  # 0th good, all others sat
    rampmodel.groupdq[0, 0, 0, 1] = 0
    rampmodel.groupdq[0, :, 0, 2] = zdq     # All groups good

    if nints > 1:
        rampmodel.data[1, :, 0, 0] = arr
        rampmodel.data[1, :, 0, 1] = arr
        rampmodel.data[1, :, 0, 2] = arr

        # All good ramps
        rampmodel.groupdq[1, :, 0, 0] = zdq
        rampmodel.groupdq[1, :, 0, 1] = zdq
        rampmodel.groupdq[1, :, 0, 2] = zdq

    rampmodel.suppress_one_group_ramps = suppress

    # Call ramp fit through the step class
    slopes, cube_model = RampFitStep.call(
        rampmodel,
        override_gain=gmodel,
        override_readnoise=rnModel,
        suppress_one_group=suppress,
        maximum_cores="none")

    return slopes, cube_model, dims


def test_one_group_not_suppressed_one_integration(setup_inputs):
    """
    This tests a one integration with three pixel ramps, with the
    one group suppression switch turned off.
    1. The fully saturated ramp should have no computed data with the
       DO_NOT_USE flag set in the final DQ.
    2. The one group ramp will have data computed and the SATURATED
       flag set in the final DQ.
    3. The good ramp is fully computed.
    """
    slopes, cube, dims = one_group_suppressed(1, False, setup_inputs)
    nints, ngroups, nrows, ncols = dims
    tol = 1e-5

    # Check slopes information
    check = np.array([[np.nan, .2, 0.2]])
    np.testing.assert_allclose(slopes.data, check, tol)

    check = np.array([[DNU | SAT, GOOD, GOOD]])
    np.testing.assert_allclose(slopes.dq, check, tol)

    check = np.array([[0., 0.04, 0.01]])
    np.testing.assert_allclose(slopes.var_poisson, check, tol)

    check = np.array([[0., 3.9999995, 0.19999999]])
    np.testing.assert_allclose(slopes.var_rnoise, check, tol)

    check = np.array([[0., 2.009975, 0.45825756]])
    np.testing.assert_allclose(slopes.err, check, tol)

    # Check slopes information
    check = np.array([[[np.nan, .2, .2]]])
    np.testing.assert_allclose(cube.data, check, tol)

    check = np.array([[[DNU | SAT, GOOD, GOOD]]])
    np.testing.assert_allclose(cube.dq, check, tol)

    check = np.array([[[0., 0.04, 0.01]]])
    np.testing.assert_allclose(cube.var_poisson, check, tol)

    check = np.array([[[0., 3.9999995, 0.19999999]]])
    np.testing.assert_allclose(cube.var_rnoise, check, tol)

    check = np.array([[[0., 2.0099752, 0.4582576]]])
    np.testing.assert_allclose(cube.err, check, tol)


def test_one_group_suppressed_one_integration(setup_inputs):
    """
    This tests a one integration with three pixel ramps, with the
    one group suppression switch turned on.
    1. The fully saturated ramp should have no computed data with the
       DO_NOT_USE flag set in the final DQ.
    2. The one group ramp will be computed as a fully saturated ramp,
       but the DO_NOT_USE flag will not be set in the final DQ.
    3. The good ramp is fully computed.
    """
    slopes, cube, dims = one_group_suppressed(1, True, setup_inputs)
    nints, ngroups, nrows, ncols = dims
    tol = 1e-5

    # Check slopes information
    check = np.array([[np.nan, np.nan, .2]])
    np.testing.assert_allclose(slopes.data, check, tol)

    check = np.array([[DNU | SAT, DNU, GOOD]])
    np.testing.assert_allclose(slopes.dq, check, tol)

    check = np.array([[0., 0., 0.01]])
    np.testing.assert_allclose(slopes.var_poisson, check, tol)

    check = np.array([[0., 0., 0.19999999]])
    np.testing.assert_allclose(slopes.var_rnoise, check, tol)

    check = np.array([[0., 0., 0.45825756]])
    np.testing.assert_allclose(slopes.err, check, tol)

    # Check slopes information
    check = np.array([[[np.nan, np.nan, .2]]])
    np.testing.assert_allclose(cube.data, check, tol)

    check = np.array([[[DNU | SAT, DNU, 0]]])
    np.testing.assert_allclose(cube.dq, check, tol)

    check = np.array([[[0., 0., 0.01]]])
    np.testing.assert_allclose(cube.var_poisson, check, tol)

    check = np.array([[[0., 0., 0.19999999]]])
    np.testing.assert_allclose(cube.var_rnoise, check, tol)

    check = np.array([[[0., 0., 0.4582576]]])
    np.testing.assert_allclose(cube.err, check, tol)


def test_one_group_not_suppressed_two_integrations(setup_inputs):
    """
    This tests three pixel ramps with two integrations and the
    one group suppression switch turned off.  The second integration
    for all ramps are good, so all pixels will have usable data in
    the final image model.
    1. A fully saturated first integration.  Usable data is obtained
       from the second integration.
    2. A one group ramp in the first integration.  Data from both
       integrations is combeined for the final image model.
    3. Good data from both integrations.
    """
    slopes, cube, dims = one_group_suppressed(2, False, setup_inputs)
    nints, ngroups, nrows, ncols = dims
    tol = 1e-5

    # Check slopes information
    check = np.array([[.2, .2, .2]])
    np.testing.assert_allclose(slopes.data, check, tol)

    check = np.array([[GOOD, GOOD, GOOD]])
    np.testing.assert_allclose(slopes.dq, check, tol)

    check = np.array([[0.005, 0.008, 0.005]])
    np.testing.assert_allclose(slopes.var_poisson, check, tol)

    check = np.array([[0.19999999, 0.19047618, 0.09999999]])
    np.testing.assert_allclose(slopes.var_rnoise, check, tol)

    check = np.array([[0.45276925, 0.44550666, 0.32403702]])
    np.testing.assert_allclose(slopes.err, check, tol)

    # Check slopes information
    check = np.array([[[np.nan, .2, .2]],
                      [[.2,     .2, .2]]])
    np.testing.assert_allclose(cube.data, check, tol)

    check = np.array([[[SAT | DNU, GOOD, GOOD]],
                      [[GOOD, GOOD, GOOD]]])
    np.testing.assert_allclose(cube.dq, check, tol)

    check = np.array([[[0.,    0.04, 0.01]],
                      [[0.005, 0.01, 0.01]]])
    np.testing.assert_allclose(cube.var_poisson, check, tol)

    check = np.array([[[0.,         3.9999995,  0.19999999]],
                      [[0.19999999, 0.19999999, 0.19999999]]])
    np.testing.assert_allclose(cube.var_rnoise, check, tol)

    check = np.array([[[0.,         2.0099752, 0.4582576]],
                      [[0.45276922, 0.4582576, 0.4582576]]])
    np.testing.assert_allclose(cube.err, check, tol)


def test_one_group_suppressed_two_integrations(setup_inputs):
    """
    This tests three pixel ramps with two integrations and the
    one group suppression switch turned on.  The key differences
    for this test are in the cube data.
    1. A fully saturated first integration.  Usable data is obtained
       from the second integration.
    2. A one group ramp in the first integration.  No data from the
       first integration is used, but DO_NOT_USE flag is not set for
       the first integration.
    3. Good data from both integrations.
    """
    slopes, cube, dims = one_group_suppressed(2, True, setup_inputs)
    nints, ngroups, nrows, ncols = dims
    tol = 1e-5

    # Check slopes information
    check = np.array([[.2, .2, .2]])
    np.testing.assert_allclose(slopes.data, check, tol)

    check = np.array([[GOOD, GOOD, GOOD]])
    np.testing.assert_allclose(slopes.dq, check, tol)

    check = np.array([[0.005, 0.01, 0.005]])
    np.testing.assert_allclose(slopes.var_poisson, check, tol)

    check = np.array([[0.19999999, 0.19999999, 0.09999999]])
    np.testing.assert_allclose(slopes.var_rnoise, check, tol)

    check = np.array([[0.45276925, 0.45825756, 0.32403702]])
    np.testing.assert_allclose(slopes.err, check, tol)

    # Check slopes information
    check = np.array([[[np.nan, np.nan, .2]],
                      [[.2, .2, .2]]])
    np.testing.assert_allclose(cube.data, check, tol)

    check = np.array([[[DNU | SAT, DNU, GOOD]],
                      [[GOOD, GOOD, GOOD]]])
    np.testing.assert_allclose(cube.dq, check, tol)

    check = np.array([[[0.,    0.,    0.01]],
                      [[0.005, 0.01, 0.01]]])
    np.testing.assert_allclose(cube.var_poisson, check, tol)

    check = np.array([[[0.,         0.,         0.19999999]],
                      [[0.19999999, 0.19999999, 0.19999999]]])
    np.testing.assert_allclose(cube.var_rnoise, check, tol)

    check = np.array([[[0.,         0.,         0.4582576]],
                      [[0.45276922, 0.4582576, 0.4582576]]])
    np.testing.assert_allclose(cube.err, check, tol)
