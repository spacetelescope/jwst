import numpy as np
import pytest
from astropy.io import fits
from stdatamodels.jwst.datamodels import GuiderRawModel, ImageModel, MaskModel, RampModel, dqflags

from jwst.dq_init import DQInitStep, dq_init_step
from jwst.dq_init.dq_initialization import do_dqinit
from jwst.refpix.refpix_step import collate_superstripes

# Set parameters for multiple runs of data
args = "xstart, ystart, xsize, ysize, nints, ngroups, instrument, exp_type"
test_data = [
    (1, 1, 2304, 2048, 2, 2, "FGS", "FGS_ID-STACK"),
    (1, 1, 2048, 2048, 2, 2, "FGS", "FGS_ID-IMAGE"),
    (1, 1, 2048, 2048, 2, 2, "NIRCAM", "NRC_IMAGE"),
    (1, 1, 1032, 1024, 1, 5, "MIRI", "MIR_IMAGE"),
]
ids = ["GuiderRawModel-Stack", "GuiderRawModel-Image", "RampModel", "RampModel"]


@pytest.mark.parametrize(args, test_data, ids=ids)
def test_dq_im_default(xstart, ystart, xsize, ysize, nints, ngroups, instrument, exp_type):
    """Check that PIXELDQ is initialized with the information from the reference file.
    test that a flagged value in the reference file flags the PIXELDQ array"""

    # create raw input data for step
    dm_ramp = make_rawramp(instrument, nints, ngroups, ysize, xsize, ystart, xstart, exp_type)

    # create a MaskModel for the dq input mask
    dq, dq_def = make_maskmodel(ysize, xsize)

    # edit reference file with known bad pixel values
    dq[100, 100] = 2  # Dead pixel
    dq[200, 100] = 4  # Hot pixel
    dq[300, 100] = 8  # Unreliable_slope
    dq[400, 100] = 16  # RC
    dq[500, 100] = 1  # Do_not_use
    dq[500, 101] = 1  # Do_not_use
    dq[100, 200] = 3  # Dead pixel + do not use
    dq[200, 200] = 5  # Hot pixel + do not use
    dq[300, 200] = 9  # Unreliable slope + do not use
    dq[400, 200] = 17  # RC + do not use

    # write mask model
    ref_data = MaskModel(dq=dq, dq_def=dq_def)
    ref_data.meta.instrument.name = instrument
    ref_data.meta.subarray.xstart = xstart
    ref_data.meta.subarray.xsize = xsize
    ref_data.meta.subarray.ystart = ystart
    ref_data.meta.subarray.ysize = ysize

    # user-supplied DQ array
    user_dq = np.zeros_like(dq)
    user_dq[0, 0] = 1  # Do_not_use
    user_dq[500, 101] = 2  # Dead pixel

    # run do_dqinit
    outfile = do_dqinit(dm_ramp, ref_data, user_dq=user_dq)

    if instrument == "FGS":
        dqdata = outfile.dq
    else:
        dqdata = outfile.pixeldq

    # assert that the pixels read back in match the mapping from ref data to science data
    assert dqdata[0, 0] == dqflags.pixel["DO_NOT_USE"]
    assert dqdata[100, 100] == dqflags.pixel["DEAD"]
    assert dqdata[200, 100] == dqflags.pixel["HOT"]
    assert dqdata[300, 100] == dqflags.pixel["UNRELIABLE_SLOPE"]
    assert dqdata[400, 100] == dqflags.pixel["RC"]
    assert dqdata[500, 100] == dqflags.pixel["DO_NOT_USE"]
    assert dqdata[500, 101] == 3
    assert dqdata[100, 200] == 1025
    assert dqdata[200, 200] == 2049
    assert dqdata[300, 200] == 16777217
    assert dqdata[400, 200] == 16385


def test_dq_im_wrong_shape():
    """Use one of the param combos from test_dq_im_default to trigger exception."""
    xstart = ystart = nints = 1
    xsize = 1032
    ysize = 1024
    ngroups = 5
    instrument = "MIRI"
    exp_type = "MIR_IMAGE"

    dm_ramp = make_rawramp(instrument, nints, ngroups, ysize, xsize, ystart, xstart, exp_type)
    dq, dq_def = make_maskmodel(ysize, xsize)

    ref_data = MaskModel(dq=dq, dq_def=dq_def)
    ref_data.meta.instrument.name = instrument
    ref_data.meta.subarray.xstart = xstart
    ref_data.meta.subarray.xsize = xsize
    ref_data.meta.subarray.ystart = ystart
    ref_data.meta.subarray.ysize = ysize

    with pytest.raises(ValueError, match="user_dq has shape"):
        do_dqinit(dm_ramp, ref_data, user_dq=np.ones((2, 2), dtype=dq.dtype))


def test_groupdq():
    """Check that GROUPDQ extension is added to the data and all values are initialized to zero."""

    # size of integration
    instrument = "MIRI"
    nints = 1
    ngroups = 5
    xsize = 1032
    ysize = 1024
    xstart = 1
    ystart = 1

    # create raw input data for step
    dm_ramp = make_rawramp(instrument, nints, ngroups, ysize, xsize, ystart, xstart)

    # create a MaskModel for the dq input mask
    dq, dq_def = make_maskmodel(ysize, xsize)

    # write mask model
    ref_data = MaskModel(dq=dq, dq_def=dq_def)
    ref_data.meta.instrument.name = instrument
    ref_data.meta.subarray.xstart = xstart
    ref_data.meta.subarray.xsize = xsize
    ref_data.meta.subarray.ystart = ystart
    ref_data.meta.subarray.ysize = ysize

    # run the correction step
    outfile = do_dqinit(dm_ramp, ref_data)

    # check that GROUPDQ was created and initialized to zero
    groupdq = outfile.groupdq

    np.testing.assert_array_equal(
        np.full((1, ngroups, ysize, xsize), 0, dtype=int),
        groupdq,
        err_msg="groupdq not initialized to zero",
    )


def test_dq_subarray():
    """Test that the pipeline properly extracts the subarray from the reference file."""
    # put dq flags in specific pixels and make sure they match in the output subarray file

    # create input data
    # create model of data with 0 value array
    ngroups = 50
    ysize = 224
    xsize = 288
    fullxsize = 1032
    fullysize = 1024

    # create the data and groupdq arrays
    csize = (1, ngroups, ysize, xsize)
    data = np.full(csize, 1.0)
    pixeldq = np.zeros((ysize, xsize), dtype=int)
    groupdq = np.zeros(csize, dtype=int)

    # create a JWST datamodel for MIRI data
    im = RampModel(data=data, pixeldq=pixeldq, groupdq=groupdq)

    im.meta.instrument.name = "MIRI"
    im.meta.instrument.detector = "MIRIMAGE"
    im.meta.instrument.filter = "F1500W"
    im.meta.instrument.band = "N/A"
    im.meta.observation.date = "2016-06-01"
    im.meta.observation.time = "00:00:00"
    im.meta.exposure.type = "MIR_IMAGE"
    im.meta.subarray.name = "MASK1550"
    im.meta.subarray.xstart = 1
    im.meta.subarray.xsize = xsize
    im.meta.subarray.ystart = 467
    im.meta.subarray.ysize = ysize

    # create full size mask model
    dq, dq_def = make_maskmodel(fullysize, fullxsize)

    # place dq flags in dq array that would be in subarray
    # MASK1550 file has colstart=1, rowstart=467
    dq[542, 100] = 2
    dq[550, 100] = 1
    dq[580, 80] = 4

    # write mask model
    ref_data = MaskModel(dq=dq, dq_def=dq_def)
    ref_data.meta.instrument.name = "MIRI"
    ref_data.meta.subarray.xstart = 1
    ref_data.meta.subarray.xsize = fullxsize
    ref_data.meta.subarray.ystart = 1
    ref_data.meta.subarray.ysize = fullysize
    ref_data.meta.description = "foo"
    ref_data.meta.reftype = "mask"
    ref_data.meta.author = "pytest"
    ref_data.meta.pedigree = "foo"
    ref_data.meta.useafter = "2000-01-01T00:00:00"

    # run correction step
    outfile = do_dqinit(im, ref_data)

    # read dq array
    outpixdq = outfile.pixeldq

    # check for dq flag in pixeldq of subarray image
    assert outpixdq[76, 100] == 1024
    assert outpixdq[84, 100] == 1
    assert outpixdq[114, 80] == 2048  # check that pixel was flagged 'NO_SAT_CHECK'


def test_dq_add1_groupdq():
    """
    Test if the dq_init code set the groupdq flag on the first
    group to 'do_not_use' by adding 1 to the flag, not overwriting to 1
    Also test whether two flags on the same pixel are added together.
    """

    # size of integration
    nints = 1
    ngroups = 5
    xsize = 1032
    ysize = 1024

    # create raw input data for step
    dm_ramp = make_rampmodel(nints, ngroups, ysize, xsize)

    # create a MaskModel for the dq input mask
    dq, dq_def = make_maskmodel(ysize, xsize)

    # write reference file with known bad pixel values

    dq[505, 505] = 1  # Do_not_use
    dq[400, 500] = 3  # do_not_use and dead pixel

    # write mask model
    ref_data = MaskModel(dq=dq, dq_def=dq_def)
    ref_data.meta.instrument.name = "MIRI"
    ref_data.meta.subarray.xstart = 1
    ref_data.meta.subarray.xsize = xsize
    ref_data.meta.subarray.ystart = 1
    ref_data.meta.subarray.ysize = ysize

    # set a flag in the pixel dq
    dm_ramp.pixeldq[505, 505] = 4

    # run correction step
    outfile = do_dqinit(dm_ramp, ref_data)

    # test if pixels in pixeldq were incremented in value by 1
    assert outfile.pixeldq[505, 505] == 5  # check that previous dq flag is added to mask value
    assert outfile.pixeldq[400, 500] == 1025  # check two flags propagate correctly


# Set parameters for multiple runs of guider data
args = "xstart, ystart, xsize, ysize, nints, ngroups, instrument, exp_type, detector"
test_data_multiple = [
    (1, 1, 2048, 2048, 2, 2, "FGS", "FGS_ID-IMAGE", "GUIDER1"),
    (1, 1, 1032, 1024, 1, 5, "MIRI", "MIR_IMAGE", "MIRIMAGE"),
]
ids = ["GuiderRawModel-Image", "RampModel"]


@pytest.mark.parametrize(args, test_data_multiple, ids=ids)
def test_fullstep_default(
    xstart, ystart, xsize, ysize, nints, ngroups, instrument, exp_type, detector
):
    """Test that the full step runs"""

    # create raw input data for step
    dm_ramp = make_rawramp(instrument, nints, ngroups, ysize, xsize, ystart, xstart, exp_type)

    dm_ramp.meta.instrument.name = instrument
    dm_ramp.meta.instrument.detector = detector
    dm_ramp.meta.observation.date = "2016-06-01"
    dm_ramp.meta.observation.time = "00:00:00"

    # run the full step
    outfile = DQInitStep.call(dm_ramp)

    # test that a pixeldq frame has been initialized
    if instrument == "FGS":
        assert outfile.dq.ndim == 2
    else:
        assert outfile.pixeldq.ndim == 2  # a 2-d pixeldq frame exists


def test_fullstep_userdq(tmp_path):
    """
    Test that the full step runs with user-supplied DQ array
    for one of the parameterized options of test_fullstep_default() above.
    """
    xstart = ystart = 1
    xsize = 1032
    ysize = 1024
    nints = 1
    ngroups = 5
    instrument = "MIRI"
    exp_type = "MIR_IMAGE"
    detector = "MIRIMAGE"

    # create raw input data for step
    dm_ramp = make_rawramp(instrument, nints, ngroups, ysize, xsize, ystart, xstart, exp_type)

    dm_ramp.meta.instrument.name = instrument
    dm_ramp.meta.instrument.detector = detector
    dm_ramp.meta.observation.date = "2016-06-01"
    dm_ramp.meta.observation.time = "00:00:00"

    # create a MaskModel for the dq input mask
    dq, dq_def = make_maskmodel(ysize, xsize)

    # write mask model
    ref_data = MaskModel(dq=dq, dq_def=dq_def)
    ref_data.meta.instrument.name = instrument
    ref_data.meta.subarray.xstart = xstart
    ref_data.meta.subarray.xsize = xsize
    ref_data.meta.subarray.ystart = ystart
    ref_data.meta.subarray.ysize = ysize

    # create user-defined DQ file
    dqarr = np.zeros((ysize, xsize), dtype=np.uint32)
    dqarr[100, 100] = 1
    hdu = fits.PrimaryHDU(data=dqarr)
    user_dq_file = str(tmp_path / "flags_of_our_user.fits")
    hdu.writeto(user_dq_file, overwrite=True)

    # run the full step
    outfile = DQInitStep.call(dm_ramp, override_mask=ref_data, user_supplied_dq=user_dq_file)

    # test that a pixeldq frame has been initialized
    assert outfile.pixeldq.ndim == 2  # a 2-d pixeldq frame exists

    # test that user DQ is propagated
    expected_dq = np.zeros_like(dqarr)
    expected_dq[100, 100] = 1
    np.testing.assert_array_equal(outfile.pixeldq, expected_dq)


def test_output_is_not_input():
    """Test that input data is not modified."""
    dm_ramp = make_rampmodel()
    result = DQInitStep.call(dm_ramp)
    assert result is not dm_ramp
    assert result.meta.cal_step.dq_init == "COMPLETE"
    assert dm_ramp.meta.cal_step.dq_init is None


def test_missing_mask():
    dm_ramp = make_rampmodel()
    result = DQInitStep.call(dm_ramp, override_mask="N/A")
    assert result is not dm_ramp
    assert result.meta.cal_step.dq_init == "SKIPPED"
    assert dm_ramp.meta.cal_step.dq_init is None


def test_open_rampmodel_error(monkeypatch, caplog):
    dm_ramp = make_rampmodel()

    # Mock an error in opening the model as a RampModel
    def mock_model(*args):
        raise ValueError("Test error")

    monkeypatch.setattr(dq_init_step.datamodels.RampModel, "__init__", mock_model)

    result = DQInitStep.call(dm_ramp)

    # No error raised, ramp is opened as a GuiderRawModel
    assert "Input opened as GuiderRawModel" in caplog.text
    assert isinstance(result, GuiderRawModel)


def test_open_double_error(monkeypatch, caplog):
    dm_ramp = make_rampmodel()

    # Mock an error in opening the model as a RampModel
    def mock_model(*args):
        raise ValueError("Test error")

    monkeypatch.setattr(dq_init_step.datamodels.RampModel, "__init__", mock_model)

    # And also in opening the model as a GuiderRawModel
    monkeypatch.setattr(dq_init_step.datamodels.GuiderRawModel, "__init__", mock_model)

    # Exception is raised
    with pytest.raises(ValueError, match="Test error"):
        DQInitStep.call(dm_ramp)


def test_open_unknown_error(monkeypatch, caplog):
    # Make a model that is not a RampModel, so the step will attempt to convert it
    dm_ramp = ImageModel()

    # Mock an unknown error in opening the model as a RampModel
    def mock_model(*args):
        raise RuntimeError("Test error")

    monkeypatch.setattr(dq_init_step.datamodels.RampModel, "__init__", mock_model)

    # Exception is raised
    with pytest.raises(RuntimeError, match="Test error"):
        DQInitStep.call(dm_ramp)


def test_superstripe():
    # NIRISS SOSS model with 10 superstripes
    model = make_superstripe_model()
    nstripe = model.meta.subarray.num_superstripe

    # full frame mask model with the same bad pixels in each stripe
    mask = make_superstripe_mask_model()

    # run the step
    result = DQInitStep.call(model, override_mask=mask)

    # groupdq matches the data shape
    assert result.groupdq.shape == result.data.shape

    # pixeldq matches stripe + image shape
    assert result.pixeldq.shape == (nstripe, *result.data.shape[-2:])

    # check that each stripe has the same pixels marked in pixeldq and
    # and DO_NOT_USE is also set in groupdq
    expected = one_stripe_mask(result.data.shape[-2:])
    for i in range(result.data.shape[0]):
        for j in range(result.data.shape[1]):
            np.testing.assert_equal(result.groupdq[i, j], expected & dqflags.pixel["DO_NOT_USE"])
    for i in range(nstripe):
        np.testing.assert_equal(result.pixeldq[i], expected)

    # reassemble the stripes into a full frame image
    reassembled = collate_superstripes(result)

    # check that the reassembled pixel DQ is the same as the full frame mask,
    # extracted to the subarray region, except for the edges which are
    # filled with the refpix flag
    ystart = model.meta.subarray.ystart - 1
    ystop = model.meta.subarray.ystart + model.meta.subarray.ysize
    np.testing.assert_equal(reassembled.pixeldq[:, 4:-4], mask.dq[ystart:ystop, 4:-4])
    np.testing.assert_equal(reassembled.pixeldq[:, :4], dqflags.pixel["REFERENCE_PIXEL"])
    np.testing.assert_equal(reassembled.pixeldq[:, -4:], dqflags.pixel["REFERENCE_PIXEL"])


def make_rawramp(instrument, nints, ngroups, ysize, xsize, ystart, xstart, exp_type=None):
    # create the data and groupdq arrays
    csize = (nints, ngroups, ysize, xsize)
    data = np.full(csize, 1.0)

    # create a JWST datamodel
    if instrument == "FGS":
        dm_ramp = GuiderRawModel(data=data)
        dm_ramp.meta.exposure.type = exp_type
    elif instrument == "MIRI":
        dm_ramp = RampModel(data=data)
    else:
        dm_ramp = RampModel(data=data)

    dm_ramp.meta.subarray.xstart = xstart
    dm_ramp.meta.subarray.xsize = xsize
    dm_ramp.meta.subarray.ystart = ystart
    dm_ramp.meta.subarray.ysize = ysize

    # add some basic metadata
    dm_ramp.meta.instrument.name = instrument
    dm_ramp.meta.observation.date = "2026-01-01"
    dm_ramp.meta.observation.time = "00:00:00"
    dm_ramp.meta.exposure.nints = nints
    dm_ramp.meta.exposure.ngroups = ngroups

    return dm_ramp


def make_rampmodel(nints=1, ngroups=5, ysize=1024, xsize=1032):
    # create the data and groupdq arrays
    csize = (nints, ngroups, ysize, xsize)
    data = np.full(csize, 1.0)
    pixeldq = np.zeros((ysize, xsize), dtype=int)
    groupdq = np.zeros(csize, dtype=int)

    # create a JWST datamodel for MIRI data
    dm_ramp = RampModel(data=data, pixeldq=pixeldq, groupdq=groupdq)

    dm_ramp.meta.instrument.name = "MIRI"
    dm_ramp.meta.instrument.detector = "MIRIMAGE"
    dm_ramp.meta.observation.date = "2018-01-01"
    dm_ramp.meta.observation.time = "00:00:00"
    dm_ramp.meta.subarray.xstart = 1
    dm_ramp.meta.subarray.xsize = xsize
    dm_ramp.meta.subarray.ystart = 1
    dm_ramp.meta.subarray.ysize = ysize

    return dm_ramp


def make_maskmodel(ysize, xsize):
    # create a mask model for the dq_init step
    csize = (ysize, xsize)
    dq = np.zeros(csize, dtype=int)
    # define a dq_def extension
    mask = MaskModel()

    dqdef = [
        (0, 1, "DO_NOT_USE", "Bad Pixel do not use"),
        (1, 2, "DEAD", "Dead Pixel"),
        (2, 4, "HOT", "Hot pixel"),
        (3, 8, "UNRELIABLE_SLOPE", "Large slope variance"),
        (4, 16, "RC", "RC pixel"),
        (5, 32, "REFERENCE_PIXEL", "Reference Pixel"),
    ]

    dq_def = np.array((dqdef), dtype=mask.get_dtype("dq_def"))

    return dq, dq_def


def make_superstripe_model():
    # Make a NIRISS SOSS superstripe raw file
    nints = 3
    nstripe = 10
    ngroups = 5
    ysize = 256
    xsize = 208
    xstart = 1
    ystart = 1793
    model = make_rawramp("NIRISS", nints * nstripe, ngroups, ysize, xsize, ystart, xstart)

    # Add some more basic metadata
    model.meta.instrument.detector = "NIS"
    model.meta.subarray.name = "SUB204STRIPE_SOSS"
    model.meta.subarray.fastaxis = -2
    model.meta.subarray.slowaxis = -1
    model.meta.exposure.integration_start = 1
    model.meta.exposure.integration_end = nints

    # Add the superstripe metadata
    model.meta.subarray.multistripe_reads1 = 4
    model.meta.subarray.multistripe_reads2 = 204
    model.meta.subarray.multistripe_skips1 = 0
    model.meta.subarray.multistripe_skips2 = 0
    model.meta.subarray.num_superstripe = nstripe
    model.meta.subarray.repeat_stripe = 1
    model.meta.subarray.interleave_reads1 = 1
    model.meta.subarray.superstripe_step = 204

    return model


def one_stripe_mask(shape):
    """
    Make a DQ mask with a range of bad pixels for a single stripe.

    Parameters
    ----------
    shape : tuple of int
        Stripe data shape.

    Returns
    -------
    dq : ndarray
        DQ array matching the stripe with bad pixels.
    """
    dq = np.zeros(shape, dtype=np.uint32)
    dq[2, 100] = dqflags.pixel["DEAD"]
    dq[4, 100] = dqflags.pixel["HOT"]
    dq[6, 100] = dqflags.pixel["UNRELIABLE_SLOPE"]
    dq[8, 100] = dqflags.pixel["RC"]
    dq[10, 100] = dqflags.pixel["DO_NOT_USE"]
    dq[2, 200] = dqflags.pixel["DO_NOT_USE"] + dqflags.pixel["DEAD"]
    dq[4, 200] = dqflags.pixel["DO_NOT_USE"] + dqflags.pixel["HOT"]
    dq[6, 200] = dqflags.pixel["DO_NOT_USE"] + dqflags.pixel["UNRELIABLE_SLOPE"]
    dq[8, 200] = dqflags.pixel["DO_NOT_USE"] + dqflags.pixel["RC"]
    return dq


def make_superstripe_mask_model():
    """Create a MaskModel with bad pixels in a superstripe region."""
    ysize = 2048
    xsize = 2048
    dq = np.zeros((ysize, xsize), dtype=np.uint32)

    # edit reference file with known bad pixel values
    substart = 1792
    subsize = 256
    nstripe = 10
    stripe = 204
    stripe_shape = (subsize, stripe)
    for i in range(nstripe):
        xstart = stripe * i + 4
        xstop = xstart + stripe
        dq[substart : substart + subsize, xstart:xstop] = one_stripe_mask(stripe_shape)

    # write mask model
    ref_data = MaskModel(dq=dq)
    ref_data.meta.instrument.name = "NIRISS"
    ref_data.meta.subarray.xstart = 1
    ref_data.meta.subarray.xsize = xsize
    ref_data.meta.subarray.ystart = 1
    ref_data.meta.subarray.ysize = ysize
    ref_data.meta.description = "Test description"
    ref_data.meta.reftype = "mask"
    ref_data.meta.author = "Test Author"
    ref_data.meta.pedigree = "test"
    ref_data.meta.useafter = "2024-01-01"

    return ref_data
