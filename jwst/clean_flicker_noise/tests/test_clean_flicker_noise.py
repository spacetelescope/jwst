import warnings

import gwcs
import numpy as np
import pytest
from astropy.utils.data import get_pkg_data_filename
from numpy.testing import assert_allclose
from stdatamodels.jwst import datamodels

from jwst.assign_wcs.tests.test_nirspec import create_nirspec_ifu_file, create_nirspec_fs_file
from jwst.msaflagopen.tests.test_msa_open import make_nirspec_mos_model
from jwst.clean_flicker_noise import clean_flicker_noise as cfn


def add_metadata(model, shape):
    model.meta.instrument.name = "MIRI"
    model.meta.instrument.detector = "MIRIMAGE"
    model.meta.instrument.filter = "F480M"
    model.meta.observation.date = "2015-10-13"
    model.meta.observation.time = "00:00:00"
    model.meta.exposure.type = "MIR_IMAGE"
    model.meta.exposure.group_time = 1.0
    model.meta.subarray.name = "FULL"
    model.meta.subarray.xstart = 1
    model.meta.subarray.ystart = 1
    model.meta.subarray.xsize = shape[3]
    model.meta.subarray.ysize = shape[2]
    model.meta.exposure.frame_time = 1.0
    model.meta.exposure.ngroups = shape[1]
    model.meta.exposure.group_time = 1.0
    model.meta.exposure.nints = shape[0]
    model.meta.exposure.nframes = 1
    model.meta.exposure.groupgap = 0
    model.meta.exposure.readpatt = "FASTR1"
    model.meta.subarray.slowaxis = 2


def make_small_ramp_model(shape=(3, 5, 10, 10)):
    rampmodel = datamodels.RampModel(shape)
    add_metadata(rampmodel, shape)

    # Make data with a constant rate
    for group in range(shape[1]):
        rampmodel.data[:, group, :, :] = group

    return rampmodel


def make_small_rate_model(shape=(3, 5, 10, 10)):
    ratemodel = datamodels.ImageModel(shape[2:])
    add_metadata(ratemodel, shape)
    ratemodel.data[:] = 1.0
    return ratemodel


def make_small_rateints_model(shape=(3, 5, 10, 10)):
    ratemodel = datamodels.CubeModel((shape[0], shape[2], shape[3]))
    add_metadata(ratemodel, shape)
    ratemodel.data[:] = 1.0
    return ratemodel


def make_flat_model(model, shape=(10, 10), value=None):
    # make a flat model with appropriate size and metadata
    flat = datamodels.FlatModel()
    if value is None:
        flat.data = np.arange(shape[0] * shape[1], dtype=float).reshape(shape)
    else:
        flat.data = np.full(shape, value)

    # add required metadata
    flat.meta.description = "test"
    flat.meta.reftype = "test"
    flat.meta.author = "test"
    flat.meta.pedigree = "test"
    flat.meta.useafter = "test"

    # copy any other matching metadata
    flat.update(model)

    # make sure shape keys match input
    flat.meta.subarray.xsize = shape[1]
    flat.meta.subarray.ysize = shape[0]

    return flat


def make_nirspec_ifu_model(shape=(2048, 2048)):
    hdul = create_nirspec_ifu_file(
        grating="PRISM", filter="CLEAR", gwa_xtil=0.35986012, gwa_ytil=0.13448857, gwa_tilt=37.1
    )
    hdul["SCI"].data = np.ones(shape, dtype=float)
    rate_model = datamodels.IFUImageModel(hdul)
    hdul.close()
    return rate_model


def make_nirspec_mos_fs_model():
    mos_model = make_nirspec_mos_model()
    mos_model.meta.instrument.msa_metadata_file = get_pkg_data_filename(
        "data/msa_fs_configuration.fits", package="jwst.assign_wcs.tests"
    )
    return mos_model


def make_nirspec_fs_model():
    hdul = create_nirspec_fs_file(grating="G140M", filter="F100LP")
    hdul["SCI"].data = np.ones((2048, 2048), dtype=float)
    rate_model = datamodels.ImageModel(hdul)
    hdul.close()

    # add the slow axis
    rate_model.meta.subarray.slowaxis = 1
    return rate_model


class MockUpdate:
    def __init__(self):
        self.seen = False

    def __call__(self, input_model, mask):
        self.seen = True
        return mask


def test_make_rate(log_watcher):
    shape = (3, 5, 10, 10)
    ramp_model = make_small_ramp_model(shape)

    watcher = log_watcher(
        "jwst.clean_flicker_noise.clean_flicker_noise", message="Creating draft rate file"
    )
    result = cfn.make_rate(ramp_model, return_cube=True)
    assert isinstance(result, datamodels.CubeModel)
    assert result.data.shape == (shape[0], shape[2], shape[3])
    assert np.all(result.data == 1.0)

    result = cfn.make_rate(ramp_model, return_cube=False)
    assert isinstance(result, datamodels.ImageModel)
    assert result.data.shape == shape[2:]
    assert np.all(result.data == 1.0)

    ramp_model.close()
    result.close()

    # Check for expected log message
    watcher.assert_seen()


def test_postprocess_rate_nirspec(log_watcher):
    rate_model = make_nirspec_mos_model()

    watcher = log_watcher("jwst.clean_flicker_noise.clean_flicker_noise", message="Assigning a WCS")
    result = cfn.post_process_rate(rate_model, assign_wcs=True)
    assert isinstance(result.meta.wcs, gwcs.WCS)
    assert np.sum(result.dq & datamodels.dqflags.pixel["MSA_FAILED_OPEN"]) == 0
    watcher.assert_seen()

    watcher.message = "Flagging failed-open"
    result = cfn.post_process_rate(result, msaflagopen=True)
    assert np.sum(result.dq & datamodels.dqflags.pixel["MSA_FAILED_OPEN"]) > 0
    watcher.assert_seen()

    rate_model.close()
    result.close()


def test_postprocess_rate_miri(log_watcher):
    rate_model = make_small_rate_model()
    assert np.sum(rate_model.dq & datamodels.dqflags.pixel["DO_NOT_USE"]) == 0

    watcher = log_watcher(
        "jwst.clean_flicker_noise.clean_flicker_noise", message="Retrieving flat DQ"
    )
    result = cfn.post_process_rate(rate_model, flat_dq=True)
    watcher.assert_seen()
    assert np.sum(result.dq & datamodels.dqflags.pixel["DO_NOT_USE"]) > 0
    assert np.all(result.data == rate_model.data)

    rate_model.close()
    result.close()


def test_mask_ifu_slices():
    rate_model = make_nirspec_ifu_model()

    rate_model = cfn.post_process_rate(rate_model, assign_wcs=True)
    mask = np.full_like(rate_model.data, True)

    # Mark IFU science regions as False: about 10% of the array
    cfn.mask_ifu_slices(rate_model, mask)
    assert np.allclose(np.sum(mask), 0.9 * mask.size, atol=0.01 * mask.size)

    rate_model.close()


@pytest.mark.parametrize("exptype,blocked", [("mos", 0.012), ("mos_fs", 0.025), ("fs", 0.05)])
def test_mask_slits(exptype, blocked):
    if exptype == "mos":
        rate_model = make_nirspec_mos_model()
    elif exptype == "mos_fs":
        rate_model = make_nirspec_mos_fs_model()
    else:
        rate_model = make_nirspec_fs_model()

    # Assign a WCS
    rate_model = cfn.post_process_rate(rate_model, assign_wcs=True)

    # Mark slits as False
    mask = np.full_like(rate_model.data, True)
    cfn.mask_slits(rate_model, mask)

    # Check that the fraction of the array blocked is as expected
    assert_allclose(np.sum(mask) / mask.size, 1 - blocked, atol=0.001)

    rate_model.close()


@pytest.mark.parametrize("fit_histogram", [True, False])
@pytest.mark.parametrize("lower_half_only", [True, False])
def test_clip_to_background(log_watcher, fit_histogram, lower_half_only):
    """Integrity checks for clipping a simple data array."""
    # Make an array with normal noise
    shape = (100, 100)
    rng = np.random.default_rng(seed=123)
    image = rng.normal(0, 0.01, size=shape)
    mask = np.full(shape, True)

    # Add a strong high outlier, strong low outlier
    image[0, 0] += 1
    image[1, 1] += -1

    # Center found is close to zero, printed when verbose=True
    watcher = log_watcher("jwst.clean_flicker_noise.clean_flicker_noise", message="center: 0.000")
    cfn.clip_to_background(
        image, mask, fit_histogram=fit_histogram, lower_half_only=lower_half_only, verbose=True
    )
    watcher.assert_seen()

    # All outliers clipped with defaults, as well as a small
    # percent of the remaining data
    assert not mask[0, 0]
    assert not mask[1, 1]
    assert np.allclose(np.sum(mask) / mask.size, 0.98, atol=0.019)

    # Upper outlier stays with large sigma_upper
    mask = np.full(shape, True)
    cfn.clip_to_background(
        image, mask, fit_histogram=fit_histogram, lower_half_only=lower_half_only, sigma_upper=1000
    )
    assert mask[0, 0]
    assert not mask[1, 1]
    assert np.allclose(np.sum(mask) / mask.size, 0.98, atol=0.019)

    # Lower outlier stays with large sigma_lower
    mask = np.full(shape, True)
    cfn.clip_to_background(
        image, mask, fit_histogram=fit_histogram, lower_half_only=lower_half_only, sigma_lower=1000
    )
    assert not mask[0, 0]
    assert mask[1, 1]
    assert np.allclose(np.sum(mask) / mask.size, 0.98, atol=0.019)


def test_clip_to_background_fit_fails(log_watcher):
    shape = (10, 10)
    watcher = log_watcher("jwst.clean_flicker_noise.clean_flicker_noise")

    # histogram failure: all data NaN
    image = np.full(shape, np.nan)
    mask = np.full(shape, True)
    watcher.message = "Histogram failed"

    cfn.clip_to_background(image, mask, fit_histogram=True, verbose=True)
    assert np.all(mask)
    watcher.assert_seen()

    # if mask is all False, warning is avoided, mask is unchanged
    mask[:] = False
    cfn.clip_to_background(image, mask, fit_histogram=True, verbose=True)
    assert not np.all(mask)

    # fit failure: data is not normal
    watcher.message = "Gaussian fit failed"
    image = np.full(shape, 0.0)
    mask = np.full(shape, True)
    image[5:, 5:] = 0.1
    cfn.clip_to_background(image, mask, fit_histogram=True, verbose=True)
    assert np.all(mask)
    watcher.assert_seen()


@pytest.mark.parametrize("exptype", ["mos", "mos_fs", "ifu"])
def test_create_mask_nirspec(monkeypatch, exptype):
    # monkeypatch local functions for speed and check that they are called
    # (actual behavior tested in separate unit tests)
    mock = MockUpdate()
    if exptype == "mos":
        # NIRSpec MOS data
        monkeypatch.setattr(cfn, "mask_slits", mock)
        rate_model = make_nirspec_mos_model()
    elif exptype == "mos_fs":
        monkeypatch.setattr(cfn, "mask_slits", mock)
        rate_model = make_nirspec_mos_fs_model()

        # also assign a WCS so the FS can be detected
        rate_model = cfn.post_process_rate(rate_model, assign_wcs=True)
    else:
        # NIRSpec IFU data
        monkeypatch.setattr(cfn, "mask_ifu_slices", mock)
        rate_model = make_nirspec_ifu_model((2048, 2038))

    rate_model.data[:] = 1.0
    mask = cfn.create_mask(rate_model)
    assert mask.shape == rate_model.data.shape

    # Nothing to mask in uniform data
    assert np.all(mask)

    # Slit function not called if mask_science_regions not set
    assert not mock.seen

    # Fixed slit region not blocked
    assert np.all(mask[cfn.NRS_FS_REGION[0] : cfn.NRS_FS_REGION[1]])

    # Call again but block science regions
    mask = cfn.create_mask(rate_model, mask_science_regions=True)
    assert mock.seen
    if exptype == "mos_fs":
        # FS region still not blocked - may need correction
        assert np.all(mask[cfn.NRS_FS_REGION[0] : cfn.NRS_FS_REGION[1]])
    else:
        assert not np.any(mask[cfn.NRS_FS_REGION[0] : cfn.NRS_FS_REGION[1]])


def test_create_mask_miri():
    # small MIRI imaging data
    rate_model = make_small_rate_model()
    rate_model.dq[5, 5] = datamodels.dqflags.pixel["DO_NOT_USE"]

    mask = cfn.create_mask(rate_model)
    assert mask.shape == rate_model.data.shape

    # Nothing to mask in uniform data
    assert np.all(mask)

    # Call again but block non-science regions
    mask = cfn.create_mask(rate_model, mask_science_regions=True)
    assert not mask[5, 5]
    assert np.sum(mask) == mask.size - 1

    rate_model.close()


def test_create_mask_from_rateints():
    # small rateints data
    rate_model = make_small_rateints_model()

    # Add an outlier in each integration
    for i in range(rate_model.data.shape[0]):
        rate_model.data[i, i, i] += 100

    mask = cfn.create_mask(rate_model)
    assert mask.shape == rate_model.data.shape
    assert mask.ndim == 3

    # Outlier is masked in each integration
    assert not mask[0, 0, 0]
    assert not mask[1, 1, 1]
    assert not mask[2, 2, 2]
    assert np.sum(mask) == mask.size - 3

    # Call again but make a single mask:
    # bad data mask is or-ed across integrations
    mask = cfn.create_mask(rate_model, single_mask=True)
    assert mask.shape == rate_model.data.shape[1:]
    assert mask.ndim == 2

    assert not mask[0, 0]
    assert not mask[1, 1]
    assert not mask[2, 2]
    assert np.sum(mask) == mask.size - 3

    rate_model.close()


def test_background_level(log_watcher):
    watcher = log_watcher("jwst.clean_flicker_noise.clean_flicker_noise")

    shape = (100, 100)
    image = np.full(shape, 1.0)
    mask = np.full(shape, True)

    # add an outlier to be clipped
    image[50, 50] = 1000

    # no background
    background = cfn.background_level(image, mask, background_method=None)
    assert background == 0.0

    # median method
    background = cfn.background_level(image, mask, background_method="median")
    assert background == 1.0

    # model method
    background = cfn.background_level(
        image, mask, background_method="model", background_box_size=(10, 10)
    )
    assert background.shape == shape
    assert np.all(background == 1.0)

    # model method with mismatched box size:
    # warns, but completes successfully
    watcher.message = "does not divide evenly"
    background = cfn.background_level(
        image, mask, background_method="model", background_box_size=(32, 32)
    )
    assert background.shape == shape
    assert np.all(background == 1.0)
    watcher.assert_seen()

    # model method with None box size: picks the largest even divisor < 32
    watcher.message = "box size [25, 25]"
    background = cfn.background_level(
        image, mask, background_method="model", background_box_size=None
    )
    assert background.shape == shape
    assert np.all(background == 1.0)
    watcher.assert_seen()

    # make image mostly bad, one good region
    image[:] = np.nan
    image[20:25, 20:25] = 1.0

    # background fit fails: falls back on simple median
    watcher.message = "Background fit failed, using median"
    background = cfn.background_level(
        image, mask, background_method="model", background_box_size=(10, 10)
    )
    assert background == 1.0
    watcher.assert_seen()


@pytest.mark.parametrize("array_type,fraction_good", [("full", 0.65), ("subarray", 0.37)])
@pytest.mark.parametrize("detector", ["NRS1", "NRS2"])
def test_fft_clean(array_type, fraction_good, detector):
    if array_type == "full":
        clean_function = cfn.fft_clean_full_frame
    else:
        clean_function = cfn.fft_clean_subarray

    shape = (80, 80)
    mask = np.full(shape, True)

    # zero image should still come out zero
    image = np.full(shape, 0.0)
    cleaned_image = clean_function(image.copy(), mask, detector)
    assert np.allclose(cleaned_image, 0.0)

    # image with regular vertical pattern, centered on 0
    high = np.full((80, 16), 0.1)
    low = np.full((80, 16), -0.1)
    image = np.hstack([high, low, high, low, high])
    cleaned_image = clean_function(image.copy(), mask, detector)

    # results should be mostly close to zero,
    # some artifacts at pattern and edge boundaries,
    # edge boundaries are currently worse for subarray correction
    good_correction = np.abs(cleaned_image) < 0.01
    assert np.allclose(np.sum(good_correction) / cleaned_image.size, fraction_good, atol=0.1)


def test_fft_clean_error(monkeypatch, log_watcher):
    def raise_error(*args, **kwargs):
        raise np.linalg.LinAlgError("Linear algebra error")

    monkeypatch.setattr(cfn.NSClean, "clean", raise_error)

    shape = (10, 10)
    mask = np.full(shape, True)
    image = np.full(shape, 1.0)

    watcher = log_watcher(
        "jwst.clean_flicker_noise.clean_flicker_noise", message="Error cleaning image"
    )
    cleaned_image = cfn.fft_clean_full_frame(image.copy(), mask, "NRS1")
    assert cleaned_image is None
    watcher.assert_seen()


def test_fft_subarray_clean_error(monkeypatch, log_watcher):
    watcher = log_watcher("jwst.clean_flicker_noise.clean_flicker_noise")

    shape = (10, 10)
    image = np.full(shape, 1.0)

    # Mask is all bad: error message, returns None
    mask = np.full(shape, False)
    watcher.message = "No good pixels"
    # RuntimeWarning: Mean of empty slice
    # RuntimeWarning: invalid value encountered in scalar divide
    with pytest.warns(RuntimeWarning):
        cleaned_image = cfn.fft_clean_subarray(image.copy(), mask, "NRS1")
    assert cleaned_image is None
    watcher.assert_seen()

    # Mask is mostly bad: warns but continues
    mask[5, 5] = True
    watcher.message = "Insufficient reference pixels"
    cleaned_image = cfn.fft_clean_subarray(image.copy(), mask, "NRS1", minfrac=0.5)
    assert np.allclose(cleaned_image, image)
    watcher.assert_seen()

    # Trigger a linear algebra error
    # This may occur when the mask is not completely bad,
    # but there is insufficient data in some region.
    def raise_error(*args, **kwargs):
        raise np.linalg.LinAlgError("Linear algebra error")

    monkeypatch.setattr(cfn.NSCleanSubarray, "clean", raise_error)

    watcher.message = "Error cleaning image"
    cleaned_image = cfn.fft_clean_subarray(image.copy(), mask, "NRS1")
    assert cleaned_image is None
    watcher.assert_seen()


def test_median_clean():
    shape = (20, 20)
    mask = np.full(shape, True)

    # zero image should still come out zero
    image = np.full(shape, 0.0)
    cleaned_image = cfn.median_clean(image, mask, 1)
    assert np.allclose(cleaned_image, 0.0)

    # image with regular vertical pattern, centered on 0
    high = np.full((20, 4), 0.1)
    low = np.full((20, 4), -0.1)
    image = np.hstack([high, low, high, low, high])

    # clean along vertical axis -
    # image should be all zero
    cleaned_image = cfn.median_clean(image, mask, 1)
    assert np.allclose(cleaned_image, 0.0)

    # clean along horizontal axis -
    # cleaning removes median value, from high stripes
    cleaned_image = cfn.median_clean(image, mask, 2)
    assert np.allclose(cleaned_image, image - 0.1)


def test_median_clean_by_channel():
    shape = (2048, 2048)
    mask = np.full(shape, True)

    # zero image should still come out zero
    image = np.full(shape, 0.0)
    cleaned_image = cfn.median_clean(image, mask, 1, fit_by_channel=True)
    assert np.allclose(cleaned_image, 0.0)

    # image with regular vertical pattern, centered on 0
    high = np.full((2048, 16), 0.1)
    low = np.full((2048, 16), -0.1)
    image = np.hstack([high, low] * 64)

    # add offset by vertical channel
    offset = np.zeros_like(image)
    offset[:512] += 0.2
    offset[512:1024] += 0.3
    offset[1024:1536] += 0.4
    offset[1536:] += 0.5
    image += offset

    # clean along vertical axis by channel -
    # image should be all zero
    cleaned_image = cfn.median_clean(image, mask, 1, fit_by_channel=True)
    assert np.allclose(cleaned_image, 0.0)

    # clean along vertical axis for the whole image -
    # vertical pattern is removed, offset by channel stays,
    # but median value is subtracted
    cleaned_image = cfn.median_clean(image, mask, 1, fit_by_channel=False)
    assert np.allclose(cleaned_image, offset - 0.35)

    # clean along horizontal axis -
    # cleaning removes offset value, leaves vertical pattern
    cleaned_image = cfn.median_clean(image, mask, 2, fit_by_channel=True)
    assert np.allclose(cleaned_image, image - offset)


@pytest.mark.parametrize("mask_science", [True, False])
def test_do_correction_miri_imaging_ramp(mask_science):
    ramp_model = make_small_ramp_model()
    cleaned, mask, bg, noise, status = cfn.do_correction(
        ramp_model, mask_science_regions=mask_science
    )

    # uniform data, correction has no effect
    assert cleaned.data.shape == ramp_model.shape
    assert np.allclose(cleaned.data, ramp_model.data)
    assert status == "COMPLETE"

    # extra models are not created
    assert mask is None
    assert bg is None
    assert noise is None

    ramp_model.close()
    cleaned.close()


@pytest.mark.parametrize("mask_science", [True, False])
def test_do_correction_nirspec_rate(mask_science):
    rate_model = make_nirspec_fs_model()
    cleaned, mask, bg, noise, status = cfn.do_correction(
        rate_model, mask_science_regions=mask_science
    )

    # uniform data, correction has no effect
    assert cleaned.data.shape == rate_model.data.shape
    assert np.allclose(cleaned.data, rate_model.data)
    assert status == "COMPLETE"

    # extra models are not created
    assert mask is None
    assert bg is None
    assert noise is None

    rate_model.close()
    cleaned.close()


@pytest.mark.parametrize("single_mask", [True, False])
def test_do_correction_rateints(single_mask):
    rate_model = make_small_rateints_model()
    cleaned, mask, bg, noise, status = cfn.do_correction(
        rate_model, single_mask=single_mask, save_mask=True
    )

    # uniform data, correction has no effect
    assert cleaned.data.shape == rate_model.shape
    assert np.allclose(cleaned.data, rate_model.data)
    assert status == "COMPLETE"

    if single_mask:
        assert mask.shape == cleaned.data.shape[-2:]
    else:
        assert mask.shape == cleaned.data.shape

    rate_model.close()
    cleaned.close()
    mask.close()


def test_do_correction_unsupported(log_watcher):
    ramp_model = make_small_ramp_model()
    ramp_model.meta.exposure.type = "MIR_MRS"

    watcher = log_watcher("jwst.clean_flicker_noise.clean_flicker_noise", message="not supported")
    cleaned, _, _, _, status = cfn.do_correction(ramp_model)
    assert cleaned is ramp_model
    assert status == "SKIPPED"
    watcher.assert_seen()

    ramp_model.close()


def test_do_correction_no_fit_by_channel(log_watcher):
    ramp_model = make_small_ramp_model()

    # fit_by_channel is only used for NIR data -
    # log a warning, but step still completes
    watcher = log_watcher(
        "jwst.clean_flicker_noise.clean_flicker_noise",
        message="can only be used for full-frame NIR",
    )
    cleaned, _, _, _, status = cfn.do_correction(ramp_model, fit_by_channel=True)
    assert status == "COMPLETE"
    watcher.assert_seen()

    ramp_model.close()
    cleaned.close()


@pytest.mark.parametrize("exptype", ["MIR_IMAGE", "NRC_IMAGE", "NIS_IMAGE"])
def test_do_correction_fft_not_allowed(log_watcher, exptype):
    ramp_model = make_small_ramp_model()
    ramp_model.meta.exposure.type = exptype

    # not allowed for MIRI, NIRCAM, NIRISS
    watcher = log_watcher(
        "jwst.clean_flicker_noise.clean_flicker_noise",
        message=f"cannot be applied to exp_type {exptype}",
    )
    cleaned, _, _, _, status = cfn.do_correction(ramp_model, fit_method="fft")
    assert cleaned is ramp_model
    assert status == "SKIPPED"
    watcher.assert_seen()

    ramp_model.close()


@pytest.mark.parametrize("none_value", [None, "none", "None"])
def test_do_correction_no_background(none_value):
    ramp_model = make_small_ramp_model()

    # Don't remove background before fitting noise
    cleaned, _, _, _, _ = cfn.do_correction(ramp_model, background_method=none_value)

    # Output data is all zero: uniform level is removed
    assert np.allclose(cleaned.data, 0.0)

    ramp_model.close()
    cleaned.close()


@pytest.mark.parametrize("ndim", [2, 3])
def test_do_correction_user_mask(tmp_path, ndim):
    ramp_model = make_small_ramp_model()

    if ndim == 3:
        mask = np.full((ramp_model.shape[0], ramp_model.shape[2], ramp_model.shape[3]), True)
        mask_model = datamodels.CubeModel(mask)
    else:
        mask = np.full(ramp_model.shape[-2:], True)
        mask_model = datamodels.ImageModel(mask)
    user_mask = str(tmp_path / "mask.fits")
    mask_model.save(user_mask)
    mask_model.close()

    cleaned, output_mask, _, _, _ = cfn.do_correction(
        ramp_model, user_mask=user_mask, save_mask=True
    )

    assert np.all(output_mask.data == mask_model.data)
    assert np.allclose(cleaned.data, ramp_model.data)

    ramp_model.close()
    cleaned.close()
    output_mask.close()


@pytest.mark.parametrize("input_type", ["rate", "rateints", "ramp"])
def test_do_correction_user_mask_mismatch(tmp_path, input_type, log_watcher):
    shape = (3, 5, 20, 20)
    if input_type == "rate":
        model = make_small_rate_model(shape)
    elif input_type == "rateints":
        model = make_small_rateints_model(shape)
    else:
        model = make_small_ramp_model(shape)

    mask = np.full((10, 10), True)
    mask_model = datamodels.ImageModel(mask)
    user_mask = str(tmp_path / "mask.fits")
    mask_model.save(user_mask)
    mask_model.close()

    watcher = log_watcher(
        "jwst.clean_flicker_noise.clean_flicker_noise", message="Mask does not match"
    )
    cleaned, output_mask, _, _, status = cfn.do_correction(
        model, user_mask=user_mask, save_mask=True
    )

    watcher.assert_seen()
    assert status == "SKIPPED"
    assert output_mask is None

    model.close()
    cleaned.close()


def test_do_correction_user_mask_mismatch_integ(tmp_path, log_watcher):
    shape = (3, 5, 20, 20)
    model = make_small_rateints_model(shape)

    mask = np.full((2, 20, 20), True)
    mask_model = datamodels.CubeModel(mask)
    user_mask = str(tmp_path / "mask.fits")
    mask_model.save(user_mask)
    mask_model.close()

    watcher = log_watcher(
        "jwst.clean_flicker_noise.clean_flicker_noise", message="Mask does not match"
    )
    cleaned, output_mask, _, _, status = cfn.do_correction(
        model, user_mask=user_mask, save_mask=True
    )

    watcher.assert_seen()
    assert status == "SKIPPED"
    assert output_mask is None

    model.close()
    cleaned.close()


@pytest.mark.parametrize("subarray", ["SUBS200A1", "ALLSLITS"])
def test_do_correction_fft_subarray(subarray):
    ramp_model = make_small_ramp_model()
    ramp_model.meta.exposure.type = "NRS_FIXEDSLIT"
    ramp_model.meta.subarray.slowaxis = 1
    ramp_model.meta.subarray.name = subarray

    cleaned, _, _, noise, _ = cfn.do_correction(ramp_model, fit_method="fft", save_noise=True)

    # uniform data, correction has no effect
    assert cleaned.data.shape == ramp_model.data.shape
    assert np.allclose(cleaned.data, ramp_model.data)
    assert np.allclose(noise.data, 0.0)

    ramp_model.close()
    cleaned.close()
    noise.close()


def test_do_correction_fft_full():
    rate_model = make_nirspec_fs_model()

    cleaned, _, _, noise, _ = cfn.do_correction(rate_model, fit_method="fft", save_noise=True)

    # uniform data, correction has no effect
    assert cleaned.data.shape == rate_model.data.shape
    assert np.allclose(cleaned.data, rate_model.data)
    assert np.allclose(noise.data, 0.0)

    rate_model.close()
    cleaned.close()
    noise.close()


def test_do_correction_clean_fails(monkeypatch, log_watcher):
    ramp_model = make_small_ramp_model()
    ramp_model.meta.exposure.type = "NRS_FIXEDSLIT"
    ramp_model.meta.subarray.slowaxis = 1

    # Add some noise so that input and output are
    # expected to be different
    rng = np.random.default_rng(seed=123)
    ramp_model.data += rng.normal(0, 0.1, size=ramp_model.data.shape)
    cleaned, _, _, _, status = cfn.do_correction(ramp_model, fit_method="fft")
    assert not np.allclose(cleaned.data, ramp_model.data)
    assert status == "COMPLETE"

    # Mock a None-value returned from cleaning function
    monkeypatch.setattr(cfn, "fft_clean_subarray", lambda *args, **kwargs: None)

    # Call again
    watcher = log_watcher("jwst.clean_flicker_noise.clean_flicker_noise", message="Cleaning failed")
    cleaned, _, _, _, status = cfn.do_correction(ramp_model, fit_method="fft")

    # Error message issued, status is skipped,
    # output data is the same as input
    watcher.assert_seen()
    assert status == "SKIPPED"
    assert np.allclose(cleaned.data, ramp_model.data)

    ramp_model.close()
    cleaned.close()


@pytest.mark.parametrize("save_type", ["noise", "background"])
@pytest.mark.parametrize("input_type", ["rate", "rateints", "ramp"])
def test_do_correction_save_intermediate(save_type, input_type):
    shape = (3, 5, 20, 20)
    if input_type == "rate":
        model = make_small_rate_model(shape)
    elif input_type == "rateints":
        model = make_small_rateints_model(shape)
    else:
        model = make_small_ramp_model(shape)

    save_bg = save_type == "background"
    save_noise = save_type == "noise"
    cleaned, _, background, noise, status = cfn.do_correction(
        model, save_background=save_bg, save_noise=save_noise
    )

    # Output background model always matches input shape
    # and datamodel type
    if save_bg:
        assert background.data.shape == model.data.shape
        assert type(background) is type(model)
        background.close()
    else:
        assert background is None
    if save_noise:
        assert noise.data.shape == model.data.shape
        assert type(noise) is type(model)
        noise.close()
    else:
        assert noise is None

    model.close()


@pytest.mark.parametrize("input_type", ["rate", "rateints", "ramp"])
def test_do_correction_with_flat_unity(tmp_path, input_type, log_watcher):
    # make input data
    shape = (3, 5, 20, 20)
    if input_type == "rate":
        model = make_small_rate_model(shape)
    elif input_type == "rateints":
        model = make_small_rateints_model(shape)
    else:
        model = make_small_ramp_model(shape)

    # make a flat image matching the input data
    flat = make_flat_model(model, shape=shape[-2:], value=1.0)
    flat_file = str(tmp_path / "flat.fits")
    flat.save(flat_file)

    watcher = log_watcher(
        "jwst.clean_flicker_noise.clean_flicker_noise", message="Dividing by flat"
    )
    cleaned, _, _, _, status = cfn.do_correction(
        model, flat_filename=flat_file, background_method=None
    )
    watcher.assert_seen()
    assert status == "COMPLETE"

    # output is flat with uniform flat, background is perfectly removed
    assert np.all(cleaned.data == 0.0)

    model.close()
    flat.close()


@pytest.mark.parametrize("apply_flat", [True, False])
@pytest.mark.parametrize("input_type", ["rate", "rateints", "ramp"])
def test_do_correction_with_flat_structure(tmp_path, log_watcher, input_type, apply_flat):
    # make input data
    shape = (3, 5, 20, 20)
    if input_type == "rate":
        model = make_small_rate_model(shape)
    elif input_type == "rateints":
        model = make_small_rateints_model(shape)
    else:
        model = make_small_ramp_model(shape)

    # make a flat image matching the input data
    flat = make_flat_model(model, shape=shape[-2:])
    if apply_flat:
        flat_file = str(tmp_path / "flat.fits")
        flat.save(flat_file)
    else:
        flat_file = None

    # multiply the data by the flat to mock real structure
    model.data *= flat.data

    watcher = log_watcher(
        "jwst.clean_flicker_noise.clean_flicker_noise", message="Dividing by flat"
    )
    cleaned, _, _, _, status = cfn.do_correction(model, flat_filename=flat_file)
    assert status == "COMPLETE"

    if apply_flat:
        watcher.assert_seen()

        # output is the same as input: flat structure is not removed
        assert np.all(cleaned.data == model.data)
    else:
        watcher.assert_not_seen()

        # output is not the same as input: flat structure is fit as background/noise
        assert not np.all(cleaned.data == model.data)

    model.close()
    flat.close()


def test_do_correction_with_flat_subarray(tmp_path, log_watcher):
    # make input data
    shape = (3, 5, 20, 20)
    model = make_small_rate_model(shape)

    # make a flat image larger than the input data
    flat_shape = (50, 50)
    flat = make_flat_model(model, shape=flat_shape)
    flat_file = str(tmp_path / "flat.fits")
    flat.save(flat_file)

    # multiply the data by the flat to mock real structure
    model.data *= flat.data[:20, :20]

    watcher = log_watcher(
        "jwst.clean_flicker_noise.clean_flicker_noise", message="Extracting matching subarray"
    )
    cleaned, _, _, _, status = cfn.do_correction(model, flat_filename=flat_file)
    assert status == "COMPLETE"
    watcher.assert_seen()

    # output is the same as input: flat structure is not removed by the
    # cleaning process
    assert np.all(cleaned.data == model.data)

    model.close()
    flat.close()
