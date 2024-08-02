import logging

import gwcs
import numpy as np
import pytest
from stdatamodels.jwst import datamodels

from jwst.assign_wcs.tests.test_nirspec import (
    create_nirspec_ifu_file, create_nirspec_fs_file)
from jwst.msaflagopen.tests.test_msa_open import make_nirspec_mos_model, get_file_path
from jwst.clean_flicker_noise import clean_flicker_noise as cfn
from jwst.tests.helpers import LogWatcher


def add_metadata(model, shape):
    model.meta.instrument.name = 'MIRI'
    model.meta.instrument.detector = 'MIRIMAGE'
    model.meta.instrument.filter = 'F480M'
    model.meta.observation.date = '2015-10-13'
    model.meta.observation.time = "00:00:00"
    model.meta.exposure.type = 'MIR_IMAGE'
    model.meta.exposure.group_time = 1.0
    model.meta.subarray.name = 'FULL'
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
    model.meta.exposure.readpatt = 'FASTR1'


def make_small_ramp_model(shape):
    rampmodel = datamodels.RampModel(shape)
    add_metadata(rampmodel, shape)

    # Make data with a constant rate
    for group in range(shape[1]):
        rampmodel.data[:, group, :, :] = group

    return rampmodel


def make_small_rate_model(shape):
    ratemodel = datamodels.ImageModel(shape[2:])
    add_metadata(ratemodel, shape)
    ratemodel.data[:] = 1.0
    return ratemodel


def make_small_rateints_model(shape):
    ratemodel = datamodels.CubeModel((shape[0], shape[2], shape[3]))
    add_metadata(ratemodel, shape)
    ratemodel.data[:] = 1.0
    return ratemodel


def make_nirspec_ifu_model(shape):
    hdul = create_nirspec_ifu_file(grating='PRISM', filter='CLEAR',
                                   gwa_xtil=0.35986012, gwa_ytil=0.13448857,
                                   gwa_tilt=37.1)
    hdul['SCI'].data = np.ones(shape, dtype=float)
    rate_model = datamodels.IFUImageModel(hdul)
    hdul.close()
    return rate_model


def make_nirspec_mos_fs_model():
    mos_model = make_nirspec_mos_model()
    mos_model.meta.instrument.msa_metadata_file = get_file_path(
        'msa_fs_configuration.fits')
    return mos_model


def make_nirspec_fs_model():
    hdul = create_nirspec_fs_file(grating="G140M", filter="F100LP")
    hdul['SCI'].data = np.ones((2048, 2048), dtype=float)
    rate_model = datamodels.ImageModel(hdul)
    hdul.close()
    return rate_model


@pytest.fixture
def log_watcher(monkeypatch):
    # Set a log watcher to check for a log message at any level
    # in the clean_flicker_noise module
    watcher = LogWatcher('')
    logger = logging.getLogger('jwst.clean_flicker_noise.clean_flicker_noise')
    for level in ['debug', 'info', 'warning', 'error']:
        monkeypatch.setattr(logger, level, watcher)
    return watcher


class MockUpdate:
    def __init__(self):
        self.seen = False

    def __call__(self, input_model, mask):
        self.seen = True
        return mask


def test_make_rate(log_watcher):
    shape = (3, 5, 10, 10)
    ramp_model = make_small_ramp_model(shape)

    log_watcher.message = 'Creating draft rate file'
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
    log_watcher.assert_seen()


def test_postprocess_rate_nirspec(log_watcher):
    rate_model = make_nirspec_mos_model()

    log_watcher.message = 'Assigning a WCS'
    result = cfn.post_process_rate(rate_model, assign_wcs=True)
    assert isinstance(result.meta.wcs, gwcs.WCS)
    assert np.sum(result.dq & datamodels.dqflags.pixel['MSA_FAILED_OPEN']) == 0
    log_watcher.assert_seen()

    log_watcher.message = 'Flagging failed-open'
    result = cfn.post_process_rate(result, msaflagopen=True)
    assert np.sum(result.dq & datamodels.dqflags.pixel['MSA_FAILED_OPEN']) > 0
    log_watcher.assert_seen()

    rate_model.close()
    result.close()


def test_postprocess_rate_miri(log_watcher):
    shape = (3, 5, 10, 10)
    rate_model = make_small_rate_model(shape)
    assert np.sum(rate_model.dq & datamodels.dqflags.pixel['NON_SCIENCE']) == 0

    log_watcher.message = 'Retrieving flat DQ'
    result = cfn.post_process_rate(rate_model, flat_dq=True)
    log_watcher.assert_seen()
    assert np.sum(result.dq & datamodels.dqflags.pixel['NON_SCIENCE']) > 0
    assert np.all(result.data == rate_model.data)

    rate_model.close()
    result.close()


@pytest.mark.slow
def test_mask_ifu_slices():
    shape = (2048, 2048)
    rate_model = make_nirspec_ifu_model(shape)

    rate_model = cfn.post_process_rate(rate_model, assign_wcs=True)
    mask = np.full_like(rate_model.data, True)

    # Mark IFU science regions as False: about 10% of the array
    cfn.mask_ifu_slices(rate_model, mask)
    assert np.allclose(np.sum(mask), 0.9 * mask.size, atol=0.01 * mask.size)


@pytest.mark.parametrize('exptype,blocked',
                         [('mos', .012), ('mos_fs', .025), ('fs', .05)])
def test_mask_slits(exptype, blocked):
    if exptype == 'mos':
        rate_model = make_nirspec_mos_model()
    elif exptype == 'mos_fs':
        rate_model = make_nirspec_mos_fs_model()
    else:
        rate_model = make_nirspec_fs_model()

    # Assign a WCS
    rate_model = cfn.post_process_rate(rate_model, assign_wcs=True)

    # Mark slits as False
    mask = np.full_like(rate_model.data, True)
    cfn.mask_slits(rate_model, mask)

    # Check that the fraction of the array blocked is as expected
    assert np.allclose(np.sum(mask) / mask.size, 1 - blocked, atol=0.001)


@pytest.mark.parametrize('fit_histogram', [True, False])
@pytest.mark.parametrize('lower_half_only', [True, False])
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
    log_watcher.message = "center: 0.000"
    cfn.clip_to_background(
        image, mask, fit_histogram=fit_histogram,
        lower_half_only=lower_half_only, verbose=True)
    log_watcher.assert_seen()

    # All outliers clipped with defaults, as well as a small
    # percent of the remaining data
    assert not mask[0, 0]
    assert not mask[1, 1]
    assert np.allclose(np.sum(mask) / mask.size, 0.98, atol=0.019)

    # Upper outlier stays with large sigma_upper
    mask = np.full(shape, True)
    cfn.clip_to_background(
        image, mask, fit_histogram=fit_histogram,
        lower_half_only=lower_half_only, sigma_upper=1000)
    assert mask[0, 0]
    assert not mask[1, 1]
    assert np.allclose(np.sum(mask) / mask.size, 0.98, atol=0.019)

    # Lower outlier stays with large sigma_lower
    mask = np.full(shape, True)
    cfn.clip_to_background(
        image, mask, fit_histogram=fit_histogram,
        lower_half_only=lower_half_only, sigma_lower=1000)
    assert not mask[0, 0]
    assert mask[1, 1]
    assert np.allclose(np.sum(mask) / mask.size, 0.98, atol=0.019)


def test_clip_to_background_fit_fails(log_watcher):
    shape = (10, 10)

    # histogram failure: all data NaN
    log_watcher.message = "Histogram failed"
    image = np.full(shape, np.nan)
    mask = np.full(shape, True)
    with pytest.warns(RuntimeWarning):
        cfn.clip_to_background(image, mask, fit_histogram=True, verbose=True)
    assert np.all(mask)
    log_watcher.assert_seen()

    # if mask is all False, warning is avoided, mask is unchanged
    mask[:] = False
    cfn.clip_to_background(image, mask, fit_histogram=True, verbose=True)
    assert not np.all(mask)

    # fit failure: data is not normal
    log_watcher.message = "Gaussian fit failed"
    image = np.full(shape, 0.0)
    mask = np.full(shape, True)
    image[5:, 5:] = 0.1
    cfn.clip_to_background(image, mask, fit_histogram=True, verbose=True)
    assert np.all(mask)
    log_watcher.assert_seen()


@pytest.mark.parametrize('exptype', ['mos', 'mos_fs', 'ifu'])
def test_create_mask_nirspec(monkeypatch, exptype):
    # monkeypatch local functions for speed and check that they are called
    # (actual behavior tested in separate unit tests)
    mock = MockUpdate()
    if exptype == 'mos':
        # NIRSpec MOS data
        monkeypatch.setattr(cfn, 'mask_slits', mock)
        rate_model = make_nirspec_mos_model()
    elif exptype == 'mos_fs':
        monkeypatch.setattr(cfn, 'mask_slits', mock)
        rate_model = make_nirspec_mos_fs_model()

        # also assign a WCS so the FS can be detected
        rate_model = cfn.post_process_rate(rate_model, assign_wcs=True)
    else:
        # NIRSpec IFU data
        monkeypatch.setattr(cfn, 'mask_ifu_slices', mock)
        rate_model = make_nirspec_ifu_model((2048, 2038))

    rate_model.data[:] = 1.0
    mask = cfn.create_mask(rate_model)
    assert mask.shape == rate_model.data.shape

    # Nothing to mask in uniform data
    assert np.all(mask)

    # Slit function not called if mask_science_regions not set
    assert not mock.seen

    # Fixed slit region not blocked
    assert np.all(mask[cfn.NRS_FS_REGION[0]:cfn.NRS_FS_REGION[1]])

    # Call again but block science regions
    mask = cfn.create_mask(rate_model, mask_science_regions=True)
    assert mock.seen
    if exptype == 'mos_fs':
        # FS region still not blocked - may need correction
        assert np.all(mask[cfn.NRS_FS_REGION[0]:cfn.NRS_FS_REGION[1]])
    else:
        assert not np.any(mask[cfn.NRS_FS_REGION[0]:cfn.NRS_FS_REGION[1]])


def test_create_mask_miri():
    # small MIRI imaging data
    shape = (3, 5, 10, 10)
    rate_model = make_small_rate_model(shape)
    rate_model.dq[5, 5] = datamodels.dqflags.pixel['NON_SCIENCE']

    mask = cfn.create_mask(rate_model)
    assert mask.shape == rate_model.data.shape

    # Nothing to mask in uniform data
    assert np.all(mask)

    # Call again but block non-science regions
    mask = cfn.create_mask(rate_model, mask_science_regions=True)
    assert not mask[5, 5]
    assert np.sum(mask) == mask.size - 1


def test_create_mask_from_rateints():
    # small rateints data
    shape = (3, 5, 10, 10)
    rate_model = make_small_rateints_model(shape)

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


def test_background_level(log_watcher):
    shape = (100, 100)
    image = np.full(shape, 1.0)
    mask = np.full(shape, True)

    # add an outlier to be clipped
    image[50, 50] = 1000

    # no background
    background = cfn.background_level(image, mask, background_method=None)
    assert background == 0.0

    # median method
    background = cfn.background_level(image, mask, background_method='median')
    assert background == 1.0

    # model method
    background = cfn.background_level(
        image, mask, background_method='model', background_box_size=(10, 10))
    assert background.shape == shape
    assert np.all(background == 1.0)

    # model method with mismatched box size:
    # warns, but completes successfully
    log_watcher.message = "does not divide evenly"
    background = cfn.background_level(
        image, mask, background_method='model', background_box_size=None)
    assert background.shape == shape
    assert np.all(background == 1.0)
    log_watcher.assert_seen()

    # make image mostly bad, one good region
    image[:] = np.nan
    image[20:25, 20:25] = 1.0

    # background fit fails: falls back on simple median
    log_watcher.message = "Background fit failed, using median"
    background = cfn.background_level(
        image, mask, background_method='model', background_box_size=(10, 10))
    assert background == 1.0
    log_watcher.assert_seen()


@pytest.mark.parametrize('array_type', ['full', 'subarray'])
@pytest.mark.parametrize('detector', ['NRS1', 'NRS2'])
def test_fft_clean(array_type, detector):
    if array_type == 'full':
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
    # some artifacts at pattern and edge boundaries
    good_correction = np.abs(cleaned_image) < 0.01
    assert np.allclose(np.sum(good_correction) / cleaned_image.size, 0.65, atol=0.1)


@pytest.mark.parametrize('array_type', ['full', 'subarray'])
def test_fft_clean_error(array_type, monkeypatch, log_watcher):
    def raise_error(*args, **kwargs):
        raise np.linalg.LinAlgError('Linear algebra error')

    if array_type == 'full':
        clean_function = cfn.fft_clean_full_frame
        monkeypatch.setattr(cfn.NSClean, 'clean', raise_error)
    else:
        clean_function = cfn.fft_clean_subarray
        monkeypatch.setattr(cfn.NSCleanSubarray, 'clean', raise_error)

    shape = (10, 10)
    mask = np.full(shape, True)
    image = np.full(shape, 1.0)

    log_watcher.message = "Error cleaning image"
    cleaned_image = clean_function(image.copy(), mask, 'NRS1')
    assert cleaned_image is None
    log_watcher.assert_seen()


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
