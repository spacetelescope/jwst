import numpy as np
import pytest

from jwst.clean_flicker_noise.background_level import background_level, clip_to_background


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
    watcher = log_watcher("jwst.clean_flicker_noise.background_level", message="center: 0.000")
    clip_to_background(
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
    clip_to_background(
        image, mask, fit_histogram=fit_histogram, lower_half_only=lower_half_only, sigma_upper=1000
    )
    assert mask[0, 0]
    assert not mask[1, 1]
    assert np.allclose(np.sum(mask) / mask.size, 0.98, atol=0.019)

    # Lower outlier stays with large sigma_lower
    mask = np.full(shape, True)
    clip_to_background(
        image, mask, fit_histogram=fit_histogram, lower_half_only=lower_half_only, sigma_lower=1000
    )
    assert not mask[0, 0]
    assert mask[1, 1]
    assert np.allclose(np.sum(mask) / mask.size, 0.98, atol=0.019)


def test_clip_to_background_fit_fails(log_watcher):
    shape = (10, 10)
    watcher = log_watcher("jwst.clean_flicker_noise.background_level")

    # histogram failure: all data NaN
    image = np.full(shape, np.nan)
    mask = np.full(shape, True)
    watcher.message = "Histogram failed"

    clip_to_background(image, mask, fit_histogram=True, verbose=True)
    assert np.all(mask)
    watcher.assert_seen()

    # if mask is all False, warning is avoided, mask is unchanged
    mask[:] = False
    clip_to_background(image, mask, fit_histogram=True, verbose=True)
    assert not np.all(mask)

    # fit failure: data is not normal
    watcher.message = "Gaussian fit failed"
    image = np.full(shape, 0.0)
    mask = np.full(shape, True)
    image[5:, 5:] = 0.1
    clip_to_background(image, mask, fit_histogram=True, verbose=True)
    assert np.all(mask)
    watcher.assert_seen()


def test_background_level(log_watcher):
    watcher = log_watcher("jwst.clean_flicker_noise.background_level")

    shape = (100, 100)
    image = np.full(shape, 1.0)
    mask = np.full(shape, True)

    # add an outlier to be clipped
    image[50, 50] = 1000

    # no background
    background = background_level(image, mask, background_method=None)
    assert background == 0.0

    # median method
    background = background_level(image, mask, background_method="median")
    assert background == 1.0

    # model method
    background = background_level(
        image, mask, background_method="model", background_box_size=(10, 10)
    )
    assert background.shape == shape
    assert np.all(background == 1.0)

    # model method with mismatched box size:
    # warns, but completes successfully
    watcher.message = "does not divide evenly"
    background = background_level(
        image, mask, background_method="model", background_box_size=(32, 32)
    )
    assert background.shape == shape
    assert np.all(background == 1.0)
    watcher.assert_seen()

    # model method with None box size: picks the largest even divisor < 32
    watcher.message = "box size [25, 25]"
    background = background_level(image, mask, background_method="model", background_box_size=None)
    assert background.shape == shape
    assert np.all(background == 1.0)
    watcher.assert_seen()

    # make image mostly bad, one good region
    image[:] = np.nan
    image[20:25, 20:25] = 1.0

    # background fit fails: falls back on simple median
    watcher.message = "Background fit failed, using median"
    background = background_level(
        image, mask, background_method="model", background_box_size=(10, 10)
    )
    assert background == 1.0
    watcher.assert_seen()
