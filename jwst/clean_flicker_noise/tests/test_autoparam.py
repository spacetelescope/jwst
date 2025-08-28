import logging

import numpy as np
import pytest

from jwst.clean_flicker_noise import autoparam
from jwst.clean_flicker_noise.tests import helpers
from jwst.lib.basic_utils import LoggingContext


def add_small_sources(image, nsource=1000, value=10.0):
    rng = np.random.default_rng(seed=42)
    ny, nx = image.shape
    y_idx = rng.choice(ny, nsource)
    x_idx = rng.choice(nx, nsource)
    image[y_idx, x_idx] += value


def add_large_source(image, center=None, radius=3.0, value=10.0):
    ny, nx = image.shape
    if center is None:
        cy, cx = ny // 2, nx // 2
    else:
        cy, cx = center
    y, x = np.mgrid[:ny, :nx]
    source = (y - cy) ** 2 + (x - cx) ** 2 < radius**2
    image[source] += value


@pytest.mark.parametrize("instrument", ["NIS", "NRC"])
def test_image_high_mask_high_slope(tmp_path, caplog, instrument):
    # Image is all 1.0
    if instrument == "NIS":
        image = helpers.make_niriss_rate_model()
    else:
        image = helpers.make_nircam_rate_model()

    # Flat has a strong gradient: a large fraction of the image will be masked
    flat = helpers.make_flat_model(image)
    flat_file = str(tmp_path / "flat.fits")
    flat.save(flat_file)

    with LoggingContext(autoparam.log, level=logging.DEBUG):
        if instrument == "NIS":
            param = autoparam.niriss_image_parameters(image, flat_file)
        else:
            param = autoparam.nircam_image_parameters(image, flat_file)
    assert "High mask, low offset, high slope" in caplog.text

    expected = {"apply_flat_field": True, "background_method": "model", "fit_by_channel": False}
    assert param == expected


@pytest.mark.parametrize("instrument", ["NIS", "NRC"])
def test_image_high_mask_low_slope(tmp_path, caplog, instrument):
    # Image is all 1.0
    if instrument == "NIS":
        image = helpers.make_niriss_rate_model()
    else:
        image = helpers.make_nircam_rate_model()

    # Add a large source at the center of the array
    add_large_source(image.data)

    # Flat is all 1.0
    flat = helpers.make_flat_model(image, value=1.0)
    flat_file = str(tmp_path / "flat.fits")
    flat.save(flat_file)

    with LoggingContext(autoparam.log, level=logging.DEBUG):
        if instrument == "NIS":
            param = autoparam.niriss_image_parameters(image, flat_file)
        else:
            param = autoparam.nircam_image_parameters(image, flat_file)
    assert "High mask, low offset, low slope" in caplog.text

    expected = {
        "apply_flat_field": True,
        "background_method": "median",
        "fit_by_channel": False,
    }
    assert param == expected


@pytest.mark.parametrize("instrument", ["NIS", "NRC"])
def test_image_low_mask(tmp_path, caplog, instrument):
    # Image is all 1.0
    if instrument == "NIS":
        image = helpers.make_niriss_rate_model()
    else:
        image = helpers.make_nircam_rate_model()

    # Flat is flat: a small fraction of the image will be masked
    flat = helpers.make_flat_model(image, value=1.0)
    flat_file = str(tmp_path / "flat.fits")
    flat.save(flat_file)

    with LoggingContext(autoparam.log, level=logging.DEBUG):
        if instrument == "NIS":
            param = autoparam.niriss_image_parameters(image, flat_file)
        else:
            param = autoparam.nircam_image_parameters(image, flat_file)
    assert "Low mask, low offset, low slope" in caplog.text

    expected = {
        "apply_flat_field": True,
        "background_method": "median",
        "fit_by_channel": False,
    }
    assert param == expected


@pytest.mark.parametrize("model", ["rate", "rateints", "ramp"])
def test_image_non_rate_models(tmp_path, caplog, model):
    if model == "ramp":
        image = helpers.make_small_ramp_model()
    elif model == "rateints":
        image = helpers.make_small_rateints_model()
    else:
        image = helpers.make_small_rate_model()

    with LoggingContext(autoparam.log, level=logging.DEBUG):
        param = autoparam.nircam_image_parameters(image, "N/A")
    assert "Low mask, low offset, low slope" in caplog.text

    expected = {"apply_flat_field": False, "background_method": "median", "fit_by_channel": False}
    assert param == expected


@pytest.mark.parametrize("flat", ["N/A", None])
def test_niriss_no_flat(flat):
    image = helpers.make_niriss_rate_model()
    param = autoparam.niriss_image_parameters(image, flat)
    assert param["apply_flat_field"] is False


def test_niriss_image_high_offset(tmp_path, caplog):
    # Image is all full frame
    image = helpers.make_niriss_rate_model(shape=(3, 5, 2048, 2048))

    # Add a channel offset
    image.data[:512] += 1.0
    image.data[512:1024] += 10.0
    image.data[1024:1536] += 100.0
    image.data[1536:] += 100.0

    # Flat is all 1.0
    flat = helpers.make_flat_model(image, shape=(2048, 2048), value=1.0)
    flat_file = str(tmp_path / "flat.fits")
    flat.save(flat_file)

    expected = {"apply_flat_field": True, "background_method": "median", "fit_by_channel": True}

    with LoggingContext(autoparam.log, level=logging.DEBUG):
        param = autoparam.niriss_image_parameters(image, flat_file)
    assert "Low mask, high offset" in caplog.text

    assert param == expected


def test_nircam_image_high_offset(tmp_path, caplog):
    # Image is all full frame
    image = helpers.make_nircam_rate_model(shape=(3, 5, 2048, 2048))

    # Add a channel offset
    image.data[:, :512] += 1.0
    image.data[:, 512:1024] += 10.0
    image.data[:, 1024:1536] += 100.0
    image.data[:, 1536:] += 100.0

    with LoggingContext(autoparam.log, level=logging.DEBUG):
        param = autoparam.nircam_image_parameters(image, "N/A")
    assert "Low mask, high offset" in caplog.text

    expected = {"apply_flat_field": False, "background_method": "median", "fit_by_channel": True}
    assert param == expected


def test_niriss_image_high_mask_high_offset_high_slope(tmp_path, caplog):
    # Image is all full frame
    image = helpers.make_niriss_rate_model(shape=(3, 5, 2048, 2048))

    # Add a channel offset
    image.data[:512] += 1.0
    image.data[512:1024] += 10.0
    image.data[1024:1536] += 100.0
    image.data[1536:] += 100.0

    # Add a large source
    add_large_source(image.data, radius=500.0, value=10000.0)

    # Flat has a strong gradient
    flat = helpers.make_flat_model(image, shape=(2048, 2048))
    flat.data /= np.nanmedian(flat.data)
    flat_file = str(tmp_path / "flat.fits")
    flat.save(flat_file)

    expected = {"apply_flat_field": True, "background_method": "model", "fit_by_channel": False}

    with LoggingContext(autoparam.log, level=logging.DEBUG):
        param = autoparam.niriss_image_parameters(image, flat_file)
    assert "High mask, high offset, high slope" in caplog.text

    assert param == expected


def test_niriss_image_high_mask_high_offset_low_slope(tmp_path, caplog):
    # Image is all full frame
    image = helpers.make_niriss_rate_model(shape=(3, 5, 2048, 2048))

    # Add a channel offset
    image.data[:512] += 10.0
    image.data[512:1024] += 10.0
    image.data[1024:1536] += 10.08
    image.data[1536:] += 10.08

    # Add a large source
    add_large_source(image.data, radius=500.0, value=10000.0, center=(520, 520))

    # Flat is flat
    flat = helpers.make_flat_model(image, shape=(2048, 2048), value=1.0)
    flat_file = str(tmp_path / "flat.fits")
    flat.save(flat_file)

    expected = {"apply_flat_field": True, "background_method": "median", "fit_by_channel": False}

    with LoggingContext(autoparam.log, level=logging.DEBUG):
        param = autoparam.niriss_image_parameters(image, flat_file)
    assert "High mask, high offset, low slope" in caplog.text

    assert param == expected


def test_niriss_image_low_mask_high_slope(tmp_path, caplog):
    # Image is all 1.0
    image = helpers.make_niriss_rate_model(shape=(3, 5, 2048, 2048))

    # Add small sources to flag
    add_small_sources(image.data, nsource=2000)

    # Add diffuse sources to make the slope high
    add_large_source(image.data, radius=700, value=0.2, center=(600, 600))
    add_large_source(image.data, radius=700, value=-0.2, center=(600, 1480))
    add_large_source(image.data, radius=700, value=0.2, center=(1480, 600))
    add_large_source(image.data, radius=700, value=-0.2, center=(1480, 1480))

    # Make channel 1 noisy
    image.data[1537:] += 0.1 * np.sin(np.arange(2048))

    # Flat is flat
    flat = helpers.make_flat_model(image, shape=(2048, 2048), value=1.0)
    flat_file = str(tmp_path / "flat.fits")
    flat.save(flat_file)

    expected = {"apply_flat_field": True, "background_method": "model", "fit_by_channel": True}

    with LoggingContext(autoparam.log, level=logging.DEBUG):
        param = autoparam.niriss_image_parameters(image, flat_file)
    assert "Low mask, low offset, high slope" in caplog.text
    assert "High channel sigma" in caplog.text

    assert param == expected


def test_niriss_image_low_mask_medium_slope_high_ratio(tmp_path, caplog):
    # Image is all 1.0
    image = helpers.make_niriss_rate_model(shape=(3, 5, 100, 100))

    # Add small sources to flag
    add_small_sources(image.data, nsource=100)

    # Add an offset to the first half of the data
    image.data[:50] += 0.001

    # Flat is flat
    flat = helpers.make_flat_model(image, shape=(100, 100), value=1.0)
    flat_file = str(tmp_path / "flat.fits")
    flat.save(flat_file)

    expected = {"apply_flat_field": True, "background_method": "median", "fit_by_channel": False}

    with LoggingContext(autoparam.log, level=logging.DEBUG):
        param = autoparam.niriss_image_parameters(image, flat_file)
    assert "Low mask, low offset, medium slope, high slope ratio" in caplog.text
    assert "Low channel sigma" in caplog.text

    assert param == expected


def test_niriss_image_low_mask_medium_slope_low_ratio(tmp_path, caplog):
    # Image is all 1.0
    image = helpers.make_niriss_rate_model(shape=(3, 5, 100, 100))

    # Add small sources to flag
    add_small_sources(image.data, nsource=100)

    # Add a couple large low-valued regions near the center
    add_large_source(image.data, radius=40, value=0.01, center=(47, 47))
    add_large_source(image.data, radius=40, value=0.01, center=(53, 53))

    # Flat is flat
    flat = helpers.make_flat_model(image, shape=(100, 100), value=1.0)
    flat_file = str(tmp_path / "flat.fits")
    flat.save(flat_file)

    expected = {"apply_flat_field": True, "background_method": "model", "fit_by_channel": False}

    with LoggingContext(autoparam.log, level=logging.DEBUG):
        param = autoparam.niriss_image_parameters(image, flat_file)
    assert "Low mask, low offset, medium slope, low slope ratio" in caplog.text
    assert "Low channel sigma" in caplog.text

    assert param == expected


def test_fit_line():
    data = np.arange(10)
    intercept, slope, sigma = autoparam._line_fit(data)
    assert np.isclose(slope, 1)
    assert np.isclose(intercept, 0)
    assert np.isclose(sigma, 0)


@pytest.mark.parametrize("intercept_only", [True, False])
def test_fit_flat_line(monkeypatch, intercept_only):
    # Sometimes, depending on architecture, fitting a flat line returns
    # a single coefficient only:
    #     https://github.com/numpy/numpy/issues/26617
    # Mock that case here to make sure it is covered.
    if intercept_only:
        monkeypatch.setattr(
            np.polynomial.Polynomial, "convert", lambda *args: np.polynomial.Polynomial(coef=[1.0])
        )

    data = np.ones(10)
    intercept, slope, sigma = autoparam._line_fit(data)
    assert np.isclose(slope, 0)
    assert np.isclose(intercept, 1)
    assert np.isclose(sigma, 0)
