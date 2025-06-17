import pytest
import numpy as np


from stdatamodels.jwst import datamodels
from jwst.background import BackgroundStep
from jwst.background.background_sub_soss import (
    find_discontinuity,
    generate_background_masks,
    subtract_soss_bkg,
    BACKGROUND_MASK_CUTOFF,
    SUBSTRIP96_ROWSTART,
)

DETECTOR_SHAPE = (256, 2048)
BKG_SCALING = 0.123
INITIAL_NAN_FRACTION = 1e-4
INITIAL_OUTLIER_FRACTION = 1e-3


@pytest.fixture(scope="module")
def generate_background_template():
    template = np.ones(DETECTOR_SHAPE)
    template[:, 710] = 1.6
    template[:, 711:] = 2.0
    return template


@pytest.fixture(scope="module")
def mock_data(generate_background_template):
    """Synthetic data with NaNs, noise, and known background structure."""

    err_scaling = 0.05
    nan_fraction = INITIAL_NAN_FRACTION

    # make random data and error arrays
    rng = np.random.default_rng(seed=42)
    data = rng.normal(0, 1, DETECTOR_SHAPE)
    # ensure all errors are positive and not too close to zero
    err = err_scaling * (1 + rng.normal(0, 1, DETECTOR_SHAPE) ** 2)

    # add NaNs
    num_nans = int(data.size * nan_fraction)
    nan_indices = np.unravel_index(rng.choice(data.size, num_nans), data.shape)
    data[nan_indices] = np.nan
    err[nan_indices] = np.nan

    # add some outliers
    num_outliers = int(data.size * INITIAL_OUTLIER_FRACTION)
    outlier_indices = np.unravel_index(rng.choice(data.size, num_outliers), data.shape)
    data[outlier_indices] = rng.normal(100, 1, num_outliers)

    data[nan_indices] = np.nan
    err[nan_indices] = np.nan

    # also add a small background to the data with same structure
    # as the known reference background to see if it will get removed
    data += generate_background_template * BKG_SCALING

    return data, err


@pytest.fixture(scope="module")
def generate_soss_image(mock_data):
    image = datamodels.ImageModel()
    image.data = mock_data[0]
    image.err = mock_data[1]
    image.dq = np.isnan(image.data)
    return image


@pytest.fixture(scope="module")
def generate_soss_cube(mock_data):
    cube = datamodels.CubeModel()
    cube.data = np.array([mock_data[0]] * 10)
    cube.err = np.array([mock_data[1]] * 10)
    cube.dq = np.isnan(cube.data)
    return cube


@pytest.fixture(scope="module")
def generate_soss_cube_substrip96(mock_data):
    cube = datamodels.CubeModel()
    cube.data = np.array([mock_data[0][SUBSTRIP96_ROWSTART : SUBSTRIP96_ROWSTART + 96, :]] * 10)
    cube.err = np.array([mock_data[1][SUBSTRIP96_ROWSTART : SUBSTRIP96_ROWSTART + 96, :]] * 10)
    cube.dq = np.isnan(cube.data)
    return cube


def test_discontinuity_finder(generate_background_template):
    disc = find_discontinuity(generate_background_template)
    # Shape consistency
    assert disc.shape == (256, 2)
    # Gradient behavior check
    assert all(disc[:, 0] == 710)


@pytest.mark.parametrize("n_repeats", [1, 10])
def test_generate_background_masks(generate_background_template, n_repeats):
    left, right = generate_background_masks(
        generate_background_template, n_repeats, for_fitting=True
    )
    # Test cutoff
    assert np.sum(right[BACKGROUND_MASK_CUTOFF:]) == 0.0
    # Check behavior across integrations
    np.testing.assert_array_equal(left[0], left[-1])
    assert right.shape == (n_repeats, *DETECTOR_SHAPE)


def test_subtract_soss_bkg(
    generate_background_template,
    generate_soss_image,
    generate_soss_cube,
    generate_soss_cube_substrip96,
):
    template_model = datamodels.SossBkgModel(
        data=np.stack((generate_background_template, generate_background_template))
    )

    # Test image-based correction.
    result = subtract_soss_bkg(generate_soss_image, template_model, 35.0, [25.0, 50.0])
    assert type(result) == type(generate_soss_image)

    # Test cube-based correction.
    result = subtract_soss_bkg(generate_soss_cube, template_model, 35.0, [25.0, 50.0])
    assert type(result) == type(generate_soss_cube)

    mock_model = generate_soss_cube_substrip96
    mock_model.meta.instrument.filter = "CLEAR"
    mock_model.meta.instrument.pupil = "GR700XD"
    mock_model.meta.exposure.type = "NIS_SOSS"

    # Test step-level call along with substrip96-shaped data.
    result = BackgroundStep.call(
        mock_model,
        override_bkg=template_model,
    )

    assert type(result) == type(generate_soss_cube_substrip96)
