import json

import numpy as np
import pytest
import stdatamodels.jwst.datamodels as dm
from astropy.modeling import models

from jwst.pathloss.pathloss import calculate_pathloss_vector
from jwst.ta_center.ta_center_step import (
    JWST_DIAMETER,
    PIXSCALE,
    TACenterStep,
    _get_wavelength,
)

X_REF_SLIT = 326.13
Y_REF_SLIT = 300.7
X_REF_SLITLESS = 38.5
Y_REF_SLITLESS = 829.0
MIRI_DETECTOR_SHAPE = (1024, 1032)  # (ny, nx) for MIRI imager


@pytest.fixture
def mock_specwcs_model(tmp_path):
    """
    Create a mock MIRI LRS specwcs reference file.

    Only contains the metadata needed for TA centering step, not a full specwcs model.

    Returns
    -------
    str
        Path to the saved mock specwcs reference file.
    """
    specwcs_model = dm.MiriLRSSpecwcsModel()

    # Set the reference positions for slit and slitless modes
    specwcs_model.meta.x_ref = X_REF_SLIT
    specwcs_model.meta.y_ref = Y_REF_SLIT
    specwcs_model.meta.x_ref_slitless = X_REF_SLITLESS
    specwcs_model.meta.y_ref_slitless = Y_REF_SLITLESS

    # Save the model to a file
    specwcs_filepath = tmp_path / "mock_specwcs.fits"
    specwcs_model.save(specwcs_filepath)

    return str(specwcs_filepath)


@pytest.fixture
def mock_pathloss_model(tmp_path):
    """
    Create a mock MIRI LRS pathloss reference file.

    Returns
    -------
    str
        Path to the file.
    """
    pathloss_model = dm.MirLrsPathlossModel()

    # Define wavelength grid over MIRI LRS wavelength range
    n_wavelengths = 50
    wavelengths = np.linspace(5.0, 14.0, n_wavelengths).astype(np.float32)

    # Define spatial grids - slit is long in x, narrow in y
    n_x = 51
    n_y = 21
    x_grid = np.linspace(-2.5, 2.5, n_x)  # arcsec
    y_grid = np.linspace(-0.5, 0.5, n_y)  # arcsec

    # Define pathloss data. This will be a slight Gaussian drop-off from the center
    # until +/- 5 pixels in y, then it'll go to zero. Flat in x (along the slit).
    pathloss_data = np.zeros((n_wavelengths, n_y, n_x), dtype=np.float32)
    y_offsets = y_grid[:, np.newaxis]  # Shape: (n_y, 1), in arcsec
    pathloss_2d = np.exp(-2 * y_offsets**2)  # Gaussian in y, broadcast to x

    # apply the cutoff in y direction
    y_pixel_scale = (y_grid[-1] - y_grid[0]) / (n_y - 1)
    cutoff_arcsec = 5 * y_pixel_scale  # ~0.25 arcsec
    pathloss_2d[np.abs(y_offsets) > cutoff_arcsec] = 0.0
    pathloss_data[:] = pathloss_2d

    # Create the pathloss table as a structured numpy array
    # The schema requires columns: wavelength (1D), pathloss (2D), pathloss_err (2D)
    dtype = [
        ("wavelength", np.float32),
        ("pathloss", np.float32, (n_y, n_x)),
        ("pathloss_err", np.float32, (n_y, n_x)),
    ]
    pathloss_table = np.zeros(n_wavelengths, dtype=dtype)
    pathloss_table["wavelength"] = wavelengths
    pathloss_table["pathloss"] = pathloss_data
    pathloss_table["pathloss_err"] = np.ones_like(pathloss_data) * 0.01
    pathloss_model.pathloss_table = pathloss_table

    # Set up wcsinfo
    pathloss_model.meta.wcsinfo.crpix1 = (n_x + 1) / 2.0  # Reference pixel in x
    pathloss_model.meta.wcsinfo.crpix2 = (n_y + 1) / 2.0  # Reference pixel in y
    pathloss_model.meta.wcsinfo.crval1 = 0.0  # Reference value in x (arcsec)
    pathloss_model.meta.wcsinfo.crval2 = 0.0  # Reference value in y (arcsec)

    # Pixel scale in arcsec/pixel for the spatial dimensions
    x_pixel_scale = (x_grid[-1] - x_grid[0]) / (n_x - 1)
    y_pixel_scale = (y_grid[-1] - y_grid[0]) / (n_y - 1)
    pathloss_model.meta.wcsinfo.cdelt1 = x_pixel_scale
    pathloss_model.meta.wcsinfo.cdelt2 = y_pixel_scale

    # Save the model to a file
    pathloss_filepath = tmp_path / "mock_pathloss.fits"
    pathloss_model.save(pathloss_filepath)

    return str(pathloss_filepath)


@pytest.fixture
def mock_filteroffset_model(tmp_path):
    """
    Create a mock MIRI filter offset reference file for testing.

    This contains the column and row offsets for the F1500W filter.

    Returns
    -------
    str
        Path to the saved mock filteroffset reference file.
    """
    filteroffset_model = dm.FilteroffsetModel()

    # Set required metadata fields for asdf file to validate
    filteroffset_model.meta.description = ""
    filteroffset_model.meta.reftype = "FILTEROFFSET"
    filteroffset_model.meta.author = ""
    filteroffset_model.meta.pedigree = ""
    filteroffset_model.meta.useafter = "2000-01-01T00:00:00"
    filteroffset_model.meta.instrument.name = "MIRI"
    filteroffset_model.meta.instrument.detector = "MIRIMAGE"

    # Create filter entry for F1500W with no offset
    filter_entry = {"filter": "F1500W", "column_offset": 0.0, "row_offset": 0.0, "pupil": "N/A"}
    filteroffset_model.filters.append(filter_entry)

    # Save to file
    filteroffset_filepath = tmp_path / "mock_filteroffset.asdf"
    filteroffset_model.save(filteroffset_filepath)
    return str(filteroffset_filepath)


@pytest.fixture
def mock_references(monkeypatch, mock_specwcs_model, mock_pathloss_model, mock_filteroffset_model):
    """
    Monkeypatch the get_reference_file method to return mock reference files.

    This allows tests to run without requiring CRDS access, and ensures that
    these tests won't start breaking if the values in the pathloss reference file
    are changed in CRDS.
    """

    def mock_get_reference_file(self, model, reftype):
        """Mock implementation of get_reference_file."""
        if reftype == "specwcs":
            return mock_specwcs_model
        elif reftype == "pathloss":
            return mock_pathloss_model
        elif reftype == "filteroffset":
            return mock_filteroffset_model
        else:
            raise ValueError(f"Unexpected reference type: {reftype}")

    monkeypatch.setattr(TACenterStep, "get_reference_file", mock_get_reference_file)


def make_slitless_data(wavelength=15.0, offset=(0, 0)):
    """
    Make a fake MIRI LRS slitless TAQ verification image.

    Parameters
    ----------
    wavelength : float, optional
        Wavelength in microns.
    offset : tuple, optional
        (x,y) offset of the source from the reference center in pixels.

    Returns
    -------
    data : 2D ndarray
        Simulated slitless TAQ image.
    """
    # Create coordinate grids using full MIRI detector size
    y, x = np.mgrid[0 : MIRI_DETECTOR_SHAPE[0], 0 : MIRI_DETECTOR_SHAPE[1]]

    # Calculate diffraction-limited Airy disk radius
    # First zero of Airy disk: theta = 1.22 * lambda / D (in radians)
    theta_rad = 1.22 * wavelength * 1e-6 / JWST_DIAMETER
    theta_arcsec = theta_rad * 206265  # radians to arcseconds
    radius_pixels = theta_arcsec / PIXSCALE  # to pixels

    # Calculate source center position using slitless reference position
    x_center = X_REF_SLITLESS + offset[0]
    y_center = Y_REF_SLITLESS + offset[1]

    # Simulate source as an Airy disk model at diffraction-limited resolution
    airy = models.AiryDisk2D(amplitude=1000.0, x_0=x_center, y_0=y_center, radius=radius_pixels)
    data = airy(x, y)

    # Add some noise with fixed seed for reproducibility
    rng = np.random.default_rng(42)
    noise = rng.normal(0, 1, data.shape)
    data += noise

    return data


def _make_taq_model(data):
    model = dm.ImageModel()
    model.data = data
    model.meta.instrument.filter = "F1500W"

    # Add subarray metadata (use FULL for now)
    model.meta.subarray.name = "FULL"
    model.meta.subarray.xstart = 1  # 1-indexed FITS convention
    model.meta.subarray.ystart = 1
    model.meta.subarray.xsize = MIRI_DETECTOR_SHAPE[1]  # nx
    model.meta.subarray.ysize = MIRI_DETECTOR_SHAPE[0]  # ny

    return model


@pytest.fixture
def slitless_taq_image(tmp_path):
    """Generate a slitless TAQ image for testing."""
    wavelength = _get_wavelength("F1500W")
    offset = (2, -3)
    data = make_slitless_data(wavelength, offset)

    # Add NaN values near the slitless source position to test that this still works ok
    data[int(Y_REF_SLITLESS) + 5, int(X_REF_SLITLESS) + 3] = np.nan
    data[int(Y_REF_SLITLESS) - 4, int(X_REF_SLITLESS) - 2] = np.inf

    model = _make_taq_model(data)

    filepath = tmp_path / "slitless_taq.fits"
    model.save(filepath)
    return str(filepath)


def make_slit_data(offset, pathloss_file):
    """
    Make a fake MIRI LRS slit TAQ verification image.

    The output data are zero everywhere except in the slit region,
    where all data are weighted according to the pathloss correction.
    An Airy disk PSF is created at the source position plus offset,
    noise is added, and then the weights are applied.

    Parameters
    ----------
    offset : tuple
        (x,y) offset of the source from the reference center in pixels.
    pathloss_file : str
        Path to the mock pathloss reference file.

    Returns
    -------
    model : ImageModel
        Simulated slit TAQ image model with PSF weighted by pathloss.
    """
    wavelength = _get_wavelength("F1500W")

    # Calculate diffraction-limited Airy disk radius
    theta_rad = 1.22 * wavelength * 1e-6 / JWST_DIAMETER
    theta_arcsec = theta_rad * 206265  # radians to arcseconds
    radius_pixels = theta_arcsec / PIXSCALE  # to pixels

    # Calculate source center position using slit reference position
    x_center = X_REF_SLIT + offset[0]
    y_center = Y_REF_SLIT + offset[1]

    # Create Airy disk model on the full detector
    y, x = np.mgrid[0 : MIRI_DETECTOR_SHAPE[0], 0 : MIRI_DETECTOR_SHAPE[1]]
    airy = models.AiryDisk2D(amplitude=1000.0, x_0=x_center, y_0=y_center, radius=radius_pixels)
    data_full = airy(x, y)

    # Retrieve the pathloss model
    pathloss_model = dm.MirLrsPathlossModel(pathloss_file)
    pathloss_table = pathloss_model.pathloss_table
    pathloss_wcs = pathloss_model.meta.wcsinfo

    # Only calculate pathloss in a region around the slit reference point
    slit_half_height = 15  # pixels in y direction (cross-dispersion)
    slit_half_width = 60  # pixels in x direction (along slit)
    y_min = max(0, int(Y_REF_SLIT - slit_half_height))
    y_max = min(MIRI_DETECTOR_SHAPE[0], int(Y_REF_SLIT + slit_half_height))
    x_min = max(0, int(X_REF_SLIT - slit_half_width))
    x_max = min(MIRI_DETECTOR_SHAPE[1], int(X_REF_SLIT + slit_half_width))

    # Calculate pathloss only in the slit region. calculate_pathloss_vector is unfortunately not vectorized.
    slit_weights = np.zeros_like(data_full)
    for i in range(y_min, y_max):
        for j in range(x_min, x_max):
            # Calculate offset from reference center in arcsec
            dx_arcsec = (j - X_REF_SLIT) * PIXSCALE
            dy_arcsec = (i - Y_REF_SLIT) * PIXSCALE

            _, pathloss_vector, _ = calculate_pathloss_vector(
                pathloss_table["pathloss"], pathloss_wcs, dx_arcsec, dy_arcsec, calc_wave=False
            )
            # Find the wavelength index closest to our wavelength
            wave_idx = np.argmin(np.abs(pathloss_table["wavelength"] - wavelength))
            slit_weights[i, j] = pathloss_vector[wave_idx]

    # Add noise
    rng = np.random.default_rng(42)
    noise = rng.normal(0, 1, data_full.shape)
    data_full += noise

    # Apply slit weights to the data (multiply by pathloss correction)
    data = data_full * slit_weights

    # Add bad values near the slit source position to test their handling
    data[int(Y_REF_SLIT) + 4, int(X_REF_SLIT) + 2] = np.nan
    data[int(Y_REF_SLIT) - 3, int(X_REF_SLIT) - 4] = np.inf

    model = _make_taq_model(data)
    return model


def _make_input_model():
    """Assign necessary metadata shared between slit and slitless modes for LRS."""
    model = dm.ImageModel()
    model.meta.instrument.name = "MIRI"
    model.meta.target.source_type = "POINT"
    model.meta.instrument.filter = "F1500W"
    model.meta.instrument.detector = "MIRIMAGE"

    model.meta.observation.date = "2025-11-10"
    model.meta.observation.time = "12:00:00"
    return model


@pytest.fixture
def input_model_slit():
    """
    Make a mock LRS dataset in slit mode to use as input to the step.

    The data itself are not used by the step; we just need enough metadata to
    retrieve the appropriate reference files.
    """
    model = _make_input_model()
    model.meta.exposure.type = "MIR_LRS-FIXEDSLIT"
    return model


@pytest.fixture
def input_model_slitless():
    """
    Make a mock LRS dataset in slitless mode to use as input to the step.

    The data itself are not used by the step; we just need enough metadata to
    retrieve the appropriate reference files.
    """
    model = _make_input_model()
    model.meta.exposure.type = "MIR_LRS-SLITLESS"
    return model


def test_ta_center_slitless(input_model_slitless, slitless_taq_image, mock_references):
    # Run the TA centering algorithm
    result = TACenterStep.call(input_model_slitless, ta_file=slitless_taq_image)
    x_center, y_center = result.source_xpos, result.source_ypos

    # Expected center position (reference position + offset)
    expected_x = X_REF_SLITLESS + 2
    expected_y = Y_REF_SLITLESS - 3

    # Check that the computed center is close to the expected position
    assert np.isclose(x_center, expected_x, atol=0.05), (
        f"X center {x_center} not close to expected {expected_x}"
    )
    assert np.isclose(y_center, expected_y, atol=0.05), (
        f"Y center {y_center} not close to expected {expected_y}"
    )


@pytest.mark.parametrize(
    "offset",
    [
        (0.0, 0.0),  # Perfectly centered
        (4.0, 0.0),  # Offset to the right
        (-3.0, 0.0),  # Offset to the left
        (0.0, 2.9),  # Offset upward
        (0.0, -3.1),  # Offset downward
        (3.0, -3.0),  # Offset right and down
        (-4.0, 3.6),  # Offset left and up
    ],
)
def test_ta_center_slit(input_model_slit, offset, tmp_path, mock_pathloss_model, mock_references):
    """
    Test TA centering for LRS slit mode with various offsets.

    Some offsets are so large the center of the PSF is slightly outside the slit.
    We need to ensure this case still works reasonably well.
    """

    # Generate slit data with the specified offset
    taq_image = make_slit_data(offset=offset, pathloss_file=mock_pathloss_model)

    # Save to file
    filepath = tmp_path / f"slit_taq_{offset[0]}_{offset[1]}.fits"
    taq_image.save(filepath)

    # Run the TA centering algorithm for slit mode
    result = TACenterStep.call(input_model_slit, ta_file=str(filepath))
    x_center, y_center = result.source_xpos, result.source_ypos

    # Expected center position (reference position + offset)
    expected_x = X_REF_SLIT + offset[0]
    expected_y = Y_REF_SLIT + offset[1]

    # Check that the computed center is close to the expected position
    # The slit truncation may cause slightly larger errors, so use slightly larger tolerance
    assert np.isclose(x_center, expected_x, atol=0.05), (
        f"Offset {offset}: X center {x_center:.2f} not close to expected {expected_x:.2f}"
    )
    assert np.isclose(y_center, expected_y, atol=0.05), (
        f"Offset {offset}: Y center {y_center:.2f} not close to expected {expected_y:.2f}"
    )


def test_skip_no_ta_file(input_model_slit):
    """Test that step is skipped when no TA file is provided."""
    result = TACenterStep.call(input_model_slit, ta_file=None)
    assert result.meta.cal_step.ta_center == "SKIPPED"
    assert not result.hasattr("source_xpos")
    assert not result.hasattr("source_ypos")


def test_skip_extended_source(input_model_slit, slitless_taq_image):
    """Test that step is skipped for extended sources."""
    input_model_slit.meta.target.source_type = "EXTENDED"
    result = TACenterStep.call(input_model_slit, ta_file=slitless_taq_image)
    assert result.meta.cal_step.ta_center == "SKIPPED"
    assert not result.hasattr("source_xpos")
    assert not result.hasattr("source_ypos")


def test_skip_wrong_exp_type(input_model_slit, slitless_taq_image):
    """Test that step is skipped for unsupported exposure types."""
    input_model_slit.meta.exposure.type = "MIR_IMAGE"
    result = TACenterStep.call(input_model_slit, ta_file=slitless_taq_image)
    assert result.meta.cal_step.ta_center == "SKIPPED"
    assert not result.hasattr("source_xpos")
    assert not result.hasattr("source_ypos")


def test_skip_unknown_filter(input_model_slit, slitless_taq_image, tmp_path):
    """Test that step is skipped for unknown filter."""
    # Create a TA model with an unknown filter
    ta_path = tmp_path / "ta_unknown_filter.fits"
    with dm.open(slitless_taq_image) as ta_model:
        ta_model.meta.instrument.filter = "N/A"
        ta_model.save(str(ta_path))

    result = TACenterStep.call(input_model_slit, ta_file=str(ta_path))
    assert result.meta.cal_step.ta_center == "SKIPPED"
    assert not result.hasattr("source_xpos")
    assert not result.hasattr("source_ypos")


def test_ta_center_from_association(
    input_model_slit, tmp_cwd, mock_pathloss_model, mock_references
):
    """Test TA centering step when run on an association with science and TA exposures."""
    # Generate slit TA data with a known offset
    offset = (2.0, -1.5)
    taq_image = make_slit_data(offset=offset, pathloss_file=mock_pathloss_model)

    # Make an association from the science and TA images
    sci_fname = "science.fits"
    input_model_slit.save(sci_fname)
    ta_fname = "ta_image.fits"
    taq_image.save(ta_fname)
    asn = {
        "asn_type": "spec2",
        "asn_id": "test_ta",
        "asn_pool": "pool_id",
        "products": [
            {
                "name": "test_product",
                "members": [
                    {"expname": sci_fname, "exptype": "science"},
                    {"expname": ta_fname, "exptype": "target_acquisition"},
                ],
            }
        ],
    }
    asn_fname = "mir_lrs_ta_asn.json"
    with open(asn_fname, "w") as f:
        json.dump(asn, f)

    # Run the step on the association
    result = TACenterStep.call(asn_fname)

    # Check that the result is the science exposure with TA centering applied
    assert result.meta.cal_step.ta_center == "COMPLETE"

    # Expected center position (reference position + offset)
    expected_x = X_REF_SLIT + offset[0]
    expected_y = Y_REF_SLIT + offset[1]

    # Check that the computed center is close to the expected position
    assert np.isclose(result.source_xpos, expected_x, atol=0.05), (
        f"X center {result.source_xpos:.2f} not close to expected {expected_x:.2f}"
    )
    assert np.isclose(result.source_ypos, expected_y, atol=0.05), (
        f"Y center {result.source_ypos:.2f} not close to expected {expected_y:.2f}"
    )
