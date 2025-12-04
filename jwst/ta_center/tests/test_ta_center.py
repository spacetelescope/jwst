import numpy as np
import pytest
import stdatamodels.jwst.datamodels as dm

from jwst.datamodels import ModelContainer
from jwst.ta_center.ta_center_step import TACenterStep, _get_wavelength
from jwst.ta_center.tests.helpers import (
    X_REF_SLIT,
    X_REF_SLITLESS,
    Y_REF_SLIT,
    Y_REF_SLITLESS,
    make_empty_lrs_model,
    make_pathloss_model,
    make_slit_data,
    make_slitless_data,
    make_ta_association,
    make_ta_model,
)


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
    # Save the model to a file
    pathloss_model = make_pathloss_model()
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


@pytest.fixture
def slitless_ta_image(tmp_path):
    """Generate a slitless TA image for testing."""
    wavelength = _get_wavelength("F1500W")
    offset = (2, -3)
    data = make_slitless_data(wavelength, offset)

    # Add NaN values near the slitless source position to test that this still works ok
    data[int(Y_REF_SLITLESS) + 5, int(X_REF_SLITLESS) + 3] = np.nan
    data[int(Y_REF_SLITLESS) - 4, int(X_REF_SLITLESS) - 2] = np.inf

    model = make_ta_model(data)

    filepath = tmp_path / "slitless_ta.fits"
    model.save(filepath)
    return str(filepath)


@pytest.fixture
def input_model_slit():
    """
    Make a mock LRS dataset in slit mode to use as input to the step.

    The data itself are not used by the step; we just need enough metadata to
    retrieve the appropriate reference files.
    """
    model = make_empty_lrs_model()
    model.meta.exposure.type = "MIR_LRS-FIXEDSLIT"
    return model


@pytest.fixture
def input_model_slitless():
    """
    Make a mock LRS dataset in slitless mode to use as input to the step.

    The data itself are not used by the step; we just need enough metadata to
    retrieve the appropriate reference files.
    """
    model = make_empty_lrs_model()
    model.meta.exposure.type = "MIR_LRS-SLITLESS"
    return model


def test_ta_center_slitless(input_model_slitless, slitless_ta_image, mock_references):
    # Run the TA centering algorithm
    result = TACenterStep.call(input_model_slitless, ta_file=slitless_ta_image)
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
def test_ta_center_slit(input_model_slit, offset, tmp_path, mock_references):
    """
    Test TA centering for LRS slit mode with various offsets.

    Some offsets are so large the center of the PSF is slightly outside the slit.
    We need to ensure this case still works reasonably well.
    """

    # Generate slit data with the specified offset
    ta_image = make_slit_data(offset=offset)

    # Save to file
    filepath = tmp_path / f"slit_ta_{offset[0]}_{offset[1]}.fits"
    ta_image.save(filepath)

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


def test_skip_extended_source(input_model_slit, slitless_ta_image):
    """Test that step is skipped for extended sources."""
    input_model_slit.meta.target.source_type = "EXTENDED"
    result = TACenterStep.call(input_model_slit, ta_file=slitless_ta_image)
    assert result.meta.cal_step.ta_center == "SKIPPED"
    assert not result.hasattr("source_xpos")
    assert not result.hasattr("source_ypos")


def test_skip_wrong_exp_type(input_model_slit, slitless_ta_image):
    """Test that step is skipped for unsupported exposure types."""
    input_model_slit.meta.exposure.type = "MIR_IMAGE"
    result = TACenterStep.call(input_model_slit, ta_file=slitless_ta_image)
    assert result.meta.cal_step.ta_center == "SKIPPED"
    assert not result.hasattr("source_xpos")
    assert not result.hasattr("source_ypos")


def test_skip_unknown_filter(input_model_slit, slitless_ta_image, tmp_path):
    """Test that step is skipped for unknown filter."""
    # Create a TA model with an unknown filter
    ta_path = tmp_path / "ta_unknown_filter.fits"
    with dm.open(slitless_ta_image) as ta_model:
        ta_model.meta.instrument.filter = "N/A"
        ta_model.save(str(ta_path))

    result = TACenterStep.call(input_model_slit, ta_file=str(ta_path))
    assert result.meta.cal_step.ta_center == "SKIPPED"
    assert not result.hasattr("source_xpos")
    assert not result.hasattr("source_ypos")


def test_skip_no_finite_pixels(input_model_slit, tmp_path, mock_references, log_watcher):
    """Test that step raises an error when no finite pixels are present in TA image."""
    # Create a TA model with all non-finite values
    data = np.full((101, 101), np.nan)
    ta_model = make_ta_model(data)

    ta_path = tmp_path / "ta_no_finite.fits"
    ta_model.save(str(ta_path))

    watcher = log_watcher(
        "jwst.ta_center.ta_center_step", message="All pixels contain non-finite values"
    )
    result = TACenterStep.call(input_model_slit, ta_file=str(ta_path))
    watcher.assert_seen()

    assert result.meta.cal_step.ta_center == "SKIPPED"
    assert not result.hasattr("source_xpos")
    assert not result.hasattr("source_ypos")


def test_ta_center_asn(input_model_slit, tmp_cwd, mock_references):
    """Test TA centering step when run on an association with science and TA exposures."""
    # Generate slit TA data with a known offset
    offset = (2.0, -1.5)
    ta_image = make_slit_data(offset=offset)

    # Create association
    asn_fname = make_ta_association(input_model_slit, ta_image)

    # Run the step on the association
    result = TACenterStep.call(asn_fname)

    assert isinstance(result, ModelContainer)
    sci_idx = result.ind_asn_type("science")
    sci_model = result[sci_idx[0]]

    # Check that the result is the science exposure with TA centering applied
    assert sci_model.meta.cal_step.ta_center == "COMPLETE"

    # Expected center position (reference position + offset)
    expected_x = X_REF_SLIT + offset[0]
    expected_y = Y_REF_SLIT + offset[1]

    # Check that the computed center is close to the expected position
    assert np.isclose(sci_model.source_xpos, expected_x, atol=0.05), (
        f"X center {sci_model.source_xpos:.2f} not close to expected {expected_x:.2f}"
    )
    assert np.isclose(sci_model.source_ypos, expected_y, atol=0.05), (
        f"Y center {sci_model.source_ypos:.2f} not close to expected {expected_y:.2f}"
    )


def test_ta_center_asn_no_ta(input_model_slit, tmp_cwd, mock_references):
    """Test TA centering step when run on an association with no TA exposure."""
    # Create association with only science exposure
    asn_fname = make_ta_association(input_model_slit, ta_model=None)

    # Run the step on the association
    result = TACenterStep.call(asn_fname)

    assert isinstance(result, ModelContainer)
    sci_idx = result.ind_asn_type("science")
    sci_model = result[sci_idx[0]]

    # Check that the step was skipped
    assert sci_model.meta.cal_step.ta_center == "SKIPPED"
    assert not sci_model.hasattr("source_xpos")
    assert not sci_model.hasattr("source_ypos")
