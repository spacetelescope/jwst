import numpy as np
import pytest
import stdatamodels.jwst.datamodels as dm

from jwst.datamodels import ModelContainer
from jwst.targ_centroid.targ_centroid_step import TargCentroidStep
from jwst.targ_centroid.tests.helpers import (
    FULL_SHAPE,
    SLIT_REF,
    SLITLESS_REF,
    make_empty_lrs_model,
    make_slit_data,
    make_slitless_data,
    make_ta_association,
    make_ta_model,
)


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
def mock_references(monkeypatch, mock_filteroffset_model):
    """
    Monkeypatch the get_reference_file method to return mock reference files.

    This allows tests to run without requiring CRDS access, and ensures that
    these tests won't start breaking if the values in the reference file
    are changed in CRDS.
    """

    def mock_get_reference_file(self, model, reftype):
        """Mock implementation of get_reference_file."""
        if reftype == "filteroffset":
            return mock_filteroffset_model
        else:
            raise ValueError(f"Unexpected reference type: {reftype}")

    monkeypatch.setattr(TargCentroidStep, "get_reference_file", mock_get_reference_file)


@pytest.fixture
def slitless_ta_image(tmp_path):
    """Generate a slitless TA image for testing."""
    wavelength = 15.0
    offset = (2, -3)
    data = make_slitless_data(wavelength, offset)

    model = make_ta_model(data, "MIR_LRS-SLITLESS")

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
    return make_empty_lrs_model("MIR_LRS-FIXEDSLIT")


@pytest.fixture
def input_model_slitless():
    """
    Make a mock LRS dataset in slitless mode to use as input to the step.

    The data itself are not used by the step; we just need enough metadata to
    retrieve the appropriate reference files.
    """
    return make_empty_lrs_model("MIR_LRS-SLITLESS")


def test_targ_centroid_slitless(input_model_slitless, slitless_ta_image, mock_references):
    # Run the TA centering algorithm
    result = TargCentroidStep.call(input_model_slitless, ta_file=slitless_ta_image)
    x_center, y_center = result.source_xpos, result.source_ypos

    # Expected center position (reference position + offset)
    expected_x = SLITLESS_REF[0] + 2
    expected_y = SLITLESS_REF[1] - 3

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
def test_targ_centroid_slit(input_model_slit, offset, tmp_path, mock_references):
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
    result = TargCentroidStep.call(input_model_slit, ta_file=str(filepath))
    x_center, y_center = result.source_xpos, result.source_ypos

    # Expected center position (reference position + offset)
    expected_x = SLIT_REF[0] + offset[0]
    expected_y = SLIT_REF[1] + offset[1]

    # Check that the computed center is close to the expected position
    # The slit truncation may cause slightly larger errors, so use slightly larger tolerance
    assert np.isclose(x_center, expected_x, atol=0.05), (
        f"Offset {offset}: X center {x_center:.2f} not close to expected {expected_x:.2f}"
    )
    assert np.isclose(y_center, expected_y, atol=0.05), (
        f"Offset {offset}: Y center {y_center:.2f} not close to expected {expected_y:.2f}"
    )


def test_targ_centroid_flat_flux(input_model_slit, tmp_path, mock_references, log_watcher):
    """
    Test TA centering on an image with flat flux (no source).

    The step should still complete without error, but the centroid position
    may be arbitrary.
    """

    # Create a flat TA image
    data = np.full(FULL_SHAPE, 1000.0)
    ta_model = make_ta_model(data, "MIR_LRS-FIXEDSLIT")

    # Save to file
    filepath = tmp_path / "flat_ta.fits"
    ta_model.save(filepath)

    # Run the TA centering algorithm
    watcher = log_watcher(
        "jwst.targ_centroid.targ_centroid_step",
        message="2D Gaussian centroid fit failed",
    )
    result = TargCentroidStep.call(input_model_slit, ta_file=str(filepath))
    watcher.assert_seen()
    _tests_for_skipped_step(result)


def test_warn_not_point_source(
    input_model_slitless, slitless_ta_image, mock_references, log_watcher
):
    """Test that a warning is logged when the source is not a point source, but the step still runs."""
    input_model_slitless.meta.target.source_type = "EXTENDED"

    watcher = log_watcher(
        "jwst.targ_centroid.targ_centroid_step",
        message="TargCentroidStep is only intended for point sources",
    )
    result = TargCentroidStep.call(input_model_slitless, ta_file=slitless_ta_image)
    watcher.assert_seen()
    assert result.meta.cal_step.targ_centroid == "COMPLETE"


def test_skip_no_ta_file(input_model_slit):
    """Test that step is skipped when no TA file is provided."""
    result = TargCentroidStep.call(input_model_slit, ta_file=None)
    _tests_for_skipped_step(result)


def test_skip_wrong_exp_type(input_model_slit, slitless_ta_image):
    """Test that step is skipped for unsupported exposure types."""
    input_model_slit.meta.exposure.type = "MIR_IMAGE"
    result = TargCentroidStep.call(input_model_slit, ta_file=slitless_ta_image)
    _tests_for_skipped_step(result)


def test_skip_bad_ta_type(input_model_slit, tmp_path):
    """Test that step is skipped when TA file is not an image."""
    # Create a non-image TA model (e.g. CubeModel)
    data = np.zeros((5, FULL_SHAPE[0], FULL_SHAPE[1]))
    ta_model = dm.CubeModel(data)
    ta_model.meta.exposure.type = "MIR_TACONFIRM"

    ta_path = tmp_path / "ta_bad_type.fits"
    ta_model.save(str(ta_path))

    result = TargCentroidStep.call(input_model_slit, ta_file=str(ta_path))
    _tests_for_skipped_step(result)


def test_skip_bad_ta_exptype(input_model_slit, tmp_path):
    """Test that step is skipped when TA file has wrong exposure type."""
    # Create a non-image TA model (e.g. CubeModel)
    data = np.zeros((FULL_SHAPE[0], FULL_SHAPE[1]))
    ta_model = dm.ImageModel(data)
    ta_model.meta.exposure.type = "MIR_IMAGE"

    ta_path = tmp_path / "ta_bad_type.fits"
    ta_model.save(str(ta_path))

    result = TargCentroidStep.call(input_model_slit, ta_file=str(ta_path))
    _tests_for_skipped_step(result)


def test_skip_mostly_nan(input_model_slit, tmp_path, mock_references, log_watcher):
    """Test that step raises an error when center-finding does not converge."""
    # Create a TA model with a source far from the reference position
    data = np.zeros(FULL_SHAPE) * np.nan
    data[int(SLIT_REF[1]) + 1, int(SLIT_REF[0]) + 1] = 1

    ta_model = make_ta_model(data, "MIR_LRS-FIXEDSLIT")

    ta_path = tmp_path / "ta_nonconverge.fits"
    ta_model.save(str(ta_path))

    watcher = log_watcher(
        "jwst.targ_centroid.targ_centroid_step", message="Not enough finite pixels in the cutout"
    )
    result = TargCentroidStep.call(input_model_slit, ta_file=str(ta_path))
    watcher.assert_seen()

    _tests_for_skipped_step(result)


def test_targ_centroid_asn(input_model_slit, tmp_cwd, mock_references):
    """Test TA centering step when run on an association with science and TA exposures."""
    # Generate slit TA data with a known offset
    offset = (2.0, -1.5)
    ta_image = make_slit_data(offset=offset)

    # Create association
    asn_fname = make_ta_association(input_model_slit, ta_image)

    # Run the step on the association
    result = TargCentroidStep.call(asn_fname)

    assert isinstance(result, ModelContainer)
    sci_idx = result.ind_asn_type("science")
    sci_model = result[sci_idx[0]]

    # Check that the result is the science exposure with TA centering applied
    assert sci_model.meta.cal_step.targ_centroid == "COMPLETE"

    # Expected center position (reference position + offset)
    expected_x = SLIT_REF[0] + offset[0]
    expected_y = SLIT_REF[1] + offset[1]

    # Check that the computed center is close to the expected position
    assert np.isclose(sci_model.source_xpos, expected_x, atol=0.05), (
        f"X center {sci_model.source_xpos:.2f} not close to expected {expected_x:.2f}"
    )
    assert np.isclose(sci_model.source_ypos, expected_y, atol=0.05), (
        f"Y center {sci_model.source_ypos:.2f} not close to expected {expected_y:.2f}"
    )


def test_skip_asn_no_ta(input_model_slit, tmp_cwd, mock_references):
    """Test TA centering step when run on an association with no TA exposure."""
    # Create association with only science exposure
    asn_fname = make_ta_association(input_model_slit, ta_model=None)

    # Run the step on the association
    result = TargCentroidStep.call(asn_fname)

    assert isinstance(result, ModelContainer)
    sci_idx = result.ind_asn_type("science")
    sci_model = result[sci_idx[0]]

    # Check that the step was skipped
    _tests_for_skipped_step(sci_model)


def _tests_for_skipped_step(model):
    assert model.meta.cal_step.targ_centroid == "SKIPPED"
    assert not model.hasattr("source_xpos")
    assert not model.hasattr("source_ypos")


class PlaceholderWCS:
    """Simple WCS that applies a linear transformation for testing."""

    def __init__(self, xscale=100.0, yscale=100.0, xoff=1.0, yoff=2.0):
        self.xscale = xscale
        self.yscale = yscale
        self.xoff = xoff
        self.yoff = yoff

    def get_transform(self, frm, to):
        def world_to_pixel(ra, dec):
            # simple linear transform for testing
            return (ra * self.xscale + self.xoff, dec * self.yscale + self.yoff)

        return world_to_pixel


def test_skip_assign_wcs_error(input_model_slit, tmp_path, monkeypatch, log_watcher):
    """If AssignWcsStep.call raises, the step should be skipped."""

    # Create a simple TA model file to pass into the step
    ta_model = make_ta_model(np.zeros((10, 10)), "MIR_LRS-FIXEDSLIT")
    ta_path = tmp_path / "ta_assign_fail.fits"
    ta_model.save(str(ta_path))

    def mock_assign_wcs(ta_model, **kwargs):
        raise RuntimeError("assign_wcs failed")

    monkeypatch.setattr(
        "jwst.targ_centroid.targ_centroid.AssignWcsStep.call",
        mock_assign_wcs,
    )

    watcher = log_watcher(
        "jwst.targ_centroid.targ_centroid_step",
        message="Error when assigning WCS",
    )

    result = TargCentroidStep.call(input_model_slit, ta_file=str(ta_path))
    watcher.assert_seen()

    _tests_for_skipped_step(result)


def test_skip_assign_wcs_skipped(input_model_slit, tmp_path, monkeypatch, log_watcher):
    """If AssignWcsStep.call returns but assign_wcs is skipped, step is skipped."""

    # Create a simple TA model file to pass into the step
    ta_model = make_ta_model(np.zeros((10, 10)), "MIR_LRS-FIXEDSLIT")
    ta_path = tmp_path / "ta_assign_incomplete.fits"
    ta_model.save(str(ta_path))

    def mock_assign_wcs(ta_model, **kwargs):
        ta_model.meta.cal_step.assign_wcs = "SKIPPED"
        return ta_model

    monkeypatch.setattr(
        "jwst.targ_centroid.targ_centroid.AssignWcsStep.call",
        mock_assign_wcs,
    )

    watcher = log_watcher("jwst.targ_centroid.targ_centroid_step", message="Failed to assign WCS")

    result = TargCentroidStep.call(input_model_slit, ta_file=str(ta_path))
    watcher.assert_seen()

    _tests_for_skipped_step(result)
