import numpy as np
import pytest
from astropy.table import QTable

from jwst.tso_photometry import tso_photometry_step as tp
from jwst.tso_photometry.tests.test_tso_photometry import mock_nircam_image

RADII_FOR_TEST = {"radius": 6.0, "radius_inner": 8.0, "radius_outer": 11.0}


def spot_check_expected_values(datamodel, catalog):
    """Spot check some output values to make sure they are as expected."""
    assert np.allclose(catalog["aperture_x"].value, datamodel.meta.wcsinfo.siaf_xref_sci - 1)
    assert np.allclose(catalog["aperture_y"].value, datamodel.meta.wcsinfo.siaf_yref_sci - 1)
    assert np.allclose(catalog["aperture_sum"].value, 1263.4778, rtol=1.0e-7)


@pytest.mark.parametrize("missing", ["xref", "yref", "bunit_data", "bunit_err"])
def test_tsophotometry_step_missing_values(missing):
    datamodel = mock_nircam_image()
    if missing == "xref":
        datamodel.meta.wcsinfo.siaf_xref_sci = None
    elif missing == "yref":
        datamodel.meta.wcsinfo.siaf_yref_sci = None
    elif missing == "bunit_data":
        datamodel.meta.bunit_data = None
    elif missing == "bunit_err":
        datamodel.meta.bunit_err = None

    with pytest.raises(ValueError, match="missing"):
        tp.TSOPhotometryStep.call(datamodel)


def test_tsophotometry_step_subarray(log_watcher):
    datamodel = mock_nircam_image()
    input_copy = datamodel.copy()

    watcher = log_watcher(
        "jwst.tso_photometry.tso_photometry_step", message="Extracting gain subarray"
    )
    catalog = tp.TSOPhotometryStep.call(datamodel, **RADII_FOR_TEST)
    watcher.assert_seen()

    # Output is a table
    assert isinstance(catalog, QTable)
    spot_check_expected_values(datamodel, catalog)

    # Input is not modified
    assert catalog is not datamodel
    np.testing.assert_allclose(datamodel.data, input_copy.data)


def test_tsophotometry_step_full_frame(log_watcher):
    datamodel = mock_nircam_image(shape=(7, 2048, 2048))
    input_copy = datamodel.copy()

    # Gain reference already matches data, no need to extract subarray
    watcher = log_watcher(
        "jwst.tso_photometry.tso_photometry_step", message="Extracting gain subarray"
    )
    catalog = tp.TSOPhotometryStep.call(datamodel, **RADII_FOR_TEST)
    watcher.assert_not_seen()

    # Results are otherwise the same as for the subarray
    assert isinstance(catalog, QTable)
    spot_check_expected_values(datamodel, catalog)

    # Input is not modified
    assert catalog is not datamodel
    np.testing.assert_allclose(datamodel.data, input_copy.data)


def test_tsophotometry_step_save_catalog(tmp_path, log_watcher):
    datamodel = mock_nircam_image()
    datamodel.meta.filename = "test_calints.fits"

    watcher = log_watcher(
        "jwst.tso_photometry.tso_photometry_step", message="Wrote TSO photometry catalog"
    )
    tp.TSOPhotometryStep.call(
        datamodel, save_catalog=True, output_dir=str(tmp_path), **RADII_FOR_TEST
    )
    watcher.assert_seen()

    # Catalog is saved
    cat_file = tmp_path / "test_phot.ecsv"
    assert cat_file.exists()

    # Spot check some values to make sure they are as expected
    catalog = QTable.read(cat_file)
    spot_check_expected_values(datamodel, catalog)


def test_tsophotometry_step_no_centroid():
    datamodel = mock_nircam_image()

    catalog = tp.TSOPhotometryStep.call(datamodel, centroid_source=False, **RADII_FOR_TEST)

    # Output is the same as the centroided result for the test data
    assert isinstance(catalog, QTable)
    spot_check_expected_values(datamodel, catalog)


def test_tsophotometry_step_moving_centroid():
    datamodel = mock_nircam_image()

    catalog = tp.TSOPhotometryStep.call(datamodel, moving_centroid=True, **RADII_FOR_TEST)

    # Output is the same as the single centroided result for the test data
    assert isinstance(catalog, QTable)
    spot_check_expected_values(datamodel, catalog)


def test_tsophotometry_step_failed_centroid(monkeypatch, log_watcher):
    datamodel = mock_nircam_image()

    def mock_centroid(*args, **kwars):
        nan_array = np.full(datamodel.shape[0], np.nan)
        return nan_array, nan_array, nan_array, nan_array, nan_array

    monkeypatch.setattr(tp, "tso_source_centroid", mock_centroid)

    watcher = log_watcher("jwst.tso_photometry.tso_photometry_step", message="Centroid fit failed")
    catalog = tp.TSOPhotometryStep.call(datamodel, centroid_source=True, **RADII_FOR_TEST)
    watcher.assert_seen()

    # Output is the same as the centroided result for the test data:
    # failed fits fall back on estimated location.
    assert isinstance(catalog, QTable)
    spot_check_expected_values(datamodel, catalog)


@pytest.mark.parametrize("box", ["search_box_width", "fit_box_width"])
def test_tsophotometry_step_odd_boxes(log_watcher, box):
    datamodel = mock_nircam_image()
    kwargs = {box: 22}
    kwargs.update(RADII_FOR_TEST)

    watcher = log_watcher(
        "jwst.tso_photometry.tso_photometry_step",
        message=f"Rounding the {box} down to 21",
        level="warning",
    )
    catalog = tp.TSOPhotometryStep.call(datamodel, **kwargs)
    watcher.assert_seen()

    # Output is the same for easy-to-fit synthetic data
    assert isinstance(catalog, QTable)
    spot_check_expected_values(datamodel, catalog)
