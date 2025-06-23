import numpy as np
import pytest
from astropy.table import QTable
from stdatamodels.jwst.datamodels import TsoPhotModel

from jwst.stpipe import Step
from jwst.tso_photometry.tests.test_tso_photometry import mock_nircam_image
from jwst.tso_photometry.tso_photometry_step import TSOPhotometryStep


@pytest.fixture()
def mock_tsophot_reffile(tmp_path, monkeypatch):
    model = TsoPhotModel()
    model.radii = [{"pupil": "ANY", "radius": 6.0, "radius_inner": 8.0, "radius_outer": 11.0}]
    required = ["description", "author", "pedigree", "useafter"]
    for key in required:
        setattr(model.meta, key, "TEST")
    model.meta.instrument.name = "NIRCAM"
    model.meta.exposure.type = "NRC_TSIMAGE"
    reffile = str(tmp_path / "tso_phot.asdf")
    model.save(reffile)

    def mock_get_reference_file(self, dataset, reference_file_type):
        if reference_file_type == "tsophot":
            return reffile
        else:
            return Step().get_reference_file(dataset, reference_file_type)

    monkeypatch.setattr(TSOPhotometryStep, "get_reference_file", mock_get_reference_file)


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
        TSOPhotometryStep.call(datamodel)


def test_tsophotometry_step_subarray(mock_tsophot_reffile, log_watcher):
    datamodel = mock_nircam_image()

    watcher = log_watcher("stpipe.TSOPhotometryStep", message="Extracting gain subarray")
    catalog = TSOPhotometryStep.call(datamodel)
    watcher.assert_seen()

    # Output is a table
    assert isinstance(catalog, QTable)

    # Spot check some values to make sure they are as expected
    assert np.isclose(catalog.meta["xcenter"], datamodel.meta.wcsinfo.siaf_xref_sci - 1, atol=0.01)
    assert np.isclose(catalog.meta["ycenter"], datamodel.meta.wcsinfo.siaf_yref_sci - 1, atol=0.01)
    assert np.allclose(catalog["aperture_sum"].value, 1263.4778, rtol=1.0e-7)


def test_tsophotometry_step_full_frame(mock_tsophot_reffile, log_watcher):
    datamodel = mock_nircam_image(shape=(7, 2048, 2048))

    # Gain reference already matches data, no need to extract subarray
    watcher = log_watcher("stpipe.TSOPhotometryStep", message="Extracting gain subarray")
    catalog = TSOPhotometryStep.call(datamodel)
    watcher.assert_not_seen()

    # Results are otherwise the same as for the subarray
    assert isinstance(catalog, QTable)

    # Spot check some values to make sure they are as expected
    assert np.isclose(catalog.meta["xcenter"], datamodel.meta.wcsinfo.siaf_xref_sci - 1, atol=0.01)
    assert np.isclose(catalog.meta["ycenter"], datamodel.meta.wcsinfo.siaf_yref_sci - 1, atol=0.01)
    assert np.allclose(catalog["aperture_sum"].value, 1263.4778, rtol=1.0e-7)


def test_tsophotometry_step_save_catalog(mock_tsophot_reffile, tmp_path, log_watcher):
    datamodel = mock_nircam_image()
    datamodel.meta.filename = "test_calints.fits"

    watcher = log_watcher("stpipe.TSOPhotometryStep", message="Wrote TSO photometry catalog")
    TSOPhotometryStep.call(
        datamodel,
        save_catalog=True,
        output_dir=str(tmp_path),
    )
    watcher.assert_seen()

    # Catalog is saved
    cat_file = tmp_path / "test_phot.ecsv"
    assert cat_file.exists()

    # Spot check some values to make sure they are as expected
    catalog = QTable.read(cat_file)
    assert np.isclose(catalog.meta["xcenter"], datamodel.meta.wcsinfo.siaf_xref_sci - 1, atol=0.01)
    assert np.isclose(catalog.meta["ycenter"], datamodel.meta.wcsinfo.siaf_yref_sci - 1, atol=0.01)
    assert np.allclose(catalog["aperture_sum"].value, 1263.4778, rtol=1.0e-7)
