from pathlib import Path

import numpy as np
import pytest
import stdatamodels.jwst.datamodels as dm
from astropy.modeling import models
from astropy.table import QTable
from numpy.testing import assert_allclose

from jwst.source_catalog import SourceCatalogStep
from jwst.source_catalog.tests.helpers import make_nircam_model, make_nircam_model_without_apcorr


@pytest.fixture
def nircam_model():
    return make_nircam_model()


@pytest.fixture
def nircam_model_without_apcorr():
    return make_nircam_model_without_apcorr()


@pytest.mark.parametrize(
    "npixels, nsources",
    (
        (5, 2),
        (1000, 1),
    ),
)
def test_source_catalog(nircam_model, npixels, nsources):
    step = SourceCatalogStep(
        snr_threshold=0.5, npixels=npixels, bkg_boxsize=50, kernel_fwhm=2.0, save_results=False
    )
    cat = step.run(nircam_model)
    assert len(cat) == nsources
    min_snr = np.min(cat["isophotal_flux"] / cat["isophotal_flux_err"])
    assert min_snr >= 100

    if npixels == 5 and nsources == 2:
        # test values of some specific computed quantities
        assert np.isclose(cat["xcentroid"][1], 19.453064764431833)
        assert np.isclose(cat["ycentroid"][1], 41.963065678485115)
        assert np.isclose(cat["aper_bkg_flux"][1].value, 1.40e-7)
        assert np.isclose(cat["aper_bkg_flux_err"][1].value, 8.52237054e-09)
        assert np.isclose(cat["CI_50_30"][1], 2.352272196434725)
        assert np.isclose(cat["sharpness"][1], 0.9102634628764403)
        assert np.isclose(cat["roundness"][1], 1.5954264)
        assert np.isclose(cat["nn_dist"][1].value, 53.07168773319497)
        assert np.isclose(cat["isophotal_flux"][1].value, 7.940753383195442e-05)
        assert cat["isophotal_flux_err"][1].unit == "Jy"
        assert np.isclose(cat["isophotal_flux_err"][1].value, 3.6102634e-07)
        assert np.isclose(cat["isophotal_abmag"][1], 19.150345729246215)
        assert np.isclose(cat["semimajor_sigma"][1].value, 18.84372169911978)
        assert np.isclose(cat["semiminor_sigma"][1].value, 7.024931388267103)
        assert np.isclose(cat["ellipticity"][1].value, 0.62720043)
        assert np.isclose(cat["orientation"][1].value, -72.78329207140818)


def test_source_catalog_no_sources(nircam_model, monkeypatch):
    npixels = 5000
    step = SourceCatalogStep(
        snr_threshold=0.5, npixels=npixels, bkg_boxsize=50, kernel_fwhm=2.0, save_results=False
    )
    cat = step.run(nircam_model)
    assert cat is None


def test_model_without_apcorr_data(nircam_model_without_apcorr):
    step = SourceCatalogStep(save_results=False)
    cat = step.run(nircam_model_without_apcorr)
    assert cat is None


def test_input_model_reset(nircam_model):
    """Changes to input model data are made in SourceCatalogStep - make sure that
    these changes are correctly reversed so the input model data/err arrays
    remain unchanged after processing (to avoid copying datamodel),
    and that the units are in MJy/sr before and after."""

    original_data = nircam_model.data.copy()
    original_data_unit = nircam_model.meta.bunit_data
    original_err = nircam_model.err.copy()
    original_err_unit = nircam_model.meta.bunit_err

    assert original_data_unit == original_err_unit == "MJy/sr"

    step = SourceCatalogStep(
        snr_threshold=0.5, npixels=5, bkg_boxsize=50, kernel_fwhm=2.0, save_results=False
    )

    step.run(nircam_model)

    assert_allclose(original_data, nircam_model.data, atol=5.0e-7)
    assert_allclose(original_err, nircam_model.err, 5.0e-7)
    assert nircam_model.meta.bunit_data == "MJy/sr"
    assert nircam_model.meta.bunit_err == "MJy/sr"


@pytest.mark.parametrize("finder", ["segmentation", "iraf", "dao"])
def test_source_catalog_point_sources(finder, nircam_model, tmp_cwd):
    """Test the three source finding algorithms with point sources."""
    data = np.random.default_rng(seed=123).normal(0, 0.5, size=(101, 101))

    # Create coordinate grids for Airy disk models
    y, x = np.mgrid[0:101, 0:101]

    # Create two Airy disk point sources
    # First source at position (33, 33)
    airy1 = models.AiryDisk2D(amplitude=10.0, x_0=33, y_0=33, radius=5.0)
    source1 = airy1(x, y)

    # Second source at position (73, 73)
    airy2 = models.AiryDisk2D(amplitude=10.0, x_0=73, y_0=73, radius=5.0)
    source2 = airy2(x, y)

    # Add the sources to the background
    data += source1 + source2

    nircam_model.data = data

    nircam_model.err = np.abs(data) / 10.0
    nircam_model.wht = np.ones(data.shape)

    step = SourceCatalogStep(
        snr_threshold=8.0,
        npixels=5,
        bkg_boxsize=50,
        kernel_fwhm=2.0,
        roundlo=-1.0,
        sharplo=0.0,
        starfinder=finder,
        save_results=True,
    )
    cat = step.run(nircam_model)

    assert len(cat) == 2, f"Expected 2 sources, found {len(cat)} with {finder} finder."

    # check output products were created
    cat_name = "step_SourceCatalogStep_cat.ecsv"
    assert Path(cat_name).exists()

    # test coverage for bug that iraf orientation did not have units
    assert "orientation" in cat.colnames
    if finder != "dao":
        assert cat["orientation"].unit == "deg"
    else:
        assert np.all(np.isnan(cat["orientation"]))

    if finder == "segmentation":
        segm_name = "step_SourceCatalogStep_segm.fits"
        assert Path(segm_name).exists()
        with dm.open(segm_name) as segm:
            assert segm.data.shape == (101, 101)
            assert segm.meta.hasattr("wcs")
            assert segm.meta.hasattr("wcsinfo")


def test_output_is_not_input(nircam_model):
    """Make sure output is not the same as input."""
    step = SourceCatalogStep(
        snr_threshold=0.5, npixels=5, bkg_boxsize=50, kernel_fwhm=2.0, save_results=False
    )
    result = step.run(nircam_model)

    # Input is an image, output is a catalog
    assert result is not nircam_model
    assert isinstance(result, QTable)
    assert isinstance(nircam_model, dm.ImageModel)
