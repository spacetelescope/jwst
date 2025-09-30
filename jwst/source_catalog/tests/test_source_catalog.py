from pathlib import Path

import numpy as np
import pytest
import stdatamodels.jwst.datamodels as dm
from astropy.table import QTable
from numpy.testing import assert_allclose
from photutils.datasets import make_gwcs

from jwst.source_catalog import SourceCatalogStep


@pytest.fixture
def nircam_model():
    rng = np.random.default_rng(seed=123)
    data = rng.normal(0, 0.5, size=(101, 101))
    data[20:80, 10:20] = 1.4
    data[20:30, 20:45] = 1.4
    data[20:80, 55:65] = 7.2
    data[70:80, 65:87] = 7.2
    data[45:55, 65:87] = 7.2
    data[20:30, 65:87] = 7.2
    data[55:75, 82:92] = 7.2
    data[25:45, 82:92] = 7.2

    wht = np.ones(data.shape)
    wht[0:10, :] = 0.0
    err = np.abs(data) / 10.0
    model = dm.ImageModel(data, wht=wht, err=err)
    model.meta.bunit_data = "MJy/sr"
    model.meta.bunit_err = "MJy/sr"
    model.meta.photometry.pixelarea_steradians = 1.0
    model.meta.wcs = make_gwcs(data.shape)
    model.meta.wcsinfo = {
        "ctype1": "RA---TAN",
        "ctype2": "DEC--TAN",
        "dec_ref": 11.99875540218638,
        "ra_ref": 22.02351763251896,
        "roll_ref": 0.005076934167039675,
        "v2_ref": 86.039011,
        "v3_ref": -493.385704,
        "v3yangle": -0.07385127,
        "vparity": -1,
        "wcsaxes": 2,
        "crpix1": 50,
        "crpix2": 50,
    }
    model.meta.instrument = {
        "channel": "LONG",
        "detector": "NRCALONG",
        "filter": "F444W",
        "lamp_mode": "NONE",
        "module": "A",
        "name": "NIRCAM",
        "pupil": "CLEAR",
    }
    model.meta.exposure.type = "NRC_IMAGE"
    model.meta.observation.date = "2021-01-01"
    model.meta.observation.time = "00:00:00"

    return model


@pytest.fixture
def nircam_model_without_apcorr():
    rng = np.random.default_rng(seed=123)
    data = rng.normal(0, 0.5, size=(101, 101))
    data[20:80, 10:20] = 1.4
    data[20:30, 20:45] = 1.4
    data[20:80, 55:65] = 7.2
    data[70:80, 65:87] = 7.2
    data[45:55, 65:87] = 7.2
    data[20:30, 65:87] = 7.2
    data[55:75, 82:92] = 7.2
    data[25:45, 82:92] = 7.2

    wht = np.ones(data.shape)
    wht[0:10, :] = 0.0
    err = np.abs(data) / 10.0
    model = dm.ImageModel(data, wht=wht, err=err)
    model.meta.bunit_data = "MJy/sr"
    model.meta.bunit_err = "MJy/sr"
    model.meta.photometry.pixelarea_steradians = 1.0
    model.meta.wcs = make_gwcs(data.shape)
    model.meta.wcsinfo = {
        "ctype1": "RA---TAN",
        "ctype2": "DEC--TAN",
        "dec_ref": 11.99875540218638,
        "ra_ref": 22.02351763251896,
        "roll_ref": 0.005076934167039675,
        "v2_ref": 86.039011,
        "v3_ref": -493.385704,
        "v3yangle": -0.07385127,
        "vparity": -1,
        "wcsaxes": 2,
        "crpix1": 50,
        "crpix2": 50,
    }
    model.meta.instrument = {
        "channel": "LONG",
        "detector": "NRCALONG",
        "filter": "F2550WR",
        "lamp_mode": "NONE",
        "module": "A",
        "name": "NIRCAM",
        "pupil": "CLEAR",
    }
    model.meta.exposure.type = "NRC_IMAGE"
    model.meta.observation.date = "2021-01-01"
    model.meta.observation.time = "00:00:00"

    return model


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
        assert np.isclose(cat["xcentroid"][1], 19.46399720865899)
        assert np.isclose(cat["ycentroid"][1], 41.95288393407728)
        assert np.isclose(cat["aper_bkg_flux"][1].value, 1400000.0)
        assert np.isclose(cat["aper_bkg_flux_err"][1].value, 85223.70700074881)
        assert np.isclose(cat["CI_50_30"][1], 2.3342599432074653)
        assert np.isclose(cat["sharpness"][1], 0.9102634628764403)
        assert np.isclose(cat["roundness"][1], 1.5954264)
        assert np.isclose(cat["nn_dist"][1].value, 53.0737632103816)
        assert np.isclose(cat["isophotal_flux"][1], 930.9999841451645)
        assert np.isclose(cat["isophotal_flux_err"][1], 3.6102633)
        assert np.isclose(cat["semimajor_sigma"][1].value, 18.847635525516534)
        assert np.isclose(cat["semiminor_sigma"][1].value, 7.031371175038476)
        assert np.isclose(cat["ellipticity"][1], 0.626936165784871)
        assert np.isclose(cat["orientation"][1].value, -72.75413766990114)


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

    # make a point source with some size that looks a bit like a psf, no need to be realistic
    point_source = np.ones((7, 7))
    point_source[1:6, 1:6] = 3.0
    point_source[2:5, 2:5] = 5.0
    point_source[3, 3] = 10.0

    data[30:37, 30:37] = point_source
    data[70:77, 70:77] = point_source

    nircam_model.data = data

    nircam_model.err = np.abs(data) / 10.0
    nircam_model.wht = np.ones(data.shape)

    step = SourceCatalogStep(
        snr_threshold=8.0,
        npixels=5,
        bkg_boxsize=50,
        kernel_fwhm=2.0,
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
