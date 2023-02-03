import pytest
import numpy as np
from numpy.testing import assert_allclose
from photutils.datasets import make_gwcs

from stdatamodels.jwst.datamodels import ImageModel

from ..source_catalog_step import SourceCatalogStep


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
    wht[0:10, :] = 0.
    err = np.abs(data) / 10.
    model = ImageModel(data, wht=wht, err=err)
    model.meta.bunit_data = 'MJy/sr'
    model.meta.bunit_err = 'MJy/sr'
    model.meta.photometry.pixelarea_steradians = 1.0
    model.meta.wcs = make_gwcs(data.shape)
    model.meta.wcsinfo = {
        'ctype1': 'RA---TAN',
        'ctype2': 'DEC--TAN',
        'dec_ref': 11.99875540218638,
        'ra_ref': 22.02351763251896,
        'roll_ref': 0.005076934167039675,
        'v2_ref': 86.039011,
        'v3_ref': -493.385704,
        'v3yangle': -0.07385127,
        'vparity': -1,
        'wcsaxes': 2,
        'crpix1': 50,
        'crpix2': 50}
    model.meta.instrument = {
        'channel': 'LONG',
        'detector': 'NRCALONG',
        'filter': 'F444W',
        'lamp_mode': 'NONE',
        'module': 'A',
        'name': 'NIRCAM',
        'pupil': 'CLEAR'}
    model.meta.exposure.type = 'NRC_IMAGE'
    model.meta.observation.date = '2021-01-01'
    model.meta.observation.time = '00:00:00'

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
    wht[0:10, :] = 0.
    err = np.abs(data) / 10.
    model = ImageModel(data, wht=wht, err=err)
    model.meta.bunit_data = 'MJy/sr'
    model.meta.bunit_err = 'MJy/sr'
    model.meta.photometry.pixelarea_steradians = 1.0
    model.meta.wcs = make_gwcs(data.shape)
    model.meta.wcsinfo = {
        'ctype1': 'RA---TAN',
        'ctype2': 'DEC--TAN',
        'dec_ref': 11.99875540218638,
        'ra_ref': 22.02351763251896,
        'roll_ref': 0.005076934167039675,
        'v2_ref': 86.039011,
        'v3_ref': -493.385704,
        'v3yangle': -0.07385127,
        'vparity': -1,
        'wcsaxes': 2,
        'crpix1': 50,
        'crpix2': 50}
    model.meta.instrument = {
        'channel': 'LONG',
        'detector': 'NRCALONG',
        'filter': 'F2550WR',
        'lamp_mode': 'NONE',
        'module': 'A',
        'name': 'NIRCAM',
        'pupil': 'CLEAR'}
    model.meta.exposure.type = 'NRC_IMAGE'
    model.meta.observation.date = '2021-01-01'
    model.meta.observation.time = '00:00:00'

    return model


@pytest.mark.parametrize('npixels, nsources', ((5, 2), (1000, 1), (5000, 0)))
def test_source_catalog(nircam_model, npixels, nsources):

    step = SourceCatalogStep(snr_threshold=0.5, npixels=npixels,
                             bkg_boxsize=50, kernel_fwhm=2.0,
                             save_results=False)
    cat = step.run(nircam_model)
    if cat is None:
        assert nsources == 0
    else:
        assert len(cat) == nsources
        min_snr = np.min(cat['isophotal_flux'] / cat['isophotal_flux_err'])
        assert min_snr >= 100


def test_model_without_apcorr_data(nircam_model_without_apcorr):
    step = SourceCatalogStep(save_results=False)
    cat = step.run(nircam_model_without_apcorr)
    assert cat is None


def test_input_model_reset(nircam_model):
    """ Changes to input model data are made in SourceCatalogStep - make sure that
        these changes are correctly reversed so the input model data/err arrays
        remain unchanged after processing (to avoid copying datamodel),
        and that the units are in MJy/sr before and after."""

    original_data = nircam_model.data.copy()
    original_data_unit = nircam_model.meta.bunit_data
    original_err = nircam_model.err.copy()
    original_err_unit = nircam_model.meta.bunit_err

    assert (original_data_unit == original_err_unit == 'MJy/sr')

    step = SourceCatalogStep(snr_threshold=0.5, npixels=5,
                             bkg_boxsize=50, kernel_fwhm=2.0,
                             save_results=False)

    step.run(nircam_model)

    assert_allclose(original_data, nircam_model.data, atol=5.e-7)
    assert_allclose(original_err, nircam_model.err, 5.e-7)
    assert (nircam_model.meta.bunit_data == 'MJy/sr')
    assert (nircam_model.meta.bunit_err == 'MJy/sr')
