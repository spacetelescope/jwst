import pytest
import numpy as np
from photutils.datasets import make_gwcs
from jwst.datamodels import ImageModel
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
    model = ImageModel(data, wht=wht)
    model.meta.bunit_data = 'MJy/sr'
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
