"""Test suite for ensuring correct FGS pointing"""
import logging
from numpy import array
from numpy import isclose
import pytest

from jwst.datamodels import Level1bModel
from jwst.lib import engdb_mast
from jwst.lib import set_telescope_pointing as stp

# Set logging for the module to be tested.
logger = logging.getLogger(stp.__name__)
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setLevel(logging.INFO)
logger.addHandler(handler)

DATA_PATH = Path(__file__).parent / 'data'
siaf_path = DATA_PATH / 'xml_data_siafxml' / 'SIAFXML'

# Define minimal model meta structure
WCS_META = {
    'meta': {
        'aperture': {
            'name': 'FGS1_FULL',
        },
        'exposure': {
            'type': 'FGS_ACQ1',
        },
        'guidestar': {
            'gs_ra': 45.1234,
            'gs_dec': -45.1234,
        },
        'observation': {
            'date_beg': '2022-05-23T00:36:08.000',
            'date_end': '2022-05-23T00:36:22.480',
            'date': '2017-01-01',
        },
    }
}


def test_fgs_pointing(mast):
    model = make_level1b()

    # Update wcs
    stp.update_wcs(model, siaf_path=siaf_db, engdb_url=engdb_mast.MAST_BASE_URL)

    # Test results
    assert isclose(model.meta.wcsinfo.pc1_1, -0.9997617158628777, atol=1e-15)
    assert isclose(model.meta.wcsinfo.pc1_2, -0.021829143247382235, atol=1e-15)
    assert isclose(model.meta.wcsinfo.pc2_1, -0.021829143247382235, atol=1e-15)
    assert isclose(model.meta.wcsinfo.pc2_2, 0.9997617158628777, atol=1e-15)

    assert isclose(model.meta.wcsinfo.crpix1, 17.80508537294, atol=1e-15)
    assert isclose(model.meta.wcsinfo.crpix2, 42.65214462281999, atol=1e-15)
    assert isclose(model.meta.wcsinfo.crval1, 45.1234, atol=1e-15)
    assert isclose(model.meta.wcsinfo.crval2, -45.1234, atol=1e-15)


# ---------
# Utilities
# ---------
@pytest.fixture
def mast():

    # See if access to MAST is available.
    try:
        engdb = engdb_mast.EngdbMast(base_url=engdb_mast.MAST_BASE_URL)
    except RuntimeError as exception:
        pytest.skip(f'Live MAST Engineering Service not available: {exception}')

    return engdb


def make_level1b():
    data = array([1.])
    data.shape = (1, 1, 1, 1)
    model = Level1bModel(data)
    model.update(WCS_META)
    return model
