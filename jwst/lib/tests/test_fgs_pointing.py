"""Test suite for ensuring correct FGS pointing"""
import logging
from numpy import array

from ...datamodels import Level1bModel
from .. import set_telescope_pointing as stp


# Set logging for the module to be tested.
logger = logging.getLogger(stp.__name__)
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setLevel(logging.INFO)
logger.addHandler(handler)

# Define minimal model meta structure
WCS_META = {
    'meta': {
        'exposure': {
            'type': 'FGS_ACQ1',
        },
        'wcsinfo': {
            'cdelt1': 1.90364333333333e-05,
            'cdelt2': 1.94607944444444e-05,
            'crpix1': 1024.5,
            'crpix2': 1024.5,
            'crval1': 0,
            'crval2': 0,
            'ctype1': 'RA---TAN',
            'ctype2': 'DEC--TAN',
            'cunit1': 'deg',
            'cunit2': 'deg',
            'dec_v1': 13.0,
            'pa_v3': 0.0,
            'pc1_1': 1.0,
            'pc1_2': 0.0,
            'pc2_1': 0.0,
            'pc2_2': 1.0,
            'ra_v1': 79.0,
            's_region': '',
            'v2_ref': 0,
            'v3_ref': 0,
            'v3yangle': -1.2508,
            'vparity': -1,
            'wcsaxes': 2
        }
    }
}


def test_fgs_pointing():
    model = make_level1b()
    stp.update_wcs(model)

    assert model.meta.wcsinfo.pc1_1 == -1.0
    assert model.meta.wcsinfo.pc1_2 == 0.0
    assert model.meta.wcsinfo.pc2_1 == 0.0
    assert model.meta.wcsinfo.pc2_2 == 1.0


# ---------
# Utilities
# ---------
def make_level1b():
    data = array([1.])
    data.shape = (1, 1, 1, 1)
    model = Level1bModel(data)
    model.update(WCS_META)
    return model
