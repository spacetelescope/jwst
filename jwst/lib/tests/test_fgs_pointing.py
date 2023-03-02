"""Test suite for ensuring correct FGS pointing"""
import logging
from numpy import array
from numpy import isclose
from pathlib import Path

import pytest

from stdatamodels.jwst.datamodels import Level1bModel

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
        'exposure': {
            'type': 'FGS_ACQ1',
        },
        'aperture': {
            'name': 'FGS1_FULL',
        },
        'observation': {
            'date': '2017-01-01',
        }
    }
}


@pytest.mark.xfail(reason='See JP-2658')
def test_fgs_pointing():
    model = make_level1b()
    stp.update_wcs(model, siaf_path=siaf_path)

    assert isclose(model.meta.wcsinfo.pc1_1, -0.9997617158628777, atol=1e-15)
    assert isclose(model.meta.wcsinfo.pc1_2, -0.021829143247382235, atol=1e-15)
    assert isclose(model.meta.wcsinfo.pc2_1, -0.021829143247382235, atol=1e-15)
    assert isclose(model.meta.wcsinfo.pc2_2, 0.9997617158628777, atol=1e-15)


# ---------
# Utilities
# ---------
def make_level1b():
    data = array([1.])
    data.shape = (1, 1, 1, 1)
    model = Level1bModel(data)
    model.update(WCS_META)
    return model
