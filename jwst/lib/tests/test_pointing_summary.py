"""Test module pointing_summary"""
import sys
from pathlib import Path

import pytest

from astropy.table import Table
from astropy.time import Time
from astropy.utils.diff import report_diff_values

import jwst.datamodels as dm
from jwst.lib import engdb_tools
import jwst.lib.pointing_summary as ps

from .engdb_mock import EngDB_Mocker

DATA_PATH = Path(__file__).parent / 'data'

# Engineering parameters
GOOD_STARTTIME = '2016-01-18'
GOOD_ENDTIME = '2016-01-19'


@pytest.fixture
def engdb():
    """Setup the service to operate through the mock service"""
    with EngDB_Mocker():
        engdb = engdb_tools.ENGDB_Service()
        yield engdb


def test_calc_pointing_deltas(engdb):
    """Test `calc_pointing_deltas` basic running"""
    truth = ('Delta(target=<SkyCoord (ICRS): (ra, dec) in deg'
             '\n    (90.75541667, -66.56055556)>, v1=<SkyCoord (ICRS): (ra, dec) in deg'
             '\n    (91.08142005, -66.60547869)>, refpoint=<SkyCoord (ICRS): (ra, dec) in deg'
             '\n    (90.70377653, -66.59540224)>, delta_v1=<Angle 0.13712727'
             ' deg>, delta_refpoint=<Angle 0.04044315 deg>)'
    )
    model = dm.ImageModel(str(DATA_PATH / 'empty.fits'))
    deltas = ps.calc_pointing_deltas(model)

    assert truth == str(deltas)
