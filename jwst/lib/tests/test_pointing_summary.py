"""Test module pointing_summary"""
import sys
from pathlib import Path

import pytest

from astropy.table import Table
from astropy.utils.diff import report_diff_values

from stdatamodels.jwst.datamodels import ImageModel

from jwst.lib import engdb_tools
import jwst.lib.pointing_summary as ps

from jwst.lib.tests.engdb_mock import EngDB_Mocker

DATA_PATH = Path(__file__).parent / 'data'

# Engineering parameters
GOOD_STARTTIME = '2016-01-18'
GOOD_ENDTIME = '2016-01-19'


@pytest.fixture
def engdb():
    """Setup the service to operate through the mock service"""
    with EngDB_Mocker():
        yield engdb_tools.ENGDB_Service(base_url='http://localhost')


@pytest.fixture(scope='module')
def data_path(jail):
    """Create data file with needed header parameters"""
    model = ImageModel()

    model.meta.target.ra = 90.75541667
    model.meta.target.dec = -66.56055556
    model.meta.pointing.ra_v1 = 91.08142005
    model.meta.pointing.dec_v1 = -66.60547869
    model.meta.wcsinfo.ra_ref = 90.70377653
    model.meta.wcsinfo.dec_ref = -66.59540224

    model.save('data.fits')
    return model.meta.filename


def test_calc_pointing_deltas(engdb, data_path):
    """Test `calc_pointing_deltas` basic running"""
    truth = ('Delta(target=<SkyCoord (ICRS): (ra, dec) in deg'
             '\n    (90.75541667, -66.56055556)>, v1=<SkyCoord (ICRS): (ra, dec) in deg'
             '\n    (91.08142005, -66.60547869)>, refpoint=<SkyCoord (ICRS): (ra, dec) in deg'
             '\n    (90.70377653, -66.59540224)>, delta_v1=<Angle 0.13712727'
             ' deg>, delta_refpoint=<Angle 0.04044315 deg>)'
             )
    with ImageModel(str(data_path)) as model:
        deltas = ps.calc_pointing_deltas(model)

    assert truth == str(deltas)


def test_calc_deltas(engdb, data_path):
    """Test `calc_deltas` basic running"""
    with ImageModel(data_path) as model:
        deltas = ps.calc_deltas([model])

    truth = Table.read(DATA_PATH / 'calc_deltas_truth.ecsv')

    # round the delta values to a reasonable level
    deltas[0][4] = round(deltas[0][4], 8)
    deltas[0][5] = round(deltas[0][5], 8)
    truth[0][4] = round(truth[0][4], 8)
    truth[0][5] = round(truth[0][5], 8)

    assert report_diff_values(truth, deltas, fileobj=sys.stderr)
