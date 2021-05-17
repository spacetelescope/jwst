"""Test module v1_calculate"""
import sys
from pathlib import Path

import pytest

from astropy.table import Table
from astropy.time import Time
from astropy.utils.diff import report_diff_values

from jwst.datamodels import ImageModel
from jwst.lib import engdb_tools
import jwst.lib.v1_calculate as v1c

from jwst.lib.tests.engdb_mock import EngDB_Mocker

DATA_PATH = Path(__file__).parent / 'data'

# Engineering parameters
GOOD_STARTTIME = '2016-01-18'
GOOD_ENDTIME = '2016-01-19'


@pytest.fixture
def engdb():
    """Setup the service to operate through the mock service"""
    with EngDB_Mocker():
        yield engdb_tools.ENGDB_Service()


def test_from_models(engdb, tmp_path):
    """Test v1_calculate_from_models for basic running"""
    model = ImageModel()
    model.meta.exposure.start_time = Time(GOOD_STARTTIME).mjd
    model.meta.exposure.end_time = Time(GOOD_ENDTIME).mjd

    v1_table = v1c.v1_calculate_from_models([model])
    v1_formatted = v1c.simplify_table(v1_table)

    # Save for post-test examination
    v1_formatted.write(tmp_path / 'test_from_models.ecsv', format='ascii.ecsv')

    truth = Table.read(DATA_PATH / 'test_from_models.ecsv')
    assert report_diff_values(truth, v1_formatted, fileobj=sys.stderr)


def test_over_time(engdb, tmp_path):
    """Test v1_calculate_over_time for basic running"""
    v1_table = v1c.v1_calculate_over_time(
        Time(GOOD_STARTTIME).mjd, Time(GOOD_ENDTIME).mjd
    )
    v1_formatted = v1c.simplify_table(v1_table)

    # Save for post-test examination
    v1_formatted.write(tmp_path / 'test_over_time.ecsv', format='ascii.ecsv')

    truth = Table.read(DATA_PATH / 'test_over_time.ecsv')
    assert report_diff_values(truth, v1_formatted, fileobj=sys.stderr)
