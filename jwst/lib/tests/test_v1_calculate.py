"""Test module v1_calculate"""
import sys
from pathlib import Path

import pytest

from astropy.table import Table
from astropy.time import Time
from astropy.utils.diff import report_diff_values

from jwst.datamodels import ImageModel
from jwst.lib import engdb_mast
from jwst.lib import engdb_tools
import jwst.lib.set_telescope_pointing as stp
import jwst.lib.v1_calculate as v1c

from jwst.lib.tests.engdb_mock import EngDB_Mocker

DATA_PATH = Path(__file__).parent / 'data'

# Engineering parameters
# Time range corresponds to simulation producing exposure jw00624028002_02101_00001_nrca1_uncal.fits
# Midpoint is about 2021-01-26T02:32:26.205
GOOD_STARTTIME = Time('59240.10349754328', format='mjd')
GOOD_ENDTIME = Time('59240.1082197338', format='mjd')

# Requires pysiaf
pytest.importorskip('pysiaf')


@pytest.fixture
def engdb_service():
    """Setup the service to operate through the mock service"""
    with EngDB_Mocker():
        yield engdb_tools.ENGDB_Service()


def test_from_models(engdb_service, tmp_path):
    """Test v1_calculate_from_models for basic running"""
    model = ImageModel()
    model.meta.exposure.start_time = GOOD_STARTTIME.mjd
    model.meta.exposure.end_time = GOOD_ENDTIME.mjd

    v1_table = v1c.v1_calculate_from_models([model], method=stp.Methods.COARSE_TR_202107)
    v1_formatted = v1c.simplify_table(v1_table)

    # Save for post-test examination
    v1_formatted.write(tmp_path / 'test_from_models_service.ecsv', format='ascii.ecsv')

    truth = Table.read(DATA_PATH / 'test_from_models_service.ecsv')
    assert report_diff_values(truth, v1_formatted, fileobj=sys.stderr)


def test_over_time(engdb_service, tmp_path):
    """Test v1_calculate_over_time for basic running"""
    v1_table = v1c.v1_calculate_over_time(GOOD_STARTTIME.mjd, GOOD_ENDTIME.mjd, method=stp.Methods.COARSE_TR_202107)
    v1_formatted = v1c.simplify_table(v1_table)

    # Save for post-test examination
    v1_formatted.write(tmp_path / 'test_over_time_service.ecsv', format='ascii.ecsv')

    truth = Table.read(DATA_PATH / 'test_over_time_service.ecsv')
    assert report_diff_values(truth, v1_formatted, fileobj=sys.stderr)


def test_from_models_mast(tmp_path):
    """Test v1_calculate_from_models for basic running"""
    model = ImageModel()
    model.meta.exposure.start_time = GOOD_STARTTIME.mjd
    model.meta.exposure.end_time = GOOD_ENDTIME.mjd

    try:
        v1_table = v1c.v1_calculate_from_models([model], method=stp.Methods.COARSE_TR_202107, engdb_url=engdb_mast.MAST_BASE_URL)
    except ValueError as exception:
        pytest.xfail(f'MAST engineering database not available, possibly no token specified: {exception}')
    v1_formatted = v1c.simplify_table(v1_table)

    # Save for post-test examination
    v1_formatted.write(tmp_path / 'test_from_models_mast.ecsv', format='ascii.ecsv')

    truth = Table.read(DATA_PATH / 'test_from_models_mast.ecsv')
    assert report_diff_values(truth, v1_formatted, fileobj=sys.stderr)


def test_over_time_mast(tmp_path):
    """Test v1_calculate_over_time for basic running"""
    try:
        v1_table = v1c.v1_calculate_over_time(GOOD_STARTTIME.mjd, GOOD_ENDTIME.mjd,
                                              method=stp.Methods.COARSE_TR_202107, engdb_url=engdb_mast.MAST_BASE_URL)
    except ValueError as exception:
        pytest.xfail(f'MAST engineering database not available, possibly no token specified: {exception}')
    v1_formatted = v1c.simplify_table(v1_table)

    # Save for post-test examination
    v1_formatted.write(tmp_path / 'test_over_time_mast.ecsv', format='ascii.ecsv')

    truth = Table.read(DATA_PATH / 'test_over_time_mast.ecsv')
    assert report_diff_values(truth, v1_formatted, fileobj=sys.stderr)
