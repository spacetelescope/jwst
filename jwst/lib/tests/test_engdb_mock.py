"""
Tests of engdb_mock
"""
import os
import pytest
import requests
from tempfile import TemporaryDirectory
import warnings

from astropy.time import Time

from . import engdb_mock
from .. import engdb_tools

GOOD_MNEMONIC = 'inrsi_gwa_y_tilt_avged'
GOOD_STARTTIME = '2016-01-23 00:00:00'
GOOD_ENDTIME = '2016-01-23 00:03:00'
EARLY_STARTTIME = '2014-01-01'
EARLY_ENDTIME = '2014-01-02'
LATE_STARTTIME = '2034-01-01'
LATE_ENDTIME = '2034-01-02'


def test_file_names():
    """
    Test the file names are constructed
    """
    fname = engdb_mock.mnemonic_data_fname(GOOD_MNEMONIC)
    assert fname == GOOD_MNEMONIC.lower() + engdb_mock.DATA


def test_cache_engdb(engdb):
    """
    Test cache creation

    Notes
    -----
    The engdb fixture is called since the cache routine requires
    the live engineering db RESTful service.
    """
    with TemporaryDirectory() as db_path:
        engdb_mock.cache_engdb(
            mnemonics=[GOOD_MNEMONIC],
            starttime=GOOD_STARTTIME,
            endtime=GOOD_ENDTIME,
            db_path=db_path
        )

        assert os.path.isfile(os.path.join(
            db_path, engdb_mock.META
        ))
        assert os.path.isfile(os.path.join(
            db_path,
            engdb_mock.mnemonic_data_fname(GOOD_MNEMONIC)
        ))


@pytest.mark.parametrize(
    "mnemonic",
    [
        GOOD_MNEMONIC,
        'CAL',
        ''
    ]
)
def test_cache_meta(db_cache, engdb, mnemonic):
    """
    Test read of the meta information
    """
    meta = db_cache.fetch_meta(mnemonic)
    live_meta = engdb.get_meta(mnemonic)
    assert meta == live_meta


def test_cache_data(db_cache, engdb):
    """
    Test read of the data information
    """
    data = db_cache.fetch_data(GOOD_MNEMONIC, GOOD_STARTTIME, GOOD_ENDTIME)
    live_data = engdb.get_records(GOOD_MNEMONIC, GOOD_STARTTIME, GOOD_ENDTIME)
    assert data == live_data


def test_cache_partial_data(db_cache, engdb):
    """
    Test read of some data.
    """
    data = db_cache.fetch_data(GOOD_MNEMONIC, GOOD_STARTTIME, GOOD_ENDTIME)
    assert data['Count'] > 4  # Just to make sure we have some data
    endtime = Time(
        engdb_tools.extract_db_time(data['Data'][1]['ObsTime']) / 1000.,
        format='unix'
    )

    data_short = db_cache.fetch_data(
        GOOD_MNEMONIC,
        GOOD_STARTTIME,
        endtime.iso
    )
    live_data_short = engdb.get_records(
        GOOD_MNEMONIC,
        GOOD_STARTTIME,
        endtime.iso
    )
    assert data_short == live_data_short


def test_cache_end_data(db_cache, engdb):
    """
    Test read of some data.
    """
    data = db_cache.fetch_data(GOOD_MNEMONIC, GOOD_STARTTIME, GOOD_ENDTIME)
    assert data['Count'] > 4  # Just to make sure we have some data

    # First test pre-data.
    data_short = db_cache.fetch_data(
        GOOD_MNEMONIC,
        EARLY_STARTTIME,
        EARLY_ENDTIME
    )
    live_data_short = engdb.get_records(
        GOOD_MNEMONIC,
        EARLY_STARTTIME,
        EARLY_ENDTIME
    )
    # We only check counts because the actual data may be different.
    assert data_short['Count'] == live_data_short['Count']

    # Filter ERFA warnings for times far into the future, LATE_STARTTIME and
    # LATE_ENDTIME.  We don't know leap seconds between now and 2034.
    warnings.filterwarnings('ignore', message='ERFA function ')
    # Now for post data
    data_short = db_cache.fetch_data(
        GOOD_MNEMONIC,
        LATE_STARTTIME,
        LATE_ENDTIME
    )
    live_data_short = engdb.get_records(
        GOOD_MNEMONIC,
        LATE_STARTTIME,
        LATE_ENDTIME
    )
    # We only check counts because the actual data may be different.
    assert data_short['Count'] == live_data_short['Count']


def test_mocker_alive(db_cache):
    with engdb_mock.EngDB_Mocker(db_path=db_cache.db_path):
        query = ''.join([
            engdb_tools.ENGDB_BASE_URL,
            engdb_tools.ENGDB_METADATA
        ])
        response = requests.get(query)
        assert response.status_code == 200


@pytest.mark.parametrize(
    'mnemonic, count',
    [
        (GOOD_MNEMONIC, 1),
        ('CAL', 14),
        ('', 2100),
    ]
)
def test_mocker_meta(db_cache, mnemonic, count):
    with engdb_mock.EngDB_Mocker(db_path=db_cache.db_path):
        query = ''.join([
            engdb_tools.ENGDB_BASE_URL,
            engdb_tools.ENGDB_METADATA,
            mnemonic
        ])
        response = requests.get(query)
        assert response.status_code == 200
        meta = response.json()
        assert meta['Count'] >= count


def test_mocker_data(db_cache, engdb):
    with engdb_mock.EngDB_Mocker(db_path=db_cache.db_path):
        query = ''.join([
            engdb_tools.ENGDB_BASE_URL,
            engdb_tools.ENGDB_DATA,
            GOOD_MNEMONIC,
            '?sTime=',
            GOOD_STARTTIME,
            '&eTime=',
            GOOD_ENDTIME
        ])
        response = requests.get(query)

    live_data = engdb.get_records(GOOD_MNEMONIC, GOOD_STARTTIME, GOOD_ENDTIME)
    assert response.json() == live_data


# #########
# Utilities
# #########
@pytest.fixture
def engdb(scope='module'):
    """
    Ensure the live engineering RESTful service is available
    """
    try:
        engdb = engdb_tools.ENGDB_Service()
    except Exception:
        pytest.skip('ENGDB service is not accessible.')
    else:
        return engdb


@pytest.fixture
def db_path():
    """
    Provide a local cache directory
    """
    with TemporaryDirectory() as db_path:
        yield db_path


@pytest.fixture
def db_cache(engdb, db_path):
    """
    Provide a local cache
    """
    engdb_mock.cache_engdb(
        mnemonics=[GOOD_MNEMONIC],
        starttime=GOOD_STARTTIME,
        endtime=GOOD_ENDTIME,
        db_path=db_path
    )

    return engdb_mock.EngDB_Local(db_path=db_path)
