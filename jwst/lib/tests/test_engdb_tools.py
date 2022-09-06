"""
Test suite for engdb_tools

Notes
-----
This file has been named specifically so it is not
automatically found by py.test. This is because, to test,
a connection to the internal engineering service is needed,
which is generally not available.
"""
import os
import pytest
import requests

from astropy.time import Time

from jwst.lib import engdb_direct, engdb_tools
from jwst.lib.tests.engdb_mock import EngDB_Mocker

GOOD_MNEMONIC = 'INRSI_GWA_Y_TILT_AVGED'
GOOD_STARTTIME = '2021-01-25'
GOOD_ENDTIME = '2021-01-26'

SHORT_STARTTIME = '2021-01-26 02:29:02.188'

BAD_SERVER = 'https://www.stsci.edu'
BAD_MNEMONIC = 'No_Such_MNEMONIC'
NODATA_STARTIME = '2014-01-01'
NODATA_ENDTIME = '2014-01-02'

ALTERNATE_HOST = 'https://twjwdmsemwebag.stsci.edu'
ALTERNATE_URL = ALTERNATE_HOST + '/JWDMSEngFqAccSide2/TlmMnemonicDataSrv.svc/'


def is_alive(url):
    """Check if a url is alive

    Parameters
    ----------
    url: str
        The URL to check.

    Returns
    -------
    is_alive: bool
        True if alive
    """
    is_alive = False
    try:
        r = requests.get(url)
        is_alive = (r.status_code == requests.codes.ok)
    except Exception:
        pass
    return is_alive


@pytest.fixture
def engdb():
    """Setup the service to operate through the mock service"""
    with EngDB_Mocker():
        engdb = engdb_tools.ENGDB_Service(base_url='http://localhost')
        yield engdb


def test_environmental(jail_environ):
    os.environ['ENG_BASE_URL'] = ALTERNATE_URL
    try:
        engdb = engdb_tools.ENGDB_Service()
    except Exception:
        pytest.skip('Alternate engineering db not available for test.')
    assert engdb.base_url == ALTERNATE_URL


def test_environmental_bad(jail_environ):
    alternate = 'https://google.com/'
    did_except = False
    os.environ['ENG_BASE_URL'] = alternate
    try:
        engdb = engdb_tools.ENGDB_Service()
    except Exception:
        did_except = True
    assert did_except, 'DB connection falsely created for {}'.format(engdb.base_url)


def test_basic(engdb):
    assert engdb._get_records(GOOD_MNEMONIC, GOOD_STARTTIME, GOOD_ENDTIME)


def test_bad_server():
    with pytest.raises(Exception):
        engdb_tools.ENGDB_Service(BAD_SERVER)


def test_db_time():
    time = 1234567890123
    stime = ''.join([
        '/Date(',
        str(time),
        '+1234',
        ')/'
    ])
    result = engdb_direct.extract_db_time(stime)
    assert result == time


def test_values(engdb):
    records = engdb._get_records(
        GOOD_MNEMONIC, SHORT_STARTTIME, SHORT_STARTTIME
    )
    assert records['Count'] == 2
    values = engdb.get_values(
        GOOD_MNEMONIC, GOOD_STARTTIME, SHORT_STARTTIME
    )
    assert len(values) == 1
    assert values[0] == 0


def test_values_with_bracket(engdb):
    records = engdb._get_records(
        GOOD_MNEMONIC, SHORT_STARTTIME, SHORT_STARTTIME
    )
    assert records['Count'] == 2
    values = engdb.get_values(
        GOOD_MNEMONIC, SHORT_STARTTIME, SHORT_STARTTIME,
        include_bracket_values=True
    )
    assert len(values) == 2
    assert values[1] == 0


def test_values_with_time(engdb):
    values = engdb.get_values(
        GOOD_MNEMONIC, GOOD_STARTTIME, SHORT_STARTTIME,
        include_obstime=True
    )
    assert len(values) >= 1
    assert isinstance(values[0], tuple)
    assert isinstance(values[0].obstime, Time)


def test_novalues(engdb):
    values = engdb.get_values(
        GOOD_MNEMONIC, NODATA_STARTIME, NODATA_ENDTIME)
    assert len(values) == 0


def test_meta(engdb):
    response = engdb.get_meta(GOOD_MNEMONIC)
    assert response['Count'] == 1
    assert response['TlmMnemonics'][0]['TlmMnemonic'] == GOOD_MNEMONIC


def test_unzip(engdb):
    """Test forunzipped versions of content"""
    values = engdb.get_values(
        GOOD_MNEMONIC, SHORT_STARTTIME, SHORT_STARTTIME,
        include_obstime=True,
        zip_results=False
    )
    assert isinstance(values, tuple)
    assert len(values.obstime) == len(values.value)
