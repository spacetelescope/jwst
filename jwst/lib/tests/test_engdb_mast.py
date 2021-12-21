"""Test the MAST Engineering interface"""
import os
import pytest
import requests

from astropy.table import Table
from astropy.time import Time
from astropy.utils.diff import report_diff_values

from jwst.lib.engdb_lib import EngDB_Value
from jwst.lib import engdb_mast

# Test query
QUERY = ('sa_zattest2', '2021-05-22T00:00:00', '2021-05-22T00:00:01')

# Expected return from query
EXPECTED_RESPONSE = ('theTime,MJD,euvalue,sqldataType\r\n'
                     '2021-05-21 23:59:59.816000,59355.9999978704,0.4641213417,real\r\n'
                     '2021-05-22 00:00:00.072000,59356.0000008333,0.4641213417,real\r\n'
                     '2021-05-22 00:00:00.328000,59356.0000037963,0.4641213417,real\r\n'
                     '2021-05-22 00:00:00.584000,59356.0000067593,0.4641213119,real\r\n'
                     '2021-05-22 00:00:00.840000,59356.0000097222,0.4641213119,real\r\n'
                     '2021-05-22 00:00:01.096000,59356.0000126852,0.4641212821,real\r\n')
EXPECTED_RECORDS = Table.read(EXPECTED_RESPONSE, format='ascii.csv')


@pytest.fixture(scope='module')
def is_alive():
    """Check if the MAST portal is accessible
    """
    is_alive = False
    try:
        r = requests.get(engdb_mast.MAST_BASE_URL)
        is_alive = (r.status_code == requests.codes.ok)
    except Exception:
        pass
    if not is_alive:
        pytest.skip(f'MAST url {engdb_mast.MAST_BASE_URL} not available. Skipping.')


@pytest.fixture(scope='module')
def engdb():
    """Open a connection"""
    try:
        engdb = engdb_mast.EngdbMast(base_url=engdb_mast.MAST_BASE_URL)
    except RuntimeError as exception:
        pytest.skip(f'Live MAST Engineering Service not available: {exception}')
    return engdb


def test_aliveness(is_alive):
    """Check connection creation

    Failure is any failure from instantiation.
    """
    engdb_mast.EngdbMast(base_url=engdb_mast.MAST_BASE_URL, token='dummytoken')


@pytest.mark.xfail(reason='Test needs updating once MAST has actual data')
def test_get_records(engdb):
    """Test getting records"""
    records = engdb._get_records(*QUERY)
    assert engdb.response.text == EXPECTED_RESPONSE
    assert report_diff_values(records, EXPECTED_RECORDS)


@pytest.mark.xfail(reason='Test needs updating once MAST has actual data')
@pytest.mark.parametrize(
    'pars, expected',
    [
        ({}, [0.4641213417, 0.4641213417, 0.4641213119, 0.4641213119]),
        ({'include_obstime': True}, [EngDB_Value(obstime=Time(59356.0000008333, scale='utc', format='mjd'), value=0.4641213417),
                                     EngDB_Value(obstime=Time(59356.0000037963, scale='utc', format='mjd'), value=0.4641213417),
                                     EngDB_Value(obstime=Time(59356.0000067593, scale='utc', format='mjd'), value=0.4641213119),
                                     EngDB_Value(obstime=Time(59356.0000097222, scale='utc', format='mjd'), value=0.4641213119)]),
        ({'include_obstime': True, 'zip_results': False}, EngDB_Value(
            obstime=[Time(59356.0000008333, scale='utc', format='mjd'),
                     Time(59356.0000037963, scale='utc', format='mjd'),
                     Time(59356.0000067593, scale='utc', format='mjd'),
                     Time(59356.0000097222, scale='utc', format='mjd')],
            value=[0.4641213417, 0.4641213417, 0.4641213119, 0.4641213119]
        )),
        ({'include_bracket_values': True}, [0.4641213417, 0.4641213417, 0.4641213417, 0.4641213119, 0.4641213119, 0.4641212821])
    ])
def test_get_values(engdb, pars, expected):
    values = engdb.get_values(*QUERY, **pars)
    assert values == expected


def test_negative_aliveness():
    """Ensure failure occurs with a bad url"""
    with pytest.raises(RuntimeError):
        engdb_mast.EngdbMast(base_url='https://127.0.0.1/_engdb_mast_test', token='dummytoken')


def test_notoken(jail_environ):
    """Check that failure occurs without a token"""
    try:
        del os.environ['MAST_API_TOKEN']
    except KeyError:
        pass

    with pytest.raises(RuntimeError):
        engdb_mast.EngdbMast()
