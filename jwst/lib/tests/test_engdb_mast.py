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
QUERY = ('sa_zattest2', '2022-02-02T22:24:58', '2022-02-02T22:24:59')

# Expected return from query
EXPECTED_RESPONSE = ('theTime,MJD,euvalue,sqldataType\r\n'
                     '2022-02-02 22:24:57.797000,59612.9340022801,-0.7914493680,real\r\n'
                     '2022-02-02 22:24:58.053000,59612.9340052431,-0.7914494276,real\r\n'
                     '2022-02-02 22:24:58.309000,59612.9340082060,-0.7914494276,real\r\n'
                     '2022-02-02 22:24:58.565000,59612.9340111690,-0.7914494276,real\r\n'
                     '2022-02-02 22:24:58.821000,59612.9340141319,-0.7914493680,real\r\n'
                     '2022-02-02 22:24:59.077000,59612.9340170949,-0.7914493680,real\r\n')
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


def test_get_records(engdb):
    """Test getting records"""
    records = engdb._get_records(*QUERY)
    assert engdb.response.text == EXPECTED_RESPONSE
    assert report_diff_values(records, EXPECTED_RECORDS)


@pytest.mark.parametrize(
    'pars, expected',
    [
        ({}, [-0.7914494276, -0.7914494276, -0.7914494276, -0.791449368]),
        ({'include_obstime': True},
         [EngDB_Value(obstime=Time('2022-02-02T22:24:58.053', scale='utc', format='isot'), value=-0.7914494276),
          EngDB_Value(obstime=Time('2022-02-02T22:24:58.309', scale='utc', format='isot'), value=-0.7914494276),
          EngDB_Value(obstime=Time('2022-02-02T22:24:58.565', scale='utc', format='isot'), value=-0.7914494276),
          EngDB_Value(obstime=Time('2022-02-02T22:24:58.821', scale='utc', format='isot'), value=-0.791449368)]),
        ({'include_obstime': True, 'zip_results': False}, EngDB_Value(
            obstime=[Time('2022-02-02T22:24:58.053', scale='utc', format='isot'),
                     Time('2022-02-02T22:24:58.309', scale='utc', format='isot'),
                     Time('2022-02-02T22:24:58.565', scale='utc', format='isot'),
                     Time('2022-02-02T22:24:58.821', scale='utc', format='isot')],
            value=[-0.7914494276, -0.7914494276, -0.7914494276, -0.791449368])),
        ({'include_bracket_values': True},
         [-0.791449368, -0.7914494276, -0.7914494276, -0.7914494276, -0.791449368, -0.791449368])
    ])
def test_get_values(engdb, pars, expected):
    values = engdb.get_values(*QUERY, **pars)
    assert str(values) == str(expected)


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
