"""
Test suite for engdb_tools

Notes
-----
This file has been named specifically so it is not
automatically found by py.test. This is because, to test,
a connection to the internal engineering service is needed,
which is generally not available.
"""

from __future__ import absolute_import

from astropy.time import Time
import pytest
import requests_mock

from .. import engdb_tools

GOOD_MNEMONIC = 'INRSI_GWA_Y_TILT_AVGED'
GOOD_STARTTIME = '2016-01-01'
GOOD_ENDTIME = '2016-12-31'

SHORT_STARTTIME = '2016-01-14'

BAD_SERVER = 'https://www.stsci.edu'
BAD_MNEMONIC = 'No_Such_MNEMONIC'


@pytest.fixture
def engdb():
    try:
        engdb = engdb_tools.ENGDB_Service()
    except:
        pytest.skip('ENGDB service is not accessible.')
    else:
        return engdb


@pytest.fixture
def engdb_mock():
    with requests_mock.Mocker() as rm:

        # Define response for aliveness
        url = ''.join([
            engdb_tools.ENGDB_BASE_URL,
            engdb_tools.ENGDB_METADATA
        ])
        rm.get(url, text='Success')

        # Define response for a mneunonic query
        url = ''.join([
            engdb_tools.ENGDB_BASE_URL,
            'Data/',
            GOOD_MNEMONIC,
            '?sTime=',
            Time(SHORT_STARTTIME).iso,
            '&eTime=',
            Time(SHORT_STARTTIME).iso
        ])
        response = {
            "TlmMnemonic": "INRSI_GWA_Y_TILT_AVGED",
            "AllPoints": 1,
            "ReqSTime": "/Date(1452729600000+0000)/",
            "ReqETime": "/Date(1452729600000+0000)/",
            "Count": 2,
            "Data": [
                {
                    "ObsTime": "/Date(1452729600000+0000)/",
                    "EUValue": 0.1968553
                }, {
                    "ObsTime": "/Date(1452731400000+0000)/",
                    "EUValue": 0.1968553
                }
            ]
        }
        rm.get(url, json=response)

        yield rm


def test_basic(engdb):
    assert engdb.get_records(GOOD_MNEMONIC, GOOD_STARTTIME, GOOD_ENDTIME)


def test_bad_server():
    with pytest.raises(Exception):
        engdb = engdb_tools.ENGDB_Service(BAD_SERVER)


def test_db_time():
    time = 1234567890123
    stime = ''.join([
        '/Date(',
        str(time),
        '+1234',
        ')/'
    ])
    result = engdb_tools.extract_db_time(stime)
    assert result == time


def test_values(engdb):
    records = engdb.get_records(
        GOOD_MNEMONIC, SHORT_STARTTIME, SHORT_STARTTIME
    )
    assert records['Count'] == 2
    values = engdb.get_values(
        GOOD_MNEMONIC, SHORT_STARTTIME, SHORT_STARTTIME
    )
    assert len(values) == 1
    assert values[0] == 0.1968553


def test_values_with_time(engdb):
    values = engdb.get_values(
        GOOD_MNEMONIC, SHORT_STARTTIME, SHORT_STARTTIME,
        include_obstime=True
    )
    assert len(values) >= 1
    assert isinstance(values[0], tuple)


def test_values_mock(engdb_mock, engdb):
    records = engdb.get_records(
        GOOD_MNEMONIC, SHORT_STARTTIME, SHORT_STARTTIME
    )
    assert records['Count'] == 2
    values = engdb.get_values(
        GOOD_MNEMONIC, SHORT_STARTTIME, SHORT_STARTTIME
    )
    assert len(values) == 1
    assert values[0] == 0.1968553
