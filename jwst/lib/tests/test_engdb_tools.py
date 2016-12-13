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

import pytest

from .. import engdb_tools
from .engdb_mock import EngDB_Mocker

GOOD_MNEMONIC = 'INRSI_GWA_Y_TILT_AVGED'
GOOD_STARTTIME = '2016-01-18'
GOOD_ENDTIME = '2016-01-19'

SHORT_STARTTIME = '2016-01-18 15:30:00'

BAD_SERVER = 'https://www.stsci.edu'
BAD_MNEMONIC = 'No_Such_MNEMONIC'


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
    assert records['Count'] == 3
    values = engdb.get_values(
        GOOD_MNEMONIC, SHORT_STARTTIME, SHORT_STARTTIME
    )
    assert len(values) == 1
    assert values[0] == 0.19687812


def test_values_with_time(engdb):
    values = engdb.get_values(
        GOOD_MNEMONIC, SHORT_STARTTIME, SHORT_STARTTIME,
        include_obstime=True
    )
    assert len(values) >= 1
    assert isinstance(values[0], tuple)


def test_meta(engdb):
    response = engdb.get_meta(GOOD_MNEMONIC)
    assert response['Count'] == 1
    assert response['TlmMnemonics'][0]['TlmMnemonic'] == GOOD_MNEMONIC


# #####################
# Utilities for testing
# #####################
@pytest.fixture
def engdb():
    with EngDB_Mocker() as mocker:
        engdb = engdb_tools.ENGDB_Service()
        yield engdb
