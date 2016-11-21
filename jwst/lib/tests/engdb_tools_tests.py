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

GOOD_MNEMONIC = 'INRSI_GWA_Y_TILT_AVGED'
GOOD_STARTTIME = '2016-01-01'
GOOD_ENDTIME = '2016-12-31'

SHORT_STARTTIME = '2016-01-14'

BAD_SERVER = 'https://www.stsci.edu'
BAD_MNEMONIC = 'No_Such_MNEMONIC'


def test_basic():
    engdb = engdb_tools.ENGDB_Service()
    assert engdb.get_records(GOOD_MNEMONIC, GOOD_STARTTIME, GOOD_ENDTIME)


def test_bad_server():
    engdb = engdb_tools.ENGDB_Service(BAD_SERVER)
    with pytest.raises(Exception):
        records = engdb.get_records(
            GOOD_MNEMONIC, GOOD_STARTTIME, GOOD_ENDTIME
        )


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


def test_values():
    engdb = engdb_tools.ENGDB_Service()
    records = engdb.get_records(GOOD_MNEMONIC, SHORT_STARTTIME, SHORT_STARTTIME)
    assert records['Count'] == 2
    values = engdb.get_values(GOOD_MNEMONIC, SHORT_STARTTIME, SHORT_STARTTIME)
    assert len(values) == 1
    assert values[0] == 0.1968553
