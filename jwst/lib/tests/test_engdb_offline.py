"""Tests for engdb_tools that do not require DB connection."""

import os

import pytest
import requests

from jwst.lib.engdb_direct import extract_db_time
from jwst.lib.engdb_tools import ENGDB_Service

BAD_SERVER = "http://localhost"


def test_environmental_bad(jail_environ):
    os.environ["ENG_BASE_URL"] = BAD_SERVER
    with pytest.raises(requests.exceptions.ConnectionError, match="localhost"):
        ENGDB_Service()


def test_bad_server():
    with pytest.raises(requests.exceptions.ConnectionError):
        ENGDB_Service(BAD_SERVER)


def test_db_time():
    time = 1234567890123
    stime = "".join(["/Date(", str(time), "+1234", ")/"])
    result = extract_db_time(stime)
    assert result == time
