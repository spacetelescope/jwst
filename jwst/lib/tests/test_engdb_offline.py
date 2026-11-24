"""Tests for engdb_tools that do not require DB connection."""

import os

import pytest
import requests

from jwst.lib.engdb_tools import ENGDB_Service

BAD_SERVER = "http://localhost"


def test_environmental_bad(jail_environ):
    os.environ["ENG_BASE_URL"] = BAD_SERVER
    with pytest.raises(requests.exceptions.ConnectionError, match="localhost"):
        ENGDB_Service()


def test_bad_server():
    with pytest.raises(requests.exceptions.ConnectionError):
        ENGDB_Service(BAD_SERVER)
