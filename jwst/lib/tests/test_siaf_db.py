"""Test the siaf db classes"""
import os
from pathlib import Path
import pytest

from jwst.lib import siaf_db


@pytest.fixture
def jail_environ():
    """Lock changes to the environment"""
    original = os.environ.copy()
    try:
        yield
    finally:
        os.environ = original


@pytest.mark.parametrize(
    'source, expected',
    [
        (None, siaf_db.SiafDbPySiaf),
        (Path(__file__).parent / 'data' / 'siaf.db', siaf_db.SiafDbSqlite),
        ('XML_DATA', siaf_db.SiafDbSqlite)
    ]
)
def test_create(source, expected, jail_environ):
    """Test the the right objects are created"""
    if source == 'XML_DATA':
        os.environ['XML_DATA'] = Path(__file__).parent / 'data'
        source = None

    siaf = siaf_db.SiafDb(source)

    assert isinstance(siaf._source, expected)
