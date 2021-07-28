"""Test the siaf db classes"""
import os
from pathlib import Path
import pytest

from jwst.lib import siafdb


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
        (None, siafdb.SiafDbPySiaf),
        (Path(__file__).parent / 'data' / 'siaf.db', siafdb.SiafDbSqlite),
        ('XML_DATA', siafdb.SiafDbSqlite)
    ]
)
def test_create(source, expected, jail_environ):
    """Test the the right objects are created"""
    if source == 'XML_DATA':
        os.environ['XML_DATA'] = Path(__file__).parent / 'data'
        source = None

    siaf = siafdb.SiafDb(source)

    assert isinstance(siaf._source, expected)


@pytest.mark.parametrize(
    'source, expected',
    [
        (None,
         siafdb.SIAF(v2_ref=206.464, v3_ref=-697.97, v3yangle=-1.25081713, vparity=1,
                     crpix1=1024.5, crpix2=1024.5, cdelt1=0.06853116, cdelt2=0.07005886,
                     vertices_idl=(-69.01, 70.294, 71.8255, -70.3952, 72.3452, 68.951, -75.8935, -70.8365))),
        (Path(__file__).parent / 'data' / 'siaf.db',
         siafdb.SIAF(v2_ref=206.464, v3_ref=-697.97, v3yangle=-1.25081713, vparity=1,
                     crpix1=1024.5, crpix2=1024.5, cdelt1=0.06853116, cdelt2=0.07005886,
                     vertices_idl=(-69.01, 70.294, 71.8255, -70.3952, 72.3452, 68.951, -75.8935, -70.8365))),
    ]
)
def test_get_wcs(source, expected):
    """Test retrieval of wcs information."""
    siaf_db = siafdb.SiafDb(source)
    siaf = siaf_db.get_wcs('FGS1_FULL_OSS')
    assert siaf == expected
