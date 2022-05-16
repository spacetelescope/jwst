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
    'source, expected, use_pysiaf',
    [
        (None, siafdb.SiafDbPySiaf, True),
        (Path(__file__).parent / 'data' / 'siaf.db', siafdb.SiafDbSqlite, False),
        ('XML_DATA', siafdb.SiafDbSqlite, False)
    ]
)
def test_create(source, expected, use_pysiaf, jail_environ):
    """Test the the right objects are created"""
    if use_pysiaf:
        pytest.importorskip('pysiaf')

    if source == 'XML_DATA':
        os.environ['XML_DATA'] = str(Path(__file__).parent / 'data')
        source = None

    with siafdb.SiafDb(source) as siaf_db:
        assert isinstance(siaf_db._db, expected)


@pytest.mark.parametrize(
    'source, expected, use_pysiaf',
    [
        (None,
         siafdb.SIAF(v2_ref=206.464, v3_ref=-697.97, v3yangle=-1.25081713, vparity=1,
                     crpix1=1024.5, crpix2=1024.5, cdelt1=0.06839158, cdelt2=0.06993081,
                     vertices_idl=(-68.8543, 70.1233, 71.5697, -70.2482, 72.1764, 68.8086, -75.5918, -70.7457)),
         True),
        (Path(__file__).parent / 'data' / 'siaf.db',
         siafdb.SIAF(v2_ref=206.464, v3_ref=-697.97, v3yangle=-1.25081713, vparity=1,
                     crpix1=1024.5, crpix2=1024.5, cdelt1=0.06853116, cdelt2=0.07005886,
                     vertices_idl=(-69.01, 70.294, 71.8255, -70.3952, 72.3452, 68.951, -75.8935, -70.8365)),
         False),
    ]
)
def test_get_wcs(source, expected, use_pysiaf):
    """Test retrieval of wcs information."""
    if use_pysiaf:
        pytest.importorskip('pysiaf')

    with siafdb.SiafDb(source) as siaf_db:
        siaf = siaf_db.get_wcs('FGS1_FULL_OSS')
    assert siaf == expected
