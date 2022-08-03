"""Test the siaf db classes"""
import os
from pathlib import Path
import pytest

from jwst.lib import siafdb

pytest.importorskip('pysiaf')
import pysiaf

# Database paths
DATA_PATH = Path(__file__).parent / 'data'
XML_DATA_SIAFXML_PATH = DATA_PATH / 'xml_data_siafxml'
SIAFXML_PATH = XML_DATA_SIAFXML_PATH / 'SIAFXML'


@pytest.fixture
def jail_environ():
    """Lock changes to the environment"""
    original = os.environ.copy()
    try:
        yield
    finally:
        os.environ = original


@pytest.mark.parametrize('source, xml_path',
                         [(None, pysiaf.JWST_PRD_DATA_ROOT),
                          (SIAFXML_PATH, SIAFXML_PATH),
                          ('XML_DATA', SIAFXML_PATH)
                          ])
def test_create(source, xml_path, jail_environ):
    """Test the the right objects are created"""
    source_actual = source
    if source == 'XML_DATA':
        os.environ['XML_DATA'] = str(XML_DATA_SIAFXML_PATH)
        source_actual = None

    siaf_db = siafdb.SiafDb(source_actual)
    assert str(siaf_db.xml_path) == str(xml_path)


@pytest.mark.parametrize(
    'aperture, expected',
    [
        ('FGS1_FULL_OSS',
         siafdb.SIAF(v2_ref=206.407, v3_ref=-697.765, v3yangle=-1.24120427, vparity=1,
                     crpix1=1024.5, crpix2=1024.5, cdelt1=0.06839158, cdelt2=0.06993081,
                     vertices_idl=(-68.8543, 70.1233, 71.5697, -70.2482, 72.1764, 68.8086, -75.5918, -70.7457))),
        ('FGS2_FULL_OSS',
         siafdb.SIAF(v2_ref=22.835, v3_ref=-699.423, v3yangle=0.2914828, vparity=1,
                     crpix1=1024.5, crpix2=1024.5, cdelt1=0.06787747, cdelt2=0.06976441,
                     vertices_idl=(-70.7843, 69.9807, 68.6042, -69.4153, -74.3558, -71.5516, 71.3065, 69.3639))),
    ]
)
def test_get_wcs(aperture, expected):
    """Test retrieval of wcs information."""
    siaf_db =  siafdb.SiafDb(SIAFXML_PATH)
    siaf = siaf_db.get_wcs(aperture)
    assert siaf == expected
