"""Test the siaf db classes"""
from contextlib import nullcontext as does_not_raise
import os
from pathlib import Path
import pytest

from jwst.lib import siafdb

pytest.importorskip('pysiaf')
import pysiaf  # noqa: 402

# Database paths
DATA_PATH = Path(__file__).parent / 'data'
PYSIAF_PRD_PATH = Path(pysiaf.JWST_PRD_DATA_ROOT).parent.parent.parent
XML_DATA_SIAFXML_PATH = DATA_PATH / 'xml_data_siafxml'
SIAFXML_PATH = XML_DATA_SIAFXML_PATH / 'SIAFXML'
OLD_PRD = 'PRDOPSSOC-053'
OLD_PRD_PATH = PYSIAF_PRD_PATH / OLD_PRD / 'SIAFXML' / 'SIAFXML'


@pytest.fixture
def jail_environ():
    """Lock changes to the environment"""
    original = os.environ.copy()
    try:
        yield
    finally:
        os.environ = original


@pytest.mark.parametrize('source, prd, xml_path, exception', [
    # Default
    (None, None, pysiaf.JWST_PRD_DATA_ROOT, does_not_raise()),
    # User-define XML path
    (SIAFXML_PATH, None, SIAFXML_PATH, does_not_raise()),
    # Use $XML_DATA
    ('XML_DATA', None, SIAFXML_PATH, does_not_raise()),
    # Non-existent user-define XML path
    ('junk_source', None, None, pytest.raises(ValueError)),
    # Latest pysiaf PRD version
    (None, pysiaf.JWST_PRD_VERSION, pysiaf.JWST_PRD_DATA_ROOT, does_not_raise()),
    # User-specified PRD version
    (None, OLD_PRD, OLD_PRD_PATH, does_not_raise()),
    # Non-existent PRD version
    (None, 'junk_prd', None, pytest.raises(ValueError)),
    # `source` overrides `prd`
    (SIAFXML_PATH, OLD_PRD, SIAFXML_PATH, does_not_raise()),
])
def test_create(source, prd, xml_path, exception, jail_environ):
    """Test the the right objects are created"""
    source_actual = source
    if source == 'XML_DATA':
        os.environ['XML_DATA'] = str(XML_DATA_SIAFXML_PATH)
        source_actual = None

    with exception:
        siaf_db = siafdb.SiafDb(source_actual, prd)
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
    siaf_db = siafdb.SiafDb(SIAFXML_PATH)
    siaf = siaf_db.get_wcs(aperture)
    assert siaf == expected


@pytest.mark.parametrize('prd, expected, exception', [
    # Valid cases
    ('PRDOPSSOC-055', 'PRDOPSSOC-055', does_not_raise()),
    ('PRDOPSSOC-054', 'PRDOPSSOC-053', does_not_raise()),
    ('PRDOPSSOC-053', 'PRDOPSSOC-053', does_not_raise()),

    # Bad specification
    ('PRDOPSSOC', 'PRDOPSSOC-053', pytest.raises(ValueError)),
    ('junk', 'PRDOPSSOC-053', pytest.raises(ValueError)),

    # Out-of-range specs
    # 999 case should produce whatever the "latest" is. Will need
    # regular updating.
    ('PRDOPSSOC-000', 'PRDOPSSOC-053', pytest.raises(ValueError)),
    ('PRDOPSSOC-999', 'PRDOPSSOC-055', does_not_raise()),
])
def test_nearest_prd(prd, expected, exception):
    """Test nearest prd finding"""

    with exception:
        prd_to_use, _ = siafdb.nearest_prd(pysiaf, prd)
        assert prd_to_use == expected
