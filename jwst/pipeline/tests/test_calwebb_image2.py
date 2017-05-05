"""Test calwebb_image2"""

from os import path
import pytest

from ...associations.asn_from_list import asn_from_list
from ...associations.lib.rules_level2_base import DMSLevel2bBase
from ...datamodels import open as dm_open
from ..calwebb_image2 import Image2Pipeline


def abspath(filepath):
    """Get the absolute file path"""
    return path.abspath(path.expanduser(path.expandvars(filepath)))


DATAPATH = abspath(
    '${TEST_BIGDATA}/miri/test_image2pipeline/'
)
EXPFILE = 'jw00001001001_01101_00001_MIRIMAGE_uncal_MiriSloperPipeline.fits'
CALFILE = 'cal.fits'

# Skip if the data is not available
pytestmark = pytest.mark.skipif(
    not path.exists(DATAPATH),
    reason='Test data not accessible'
)


def test_asn(tmpdir):
    exppath = path.join(DATAPATH, EXPFILE)
    lv2_meta = {
        'program': 'test',
        'target': 'test',
        'asn_pool': 'test',
    }
    asn = asn_from_list([exppath], rule=DMSLevel2bBase, meta=lv2_meta)
    asn_file, serialized = asn.dump()
    with tmpdir.join(asn_file).open('w') as fp:
        fp.write(serialized)
    with tmpdir.as_cwd():
        Image2Pipeline.call(asn_file, output_file=CALFILE)
        assert path.isfile(CALFILE)


def test_datamodel(tmpdir):
    model = dm_open(path.join(DATAPATH, EXPFILE))
    with tmpdir.as_cwd():
        Image2Pipeline.call(model, output_file=CALFILE)
        assert path.isfile(CALFILE)


def test_file(tmpdir):
    exppath = path.join(DATAPATH, EXPFILE)
    with tmpdir.as_cwd():
        Image2Pipeline.call(exppath, output_file=CALFILE)
        assert path.isfile(CALFILE)
