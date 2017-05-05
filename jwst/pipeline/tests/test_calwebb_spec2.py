"""Test calwebb_spec2"""

from contextlib import contextmanager
from os import path
import pytest

from ...associations.asn_from_list import asn_from_list
from ...associations.lib.rules_level2_base import DMSLevel2bBase
from ...datamodels import open as dm_open
from ..calwebb_spec2 import Spec2Pipeline


def abspath(filepath):
    """Get the absolute file path"""
    return path.abspath(path.expanduser(path.expandvars(filepath)))


DATAPATH = abspath(
    '$TEST_BIGDATA/pipelines'
)
EXPFILE = 'jw00035001001_01101_00001_mirimage_rate.fits'
CALFILE = EXPFILE.replace('_rate', '_cal')

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
        Spec2Pipeline.call(asn_file)
        assert path.isfile(CALFILE)


def test_datamodel(tmpdir):
    model = dm_open(path.join(DATAPATH, EXPFILE))
    with tmpdir.as_cwd():
        Spec2Pipeline.call(model)
        assert path.isfile(CALFILE)


def test_file(tmpdir):
    exppath = path.join(DATAPATH, EXPFILE)
    with tmpdir.as_cwd():
        Spec2Pipeline.call(exppath)
        assert path.isfile(CALFILE)
