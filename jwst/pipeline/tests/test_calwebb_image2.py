"""Test calwebb_image2"""

import os
from os import path
import pytest
import tempfile

from ...associations.asn_from_list import asn_from_list
from ...associations.lib.rules_level2_base import DMSLevel2bBase
from ...datamodels import open as dm_open
from ...stpipe.step import Step
from ..calwebb_image2 import Image2Pipeline


def abspath(filepath):
    """Get the absolute file path"""
    return path.abspath(path.expanduser(path.expandvars(filepath)))


DATAPATH = abspath(
    '${TEST_BIGDATA}/miri/test_image2pipeline/'
)
EXPFILE = 'jw00001001001_01101_00001_MIRIMAGE_uncal_MiriSloperPipeline.fits'
CALFILE = 'jw00001001001_01101_00001_mirimage_uncal_mirisloperpipeline_cal.fits'

# Skip if the data is not available
pytestmark = pytest.mark.skipif(
    not path.exists(DATAPATH),
    reason='Test data not accessible'
)


@pytest.fixture
def mk_tmp_dirs():
    tmp_current_path = tempfile.mkdtemp()
    tmp_data_path = tempfile.mkdtemp()
    tmp_config_path = tempfile.mkdtemp()

    old_path = os.getcwd()
    try:
        os.chdir(tmp_current_path)
        yield (tmp_current_path, tmp_data_path, tmp_config_path)
    finally:
        os.chdir(old_path)


def test_asn(mk_tmp_dirs):
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs
    exppath = path.join(DATAPATH, EXPFILE)
    lv2_meta = {
        'program': 'test',
        'target': 'test',
        'asn_pool': 'test',
    }
    asn = asn_from_list([exppath], rule=DMSLevel2bBase, meta=lv2_meta)
    asn_file, serialized = asn.dump()
    with open(asn_file, 'w') as fp:
        fp.write(serialized)

    args = [
        path.join(path.dirname(__file__), 'calwebb_image2_save.cfg'),
        asn_file,
    ]
    Step.from_cmdline(args)

    assert path.isfile(CALFILE)


def test_datamodel(tmpdir):
    model = dm_open(path.join(DATAPATH, EXPFILE))
    cfg = path.join(path.dirname(__file__), 'calwebb_image2_save.cfg')
    Image2Pipeline.call(model, config_file=cfg)
    assert path.isfile('jw00001001001_01101_00001_MIRIMAGE_uncal_cal.fits')


def test_file(tmpdir):
    exppath = path.join(DATAPATH, EXPFILE)
    cfg = path.join(path.dirname(__file__), 'calwebb_image2_save.cfg')
    Image2Pipeline.call(exppath, config_file=cfg)
    assert path.isfile(CALFILE)
