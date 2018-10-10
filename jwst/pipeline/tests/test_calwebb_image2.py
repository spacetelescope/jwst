"""Test calwebb_image2"""
from os import path

import pytest

from .helpers import (
    SCRIPT_DATA_PATH,
    abspath,
)

from ...associations.asn_from_list import asn_from_list
from ...associations.lib.rules_level2_base import DMSLevel2bBase
from ...datamodels import DataModel
from ...datamodels import open as dm_open
from ...stpipe.step import Step
from ..calwebb_image2 import Image2Pipeline

DATAPATH = abspath(
    '${TEST_BIGDATA}/miri/test_image2pipeline/'
)
EXPFILE = 'jw00001001001_01101_00001_mirimage_rate.fits'
CALFILE = 'jw00001001001_01101_00001_mirimage_cal.fits'

# Skip if the data is not available
pytestmark = pytest.mark.skipif(
    not path.exists(DATAPATH),
    reason='Test data not accessible'
)


@pytest.mark.bigdata
def test_result_return(mk_tmp_dirs):
    """Ensure that a result is returned programmatically"""
    exppath = path.join(DATAPATH, EXPFILE)
    cfg = path.join(SCRIPT_DATA_PATH, 'calwebb_image2_save.cfg')
    results = Image2Pipeline.call(exppath, config_file=cfg)
    assert isinstance(results[0], DataModel)


@pytest.mark.bigdata
def test_no_cfg(mk_tmp_dirs):
    """What happens when the pipeline is run without a config"""
    exppath = path.join(DATAPATH, EXPFILE)
    args = [
        path.join(SCRIPT_DATA_PATH, 'calwebb_image2_empty.cfg'),
        exppath,
    ]
    Step.from_cmdline(args)

    assert path.isfile(CALFILE)


@pytest.mark.bigdata
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
        path.join(SCRIPT_DATA_PATH, 'calwebb_image2_save.cfg'),
        asn_file,
    ]
    Step.from_cmdline(args)

    assert path.isfile(CALFILE)


@pytest.mark.bigdata
def test_datamodel(mk_tmp_dirs):
    model = dm_open(path.join(DATAPATH, EXPFILE))
    cfg = path.join(SCRIPT_DATA_PATH, 'calwebb_image2_save.cfg')
    Image2Pipeline.call(model, config_file=cfg)
    assert path.isfile(CALFILE)


@pytest.mark.bigdata
def test_file(mk_tmp_dirs):
    exppath = path.join(DATAPATH, EXPFILE)
    cfg = path.join(SCRIPT_DATA_PATH, 'calwebb_image2_save.cfg')
    Image2Pipeline.call(exppath, config_file=cfg)
    assert path.isfile(CALFILE)


@pytest.mark.bigdata
def test_file_outputdir(mk_tmp_dirs):
    """Test putting results in another folder"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    exppath = path.join(DATAPATH, EXPFILE)
    outfile = 'junk.fits'

    args = [
        path.join(SCRIPT_DATA_PATH, 'calwebb_image2_save.cfg'),
        exppath,
        '--output_file=' + outfile,
        '--output_dir=' + tmp_data_path,
    ]

    Step.from_cmdline(args)

    result_name, result_ext = path.splitext(outfile)
    result_path = path.join(tmp_data_path, result_name + '_cal.fits')
    assert path.isfile(result_path)
