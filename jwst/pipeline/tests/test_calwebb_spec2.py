"""Test calwebb_spec2"""

import os
from os import path
import pytest
import tempfile

from ...associations.asn_from_list import asn_from_list
from ...associations.lib.rules_level2_base import DMSLevel2bBase
from ...datamodels import open as dm_open
from ..calwebb_spec2 import Spec2Pipeline
from ...stpipe.step import Step


def abspath(filepath):
    """Get the absolute file path"""
    return path.abspath(path.expanduser(path.expandvars(filepath)))


DATAPATH = abspath(
    '$TEST_BIGDATA/pipelines'
)
EXPFILE = 'jw00035001001_01101_00001_mirimage_rate.fits'
CALFILE = EXPFILE.replace('_rate', '_cal')
BSUBFILE = EXPFILE.replace('_rate', '_bsub')

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


def test_asn_with_bkg(mk_tmp_dirs):
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs
    exppath = path.join(DATAPATH, EXPFILE)
    lv2_meta = {
        'program': 'test',
        'target': 'test',
        'asn_pool': 'test',
    }
    asn = asn_from_list([exppath], rule=DMSLevel2bBase, meta=lv2_meta)
    asn['products'][0]['members'].append({
        'expname': exppath, 'exptype': 'BACKGROUND'
    })
    asn_file, serialized = asn.dump()
    with open(asn_file, 'w') as fp:
        fp.write(serialized)

    args = [
        path.join(path.dirname(__file__), 'calwebb_spec2.cfg'),
        asn_file,
        '--steps.bkg_subtract.save_results=true'
    ]

    Step.from_cmdline(args)

    assert path.isfile(CALFILE)
    assert path.isfile(BSUBFILE)


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
        path.join(path.dirname(__file__), 'calwebb_spec2.cfg'),
        asn_file,
    ]

    Step.from_cmdline(args)

    assert path.isfile(CALFILE)


def test_asn_multiple_products(mk_tmp_dirs):
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs
    exppath = path.join(DATAPATH, EXPFILE)
    lv2_meta = {
        'program': 'test',
        'target': 'test',
        'asn_pool': 'test',
    }
    asn = asn_from_list([exppath, exppath], rule=DMSLevel2bBase, meta=lv2_meta)
    asn['products'][0]['name'] = 'product1'
    asn['products'][1]['name'] = 'product2'
    asn_file, serialized = asn.dump()
    with open(asn_file, 'w') as fp:
        fp.write(serialized)

    args = [
        path.join(path.dirname(__file__), 'calwebb_spec2.cfg'),
        asn_file,
    ]

    Step.from_cmdline(args)

    assert path.isfile('product1_cal.fits')
    assert path.isfile('product2_cal.fits')


def test_datamodel(mk_tmp_dirs):
    model = dm_open(path.join(DATAPATH, EXPFILE))
    cfg = path.join(path.dirname(__file__), 'calwebb_spec2_save.cfg')
    Spec2Pipeline.call(model, config_file=cfg)
    assert path.isfile(CALFILE)


def test_file(mk_tmp_dirs):
    exppath = path.join(DATAPATH, EXPFILE)
    cfg = path.join(path.dirname(__file__), 'calwebb_spec2_save.cfg')
    Spec2Pipeline.call(exppath, config_file=cfg)
    assert path.isfile(CALFILE)
