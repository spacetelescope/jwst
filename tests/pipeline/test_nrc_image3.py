"""Test calwebb_image3 with NRC"""

from os import path
import re

from .helpers import (
    SCRIPT_PATH,
    SCRIPT_DATA_PATH,
    abspath,
    mk_tmp_dirs,
    require_bigdata,
    runslow,
    update_asn_basedir,
)

from ...associations import load_asn
from ...stpipe.step import (Step, remove_suffix)
from ..collect_pipeline_cfgs import collect_pipeline_cfgs

DATAPATH = abspath(
    path.join('$TEST_BIGDATA', 'pipelines', 'nircam_calimage3')
)


@runslow
@require_bigdata
def test_run_full(mk_tmp_dirs):
    """Test a full run"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    cfg_dir = path.join(tmp_data_path, 'cfgs')
    collect_pipeline_cfgs(cfg_dir)

    asn_path = path.join(DATAPATH, 'mosaic_long_asn.json')
    args = [
        path.join(cfg_dir, 'calwebb_image3.cfg'),
        asn_path,
    ]

    Step.from_cmdline(args)

    # Check for the CRF files
    with open(asn_path) as fh:
        asn = load_asn(fh)
    expfilenames = [
        path.split(path.splitext(member['expname'])[0])[1]
        for member in asn['products'][0]['members']
    ]
    crffilenames = []
    for expfilename in expfilenames:
        name = remove_suffix(path.splitext(expfilename)[0])[0]
        crffilenames.append(name + '_a3001_crf.fits')
    for crffilename in crffilenames:
        assert path.isfile(crffilename)

    # Check for the level3 products
    product_name = asn['products'][0]['name']
    assert path.isfile(product_name + '_cat.ecsv')
    assert path.isfile(product_name + '_i2d.fits')


@runslow
@require_bigdata
def test_single_image(mk_tmp_dirs):
    """Test a full run"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    cfg_dir = path.join(tmp_data_path, 'cfgs')
    collect_pipeline_cfgs(cfg_dir)

    input_file = 'nrca5_47Tuc_subpix_dither1_newpos_cal.fits'
    args = [
        path.join(cfg_dir, 'calwebb_image3.cfg'),
        path.join(
            DATAPATH,
            input_file
        )
    ]

    Step.from_cmdline(args)

    # Check for the level3 products
    name = remove_suffix(path.splitext(input_file)[0])[0]
    assert path.isfile(name + '_cat.ecsv')
    assert path.isfile(name + '_i2d.fits')
