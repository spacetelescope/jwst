"""Test calwebb_ami3 with NIR"""

from collections import defaultdict
from glob import glob
from os import path

import pytest

from .helpers import (
    abspath,
    update_asn_basedir,
)

from ...associations import load_asn
from ...stpipe.step import (Step, remove_suffix)
from ..collect_pipeline_cfgs import collect_pipeline_cfgs

DATAPATH = abspath(
    path.join('$TEST_BIGDATA', 'niriss', 'test_ami_pipeline')
)


@pytest.mark.slow
@pytest.mark.bigdata
def test_run_full(mk_tmp_dirs):
    """Test a full run"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    cfg_dir = path.join(tmp_data_path, 'cfgs')
    collect_pipeline_cfgs(cfg_dir)

    asn_path = update_asn_basedir(
        path.join(DATAPATH, 'test_lg1_asn.json'),
        root=DATAPATH
    )
    args = [
        path.join(cfg_dir, 'calwebb_ami3.cfg'),
        asn_path,
    ]

    Step.from_cmdline(args)

    # Now test for file existence. Get the association
    with open(asn_path) as fh:
        asn = load_asn(fh)
    acid = asn['asn_id']
    product = asn['products'][0]
    product_name = product['name']
    members_by_type = defaultdict(list)
    for member in product['members']:
        expname = path.split(member['expname'])[1]
        members_by_type[member['exptype'].lower()].append(expname)

    output_files = glob('*')
    print('Created files ares: {}'.format(output_files))

    # Check Level3 products
    product_name_file = product_name + '_amiavg.fits'
    assert product_name_file in output_files
    output_files.remove(product_name_file)

    product_name_file = product_name + '_psf-amiavg.fits'
    assert product_name_file in output_files
    output_files.remove(product_name_file)

    product_name_file = product_name + '_aminorm.fits'
    assert product_name_file in output_files
    output_files.remove(product_name_file)

    # Check Level2 products
    for member in members_by_type['psf']:
        name, ext = path.splitext(path.split(member)[1])
        name, separator = remove_suffix(name)
        name = name + separator + acid + separator + 'ami' + ext
        assert name in output_files
        output_files.remove(name)

    for member in members_by_type['science']:
        name, ext = path.splitext(path.split(member)[1])
        name, separator = remove_suffix(name)
        name = name + separator + acid + separator + 'ami' + ext
        assert name in output_files
        output_files.remove(name)

    # If there are files left, this is an error
    assert len(output_files) == 0


@pytest.mark.slow
@pytest.mark.bigdata
def test_run_single(mk_tmp_dirs):
    """Test a full run"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    cfg_dir = path.join(tmp_data_path, 'cfgs')
    collect_pipeline_cfgs(cfg_dir)

    asn_path = update_asn_basedir(
        path.join(DATAPATH, 'test_lg1_single_asn.json'),
        root=DATAPATH
    )
    args = [
        path.join(cfg_dir, 'calwebb_ami3.cfg'),
        asn_path,
    ]

    Step.from_cmdline(args)

    # Now test for file existence. Get the association
    with open(asn_path) as fh:
        asn = load_asn(fh)
    acid = asn['asn_id']
    product = asn['products'][0]
    product_name = product['name']
    members_by_type = defaultdict(list)
    for member in product['members']:
        expname = path.split(member['expname'])[1]
        members_by_type[member['exptype'].lower()].append(expname)

    output_files = glob('*')
    print('Created files ares: {}'.format(output_files))

    # Check Level3 products
    product_name_file = product_name + '_amiavg.fits'
    assert product_name_file in output_files
    output_files.remove(product_name_file)

    # Check Level2 products
    for member in members_by_type['psf']:
        name, ext = path.splitext(path.split(member)[1])
        name, separator = remove_suffix(name)
        name = name + separator + acid + separator + 'ami' + ext
        assert name in output_files
        output_files.remove(name)

    for member in members_by_type['science']:
        name, ext = path.splitext(path.split(member)[1])
        name, separator = remove_suffix(name)
        name = name + separator + acid + separator + 'ami' + ext
        assert name in output_files
        output_files.remove(name)

    # If there are files left, this is an error
    assert len(output_files) == 0
