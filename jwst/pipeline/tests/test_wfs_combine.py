"""Test wfs_combine"""

from glob import glob
from os import path

import pytest

from .helpers import abspath

from ...associations import load_asn
from ...associations.lib.rules_level3_base import format_product
from ...pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from ...stpipe.step import Step

DATAPATH = abspath(
    path.join('$TEST_BIGDATA', 'nircam', 'test_wfs_combine')
)


@pytest.mark.slow
@pytest.mark.bigdata
def test_run_full(mk_tmp_dirs):
    """Test a full run"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    collect_pipeline_cfgs(tmp_config_path)

    asn_path = path.join(DATAPATH, 'wfs_3sets_asn.json')

    args = [
        path.join(tmp_config_path, 'calwebb_wfs-image3.cfg'),
        asn_path
    ]

    Step.from_cmdline(args)

    # Get the association info.
    with open(asn_path) as fh:
        asn = load_asn(fh)

    # Check for file existence
    output_files = glob('*')
    print('output_files = {}'.format(output_files))

    for product in asn['products']:
        prod_name = product['name']
        prod_name = format_product(prod_name, suffix='wfscmb')
        prod_name += '.fits'
        assert prod_name in output_files
        output_files.remove(prod_name)

    # There should be no more files
    assert len(output_files) == 0
