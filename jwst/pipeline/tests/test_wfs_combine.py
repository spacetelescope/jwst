"""Test wfs_combine"""

from glob import glob
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
from ...stpipe.step import Step

DATAPATH = abspath(
    path.join('$TEST_BIGDATA', 'nircam', 'test_wfs_combine')
)


@runslow
@require_bigdata
def test_run_full(mk_tmp_dirs):
    """Test a full run"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    asn_path = update_asn_basedir(
        path.join(DATAPATH, 'wfs_3sets_asn.json'),
        root=DATAPATH
    )
    args = [
        path.join(SCRIPT_DATA_PATH, 'cfgs', 'wfs_combine.cfg'),
        asn_path,
    ]

    Step.from_cmdline(args)

    # Get the association info.
    with open(asn_path) as fh:
        asn = load_asn(fh)

    # Check for file existenc
    output_files = glob('*')
    for product in asn['products']:
        prod_name = path.splitext(product['name'])[0] + '_wfscmb.fits'
        assert prod_name in output_files
        output_files.remove(prod_name)

    # There should be no more files
    assert len(output_files) == 0
