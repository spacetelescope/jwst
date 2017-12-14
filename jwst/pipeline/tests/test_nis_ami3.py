"""Test calwebb_ami3 with NIR"""

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
from ...stpipe.step import (REMOVE_SUFFIX, Step)

DATAPATH = abspath(
    path.join('$TEST_BIGDATA', 'niriss', 'test_ami_pipeline')
)


@runslow
@require_bigdata
def test_run_full(mk_tmp_dirs):
    """Test a full run"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    asn_path = update_asn_basedir(
        path.join(DATAPATH, 'test_lg1_asn.json'),
        root=DATAPATH
    )
    args = [
        path.join(SCRIPT_DATA_PATH, 'cfgs', 'calwebb_ami3.cfg'),
        asn_path,
    ]

    Step.from_cmdline(args)

    assert False
