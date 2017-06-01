"""Test calwebb_spec3"""

from os import path
import pytest

from .helpers import (
    SCRIPT_PATH,
    abspath,
    mk_tmp_dirs,
    update_asn_basedir,
)

from ...associations import load_asn
from ...stpipe.step import Step

DATAPATH = abspath(
    '$TEST_BIGDATA/miri/test_spec3pipeline'
)

# Skip if the data is not available
pytestmark = pytest.mark.skipif(
    not path.exists(DATAPATH),
    reason='Test data not accessible'
)


def test_run_nosteps(mk_tmp_dirs):
    """Test a basic run"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    asn_path = update_asn_basedir(
        path.join(SCRIPT_PATH, 'basic_spec3_asn.json'),
        root=DATAPATH
    )
    args = [
        path.join(SCRIPT_PATH, 'calwebb_spec3_full.cfg'),
        asn_path,
        '--steps.skymatch.skip=true',
        '--steps.outlier_detection.skip=true',
        '--steps.resample_spec.skip=true',
        '--steps.cube_build.skip=true',
        '--steps.extract_1d.skip=true',
    ]

    Step.from_cmdline(args)


def test_run_extract_1d_only(mk_tmp_dirs):
    """Test a basic run"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    asn_path = update_asn_basedir(
        path.join(SCRIPT_PATH, 'basic_spec3_asn.json'),
        root=DATAPATH
    )
    args = [
        path.join(SCRIPT_PATH, 'calwebb_spec3_full.cfg'),
        asn_path,
        '--steps.skymatch.skip=true',
        '--steps.outlier_detection.skip=true',
        '--steps.resample_spec.skip=true',
        '--steps.cube_build.skip=true',
    ]

    Step.from_cmdline(args)
