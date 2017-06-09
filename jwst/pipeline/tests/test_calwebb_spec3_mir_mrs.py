"""Test calwebb_spec3 with MIRI MRS

Notes
-----
The test data has been arranged to match that of the pandokia jwst_test_data
file structure. As such the environmental variable TEST_BIGDATA points to
the top of the example data tree.
"""

from os import path
import pytest

from .helpers import (
    abspath,
    mk_tmp_dirs,
    update_asn_basedir,
)

from ...stpipe.step import Step

DATAPATH = abspath(
    '$TEST_BIGDATA/miri/test_datasets/mrs/simulated'
)

# Skip if the data is not available
pytestmark = pytest.mark.skipif(
    not path.exists(DATAPATH),
    reason='Test data not accessible'
)


@pytest.mark.xfail(
    reason='step not prepared yet'
)
def test_run_skymatch_only(mk_tmp_dirs):
    """Test a basic run"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    asn_path = update_asn_basedir(
        path.join(DATAPATH, 'single_spec3_asn.json'),
        root=path.join(DATAPATH, 'level2b')
    )
    args = [
        path.join(DATAPATH, 'calwebb_spec3_full.cfg'),
        asn_path,
        '--steps.outlier_detection.skip=true',
        '--steps.resample_spec.skip=true',
        '--steps.cube_build.skip=true',
        '--steps.extract_1d.skip=true',
    ]

    Step.from_cmdline(args)
