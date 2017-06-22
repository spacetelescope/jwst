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
    SCRIPT_DATA_PATH,
    abspath,
    mk_tmp_dirs,
    update_asn_basedir,
)

from ...associations import load_asn
from ...stpipe.step import Step

DATAPATH = abspath(
    '$TEST_BIGDATA/miri/test_datasets/mrs/simulated'
)

# Skip if the data is not available
pytestmark = pytest.mark.skipif(
    not path.exists(DATAPATH),
    reason='Test data not accessible'
)


def test_run_extract_1d_only(mk_tmp_dirs):
    """Test only the extraction step. Should produce nothing
    because extraction requires resampling
    """
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    asn_path = update_asn_basedir(
        path.join(DATAPATH, 'single_spec3_asn.json'),
        root=path.join(DATAPATH, 'level2b')
    )
    args = [
        path.join(SCRIPT_DATA_PATH, 'calwebb_spec3_default.cfg'),
        asn_path,
        '--steps.skymatch.skip=true',
        '--steps.outlier_detection.skip=true',
    ]

    Step.from_cmdline(args)

    with open(asn_path) as fd:
        asn = load_asn(fd)
    product_name_base = asn['products'][0]['name']
    product_name = product_name_base + '_ch3-4-long_s3d.fits'
    assert path.isfile(product_name)
    product_name = product_name_base + '_x1d.fits'
    assert path.isfile(product_name)


def test_run_resample_only(mk_tmp_dirs):
    """Test resample step only."""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    asn_path = update_asn_basedir(
        path.join(DATAPATH, 'single_spec3_asn.json'),
        root=path.join(DATAPATH, 'level2b')
    )
    args = [
        path.join(SCRIPT_DATA_PATH, 'calwebb_spec3_default.cfg'),
        asn_path,
        '--steps.skymatch.skip=true',
        '--steps.outlier_detection.skip=true',
        '--steps.extract_1d.skip=true',
    ]

    Step.from_cmdline(args)

    with open(asn_path) as fd:
        asn = load_asn(fd)
    product_name_base = asn['products'][0]['name']
    product_name = product_name_base + '_ch3-4-long_s3d.fits'
    assert path.isfile(product_name)


@pytest.mark.xfail(
    reason='outlier step not ready'
)
def test_run_skymatch_only(mk_tmp_dirs):
    """Test a basic run"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    asn_path = update_asn_basedir(
        path.join(DATAPATH, 'single_spec3_asn.json'),
        root=path.join(DATAPATH, 'level2b')
    )
    args = [
        path.join(SCRIPT_DATA_PATH, 'calwebb_spec3_default.cfg'),
        asn_path,
        '--steps.outlier_detection.skip=true',
        '--steps.resample_spec.skip=true',
        '--steps.cube_build.skip=true',
        '--steps.extract_1d.skip=true',
    ]

    Step.from_cmdline(args)
