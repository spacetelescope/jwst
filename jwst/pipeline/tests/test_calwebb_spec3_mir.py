"""Test calwebb_spec3"""

from os import path
import pytest

from .helpers import (
    SCRIPT_PATH,
    SCRIPT_DATA_PATH,
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


@pytest.mark.xfail(
    reason='No available testdata'
)
def test_run_skymatch_only(mk_tmp_dirs):
    """Test a basic run"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    asn_path = update_asn_basedir(
        path.join(SCRIPT_PATH, 'mir_mrs_spec3_asn.json'),
        root=DATAPATH
    )
    args = [
        path.join(SCRIPT_DATA_PATH, 'calwebb_spec3_full.cfg'),
        asn_path,
        '--steps.outlier_detection.skip=true',
        '--steps.resample_spec.skip=true',
        '--steps.cube_build.skip=true',
        '--steps.extract_1d.skip=true',
    ]

    Step.from_cmdline(args)


@pytest.mark.xfail(
    reason='due to bug in resample.utils'
)
def test_run_outlier_only(mk_tmp_dirs):
    """Test a basic run"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    asn_path = update_asn_basedir(
        path.join(DATAPATH, 'mir_lrs_fixedslit_spec3_asn.json'),
        root=DATAPATH
    )
    args = [
        path.join(SCRIPT_DATA_PATH, 'calwebb_spec3_full.cfg'),
        asn_path,
        '--steps.skymatch.skip=true',
        '--steps.resample_spec.skip=true',
        '--steps.cube_build.skip=true',
        '--steps.extract_1d.skip=true',
    ]

    Step.from_cmdline(args)
    assert False


def test_run_resample_only(mk_tmp_dirs):
    """Test resample step only."""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    asn_path = update_asn_basedir(
        path.join(DATAPATH, 'mir_lrs_fixedslit_spec3_asn.json'),
        root=DATAPATH
    )
    args = [
        path.join(SCRIPT_DATA_PATH, 'calwebb_spec3_full.cfg'),
        asn_path,
        '--steps.skymatch.skip=true',
        '--steps.outlier_detection.skip=true',
        '--steps.cube_build.skip=true',
        '--steps.extract_1d.skip=true',
    ]

    Step.from_cmdline(args)

    with open(asn_path) as fd:
        asn = load_asn(fd)
    product_name_base = asn['products'][0]['name']
    product_name = product_name_base + '_s2d.fits'
    assert path.isfile(product_name)


def test_run_extract_1d_only(mk_tmp_dirs):
    """Test only the extraction step. Should produce nothing
    because extraction requires resampling
    """
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    asn_path = update_asn_basedir(
        path.join(DATAPATH, 'mir_lrs_fixedslit_spec3_asn.json'),
        root=DATAPATH
    )
    args = [
        path.join(SCRIPT_DATA_PATH, 'calwebb_spec3_full.cfg'),
        asn_path,
        '--steps.skymatch.skip=true',
        '--steps.outlier_detection.skip=true',
        '--steps.resample_spec.skip=true',
        '--steps.cube_build.skip=true',
    ]

    Step.from_cmdline(args)


def test_run_nosteps(mk_tmp_dirs):
    """Test where no steps execute"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    asn_path = update_asn_basedir(
        path.join(DATAPATH, 'mir_lrs_fixedslit_spec3_asn.json'),
        root=DATAPATH
    )
    args = [
        path.join(SCRIPT_DATA_PATH, 'calwebb_spec3_full.cfg'),
        asn_path,
        '--steps.skymatch.skip=true',
        '--steps.outlier_detection.skip=true',
        '--steps.resample_spec.skip=true',
        '--steps.cube_build.skip=true',
        '--steps.extract_1d.skip=true',
    ]

    Step.from_cmdline(args)


@pytest.mark.xfail(
    reason='None of the steps operate'
)
def test_run_mir_lrs_fixedslit(mk_tmp_dirs):
    """Test a full run"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    asn_path = update_asn_basedir(
        path.join(DATAPATH, 'mir_lrs_fixedslit_spec3_asn.json'),
        root=DATAPATH
    )
    args = [
        path.join(SCRIPT_DATA_PATH, 'calwebb_spec3_full.cfg'),
        asn_path,
    ]

    Step.from_cmdline(args)
    assert False
