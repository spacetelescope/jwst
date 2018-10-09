"""Test calwebb_spec3 with MIRI LRS-Fixedslit

Notes
-----
The test data has been arranged to match that of the pandokia jwst_test_data
file structure. As such the environmental variable TEST_BIGDATA points to
the top of the example data tree.
"""

from glob import glob
from os import path

import pytest

from .helpers import (
    SCRIPT_DATA_PATH,
    abspath,
    update_asn_basedir,
)

from ...associations import load_asn
from ...stpipe.step import Step

DATAPATH = abspath(
    '$TEST_BIGDATA/miri/test_datasets/lrs_fixedslit'
)


@pytest.mark.xfail(
    reason='Fails due to issue #947',
    run=False,
)
@pytest.mark.slow
@pytest.mark.bigdata
def test_run_outlier_only(mk_tmp_dirs):
    """Test a basic run"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    asn_path = update_asn_basedir(
        path.join(DATAPATH, 'jw80600-a3002_20170227t160430_spec3_001_asn.json'),
        root=path.join(DATAPATH, 'level2b')
    )
    args = [
        path.join(SCRIPT_DATA_PATH, 'calwebb_spec3_default.cfg'),
        asn_path,
        '--steps.mrs_imatch.skip=true',
        '--steps.resample_spec.skip=true',
        '--steps.cube_build.skip=true',
        '--steps.extract_1d.skip=true',
    ]

    Step.from_cmdline(args)
    assert False


@pytest.mark.bigdata
def test_run_resample_mock_only(mk_tmp_dirs):
    """Test resample step only."""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    asn_path = update_asn_basedir(
        path.join(DATAPATH, 'jw80600-a3002_20170227t160430_spec3_001_asn.json'),
        root=path.join(DATAPATH, 'level2b')
    )
    args = [
        path.join(SCRIPT_DATA_PATH, 'calwebb_spec3_mock.cfg'),
        asn_path,
        '--steps.mrs_imatch.skip=true',
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


@pytest.mark.bigdata
def test_run_extract_1d_resample_mock_only(mk_tmp_dirs):
    """Test resample step only."""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    asn_path = update_asn_basedir(
        path.join(DATAPATH, 'jw80600-a3002_20170227t160430_spec3_001_asn.json'),
        root=path.join(DATAPATH, 'level2b')
    )
    args = [
        path.join(SCRIPT_DATA_PATH, 'calwebb_spec3_mock.cfg'),
        asn_path,
        '--steps.mrs_imatch.skip=true',
        '--steps.outlier_detection.skip=true',
        '--steps.cube_build.skip=true',
    ]

    Step.from_cmdline(args)

    with open(asn_path) as fd:
        asn = load_asn(fd)
    product_name_base = asn['products'][0]['name']
    product_name = product_name_base + '_s2d.fits'
    assert path.isfile(product_name)
    product_name = product_name_base + '_x1d.fits'
    assert path.isfile(product_name)


@pytest.mark.xfail(
    reason='Failure documented in issue #1006',
    run=False,
)
@pytest.mark.slow
@pytest.mark.bigdata
def test_run_resample_only(mk_tmp_dirs):
    """Test resample step only."""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    asn_path = update_asn_basedir(
        path.join(DATAPATH, 'jw80600-a3002_20170227t160430_spec3_001_asn.json'),
        root=path.join(DATAPATH, 'level2b')
    )
    args = [
        path.join(SCRIPT_DATA_PATH, 'calwebb_spec3_default.cfg'),
        asn_path,
        '--steps.mrs_imatch.skip=true',
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


@pytest.mark.bigdata
def test_run_extract_1d_only(mk_tmp_dirs):
    """Test only the extraction step. Should produce nothing
    because extraction requires resampling
    """
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    asn_path = update_asn_basedir(
        path.join(DATAPATH, 'jw80600-a3002_20170227t160430_spec3_001_asn.json'),
        root=path.join(DATAPATH, 'level2b')
    )
    args = [
        path.join(SCRIPT_DATA_PATH, 'calwebb_spec3_default.cfg'),
        asn_path,
        '--steps.mrs_imatch.skip=true',
        '--steps.outlier_detection.skip=true',
        '--steps.resample_spec.skip=true',
        '--steps.cube_build.skip=true',
    ]

    Step.from_cmdline(args)

    # Check that no other products built
    files = glob('*s3d*')
    files.extend(glob('*s2d*'))
    files.extend(glob('*x1d*'))
    assert not files


@pytest.mark.bigdata
def test_run_nosteps(mk_tmp_dirs):
    """Test where no steps execute"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    asn_path = update_asn_basedir(
        path.join(DATAPATH, 'jw80600-a3002_20170227t160430_spec3_001_asn.json'),
        root=path.join(DATAPATH, 'level2b')
    )
    args = [
        path.join(SCRIPT_DATA_PATH, 'calwebb_spec3_default.cfg'),
        asn_path,
        '--steps.mrs_imatch.skip=true',
        '--steps.outlier_detection.skip=true',
        '--steps.resample_spec.skip=true',
        '--steps.cube_build.skip=true',
        '--steps.extract_1d.skip=true',
    ]

    Step.from_cmdline(args)

    # Check that no other products built
    files = glob('*s3d*')
    files.extend(glob('*s2d*'))
    files.extend(glob('*x1d*'))
    assert not files


@pytest.mark.xfail(
    reason='None of the steps operate',
    run=False,
)
@pytest.mark.slow
@pytest.mark.bigdata
def test_run_mir_lrs_fixedslit(mk_tmp_dirs):
    """Test a full run"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    asn_path = update_asn_basedir(
        path.join(DATAPATH, 'jw80600-a3002_20170227t160430_spec3_001_asn.json'),
        root=path.join(DATAPATH, 'level2b')
    )
    args = [
        path.join(SCRIPT_DATA_PATH, 'calwebb_spec3_default.cfg'),
        asn_path,
    ]

    Step.from_cmdline(args)
    assert False
