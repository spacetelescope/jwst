"""Test calwebb_spec3"""

from glob import glob
from os import path
from shutil import copy as file_copy

import pytest

from .helpers import (
    SCRIPT_DATA_PATH,
    abspath,
)

from ...associations import load_asn
from ...stpipe.step import Step

DATAPATH = abspath(
    '$TEST_BIGDATA/nirspec/test_datasets/msa/simulated-3nod'
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

    # Copy necessary data to the tmp_data_path
    file_copy(
        path.join(DATAPATH, 'two_member_spec3_asn.json'),
        tmp_data_path
    )
    asn_path = path.join(tmp_data_path, 'two_member_spec3_asn.json')

    with open(asn_path) as fp:
        asn = load_asn(fp)
    for product in asn['products']:
        for member in product['members']:
            file_copy(
                path.join(DATAPATH, 'level2b_twoslit', member['expname']),
                tmp_data_path
            )

    args = [
        path.join(SCRIPT_DATA_PATH, 'calwebb_spec3_default.cfg'),
        asn_path,
        '--steps.mrs_imatch.skip=true',
        '--steps.outlier_detection.save_results=true',
        '--steps.resample_spec.skip=true',
        '--steps.cube_build.skip=true',
        '--steps.extract_1d.skip=true',
    ]

    Step.from_cmdline(args)
    assert False


@pytest.mark.bigdata
def test_run_outlier_only_mock(mk_tmp_dirs):
    """Test a basic run"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    # Copy necessary data to the tmp_data_path
    file_copy(
        path.join(DATAPATH, 'two_member_spec3_asn.json'),
        tmp_data_path
    )
    asn_path = path.join(tmp_data_path, 'two_member_spec3_asn.json')

    with open(asn_path) as fp:
        asn = load_asn(fp)
    for product in asn['products']:
        for member in product['members']:
            file_copy(
                path.join(DATAPATH, 'level2b_twoslit', member['expname']),
                tmp_data_path
            )

    args = [
        path.join(SCRIPT_DATA_PATH, 'calwebb_spec3_mock.cfg'),
        asn_path,
        '--steps.mrs_imatch.skip=true',
        '--steps.resample_spec.skip=true',
        '--steps.cube_build.skip=true',
        '--steps.extract_1d.skip=true',
    ]

    Step.from_cmdline(args)

    # Check for the Source-based cal name.
    with open(asn_path) as fp:
        asn = load_asn(fp)
    product_name_template = asn['products'][0]['name']
    product_name_glob = product_name_template.format(
        source_id='s0000[14]',
    ) + '_cal.fits'
    assert len(glob(product_name_glob)) == 2

    # Check for the outlier resutls
    product_name_glob = product_name_template.format(
        source_id='s0000[14]',
    ) + '_crj.fits'
    assert len(glob(product_name_glob)) == 2


@pytest.mark.xfail(
    reason='Fails as documented in issue #1005',
    run=False,
)
@pytest.mark.slow
@pytest.mark.bigdata
def test_run_resample_only(mk_tmp_dirs):
    """Test resample step only."""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    # Copy necessary data to the tmp_data_path
    file_copy(
        path.join(DATAPATH, 'two_member_spec3_asn.json'),
        tmp_data_path
    )
    asn_path = path.join(tmp_data_path, 'two_member_spec3_asn.json')

    with open(asn_path) as fp:
        asn = load_asn(fp)
    for product in asn['products']:
        for member in product['members']:
            file_copy(
                path.join(DATAPATH, 'level2b_twoslit', member['expname']),
                tmp_data_path
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
    product_name_template = asn['products'][0]['name']
    product_name_glob = product_name_template.format(
        source_id='s0000[14]',
    ) + '_cal.fits'
    assert len(glob(product_name_glob)) == 2

    # Check for resample results
    product_name_glob = product_name_template.format(
        source_id='s0000[14]',
    ) + '_s2d.fits'
    assert len(glob(product_name_glob)) == 2


@pytest.mark.slow
@pytest.mark.bigdata
def test_run_resample_mock_only(mk_tmp_dirs):
    """Test resample step only."""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    # Copy necessary data to the tmp_data_path
    file_copy(
        path.join(DATAPATH, 'two_member_spec3_asn.json'),
        tmp_data_path
    )
    asn_path = path.join(tmp_data_path, 'two_member_spec3_asn.json')

    with open(asn_path) as fp:
        asn = load_asn(fp)
    for product in asn['products']:
        for member in product['members']:
            file_copy(
                path.join(DATAPATH, 'level2b_twoslit', member['expname']),
                tmp_data_path
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
    product_name_template = asn['products'][0]['name']
    product_name_glob = product_name_template.format(
        source_id='s0000[14]',
    ) + '_cal.fits'
    assert len(glob(product_name_glob)) == 2

    # Check for resample results
    product_name_glob = product_name_template.format(
        source_id='s0000[14]',
    ) + '_s2d.fits'
    assert len(glob(product_name_glob)) == 2


@pytest.mark.bigdata
def test_run_cube_build(mk_tmp_dirs):
    """NRS MSA data is not cube data. Nothing should happen"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    # Copy necessary data to the tmp_data_path
    file_copy(
        path.join(DATAPATH, 'two_member_spec3_asn.json'),
        tmp_data_path
    )
    asn_path = path.join(tmp_data_path, 'two_member_spec3_asn.json')

    with open(asn_path) as fp:
        asn = load_asn(fp)
    for product in asn['products']:
        for member in product['members']:
            file_copy(
                path.join(DATAPATH, 'level2b_twoslit', member['expname']),
                tmp_data_path
            )

    args = [
        path.join(SCRIPT_DATA_PATH, 'calwebb_spec3_default.cfg'),
        asn_path,
        '--steps.mrs_imatch.skip=true',
        '--steps.outlier_detection.skip=true',
        '--steps.resample_spec.skip=true',
        '--steps.extract_1d.skip=true',
    ]

    Step.from_cmdline(args)

    # Check for the Source-based cal name.
    with open(asn_path) as fp:
        asn = load_asn(fp)
    product_name_template = asn['products'][0]['name']
    product_name_glob = product_name_template.format(
        source_id='s0000[14]',
    ) + '_cal.fits'
    assert len(glob(product_name_glob)) == 2

    # Assert that no cubes were built.
    cube_files = glob('*s3d*')
    assert not cube_files

@pytest.mark.bigdata
def test_run_extract_1d_only(mk_tmp_dirs):
    """Test only the extraction step. Should produce nothing
    because extraction requires resampling
    """
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    # Copy necessary data to the tmp_data_path
    file_copy(
        path.join(DATAPATH, 'two_member_spec3_asn.json'),
        tmp_data_path
    )
    asn_path = path.join(tmp_data_path, 'two_member_spec3_asn.json')

    with open(asn_path) as fp:
        asn = load_asn(fp)
    for product in asn['products']:
        for member in product['members']:
            file_copy(
                path.join(DATAPATH, 'level2b_twoslit', member['expname']),
                tmp_data_path
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

    # Though the calibration is not run, the conversion to
    # source base has occured. Check
    with open(asn_path) as fd:
        asn = load_asn(fd)
    product_name_template = asn['products'][0]['name']
    product_name_glob = product_name_template.format(
        source_id='s0000[14]',
    ) + '_cal.fits'
    assert len(glob(product_name_glob)) == 2

    # Check that no other products built
    files = glob('*s3d*')
    files.extend(glob('*s2d*'))
    files.extend(glob('*x1d*'))
    assert not files


@pytest.mark.bigdata
def test_run_extract_1d_resample_mock(mk_tmp_dirs):
    """Test only the extraction step. Should produce nothing
    because extraction requires resampling
    """
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    # Copy necessary data to the tmp_data_path
    file_copy(
        path.join(DATAPATH, 'two_member_spec3_asn.json'),
        tmp_data_path
    )
    asn_path = path.join(tmp_data_path, 'two_member_spec3_asn.json')

    with open(asn_path) as fp:
        asn = load_asn(fp)
    for product in asn['products']:
        for member in product['members']:
            file_copy(
                path.join(DATAPATH, 'level2b_twoslit', member['expname']),
                tmp_data_path
            )

    args = [
        path.join(SCRIPT_DATA_PATH, 'calwebb_spec3_mock.cfg'),
        asn_path,
        '--steps.mrs_imatch.skip=true',
        '--steps.outlier_detection.skip=true',
        '--steps.cube_build.skip=true',
    ]

    Step.from_cmdline(args)

    # Though the calibration is not run, the conversion to
    # source base has occured. Check
    with open(asn_path) as fd:
        asn = load_asn(fd)
    product_name_template = asn['products'][0]['name']
    product_name_glob = product_name_template.format(
        source_id='s0000[14]',
    ) + '_cal.fits'
    assert len(glob(product_name_glob)) == 2

    product_name_glob = product_name_template.format(
        source_id='s0000[14]',
    ) + '_s2d.fits'
    assert len(glob(product_name_glob)) == 2

    product_name_glob = product_name_template.format(
        source_id='s0000[14]',
    ) + '_x1d.fits'
    assert len(glob(product_name_glob)) == 2


@pytest.mark.bigdata
def test_run_nosteps(mk_tmp_dirs):
    """Test where no steps execute"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    # Copy necessary data to the tmp_data_path
    file_copy(
        path.join(DATAPATH, 'two_member_spec3_asn.json'),
        tmp_data_path
    )
    asn_path = path.join(tmp_data_path, 'two_member_spec3_asn.json')

    with open(asn_path) as fp:
        asn = load_asn(fp)
    for product in asn['products']:
        for member in product['members']:
            file_copy(
                path.join(DATAPATH, 'level2b_twoslit', member['expname']),
                tmp_data_path
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

    # Check for the Source-based cal name.
    with open(asn_path) as fp:
        asn = load_asn(fp)
    product_name_template = asn['products'][0]['name']
    product_name_glob = product_name_template.format(
        source_id='s0000[14]',
    ) + '_cal.fits'
    assert len(glob(product_name_glob)) == 2

    # Check that no other products built
    files = glob('*s3d*')
    files.extend(glob('*s2d*'))
    files.extend(glob('*x1d*'))
    assert not files

@pytest.mark.xfail(
    reason='Many individual steps have failures',
    run=False,
)
@pytest.mark.slow
@pytest.mark.bigdata
def test_run_full(mk_tmp_dirs):
    """Test a full run"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    # Copy necessary data to the tmp_data_path
    file_copy(
        path.join(DATAPATH, 'two_member_spec3_asn.json'),
        tmp_data_path
    )
    asn_path = path.join(tmp_data_path, 'two_member_spec3_asn.json')

    with open(asn_path) as fp:
        asn = load_asn(fp)
    for product in asn['products']:
        for member in product['members']:
            file_copy(
                path.join(DATAPATH, 'level2b_twoslit', member['expname']),
                tmp_data_path
            )

    args = [
        path.join(SCRIPT_DATA_PATH, 'calwebb_spec3_default.cfg'),
        asn_path,
    ]

    Step.from_cmdline(args)
    assert False
