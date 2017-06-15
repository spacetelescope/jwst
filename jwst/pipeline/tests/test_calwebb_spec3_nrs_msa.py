"""Test calwebb_spec3"""

from os import path
import pytest

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
    '$TEST_BIGDATA/nirspec/test_datasets/msa/simulated-3nod'
)


@runslow
@require_bigdata
def test_run_outlier_only(mk_tmp_dirs):
    """Test a basic run"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    asn_path = update_asn_basedir(
        path.join(DATAPATH, 'two_member_spec3_asn.json'),
        root=path.join(DATAPATH, 'level2b')
    )
    args = [
        path.join(SCRIPT_DATA_PATH, 'calwebb_spec3_default.cfg'),
        asn_path,
        '--steps.skymatch.skip=true',
        '--steps.resample_spec.skip=true',
        '--steps.cube_build.skip=true',
        '--steps.extract_1d.skip=true',
    ]

    Step.from_cmdline(args)
    assert False


@pytest.mark.skip(
    reason='Dies with crds error no META.INSTRUMENT.NAME'
)
@runslow
@require_bigdata
def test_run_resample_only(mk_tmp_dirs):
    """Test resample step only."""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    asn_path = update_asn_basedir(
        path.join(DATAPATH, 'mos_udf_g235M_spec3_asn.json'),
        root=path.join(DATAPATH, 'level2b')
    )
    args = [
        path.join(SCRIPT_DATA_PATH, 'calwebb_spec3_default.cfg'),
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


@runslow
@require_bigdata
def test_run_extract_1d_only(mk_tmp_dirs):
    """Test only the extraction step. Should produce nothing
    because extraction requires resampling
    """
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    asn_path = update_asn_basedir(
        path.join(DATAPATH, 'mos_udf_g235M_spec3_asn.json'),
        root=path.join(DATAPATH, 'level2b')
    )
    args = [
        path.join(SCRIPT_DATA_PATH, 'calwebb_spec3_default.cfg'),
        asn_path,
        '--steps.skymatch.skip=true',
        '--steps.outlier_detection.skip=true',
        '--steps.resample_spec.skip=true',
        '--steps.cube_build.skip=true',
    ]

    Step.from_cmdline(args)


@runslow
@require_bigdata
def test_run_nosteps(mk_tmp_dirs):
    """Test where no steps execute"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    asn_path = update_asn_basedir(
        path.join(DATAPATH, 'two_member_spec3_asn.json'),
        root=path.join(DATAPATH, 'level2b')
    )
    args = [
        path.join(SCRIPT_DATA_PATH, 'calwebb_spec3_default.cfg'),
        asn_path,
        '--steps.skymatch.skip=true',
        '--steps.outlier_detection.skip=true',
        '--steps.resample_spec.skip=true',
        '--steps.cube_build.skip=true',
        '--steps.extract_1d.skip=true',
    ]

    Step.from_cmdline(args)

    # Check for the Source-based cal name.
    with open(asn_path) as fp:
        asn = load_asn(fp)
    product_name = asn['products'][0]['name']
    assert path.isfile(product_name + '_cal.fits')


@pytest.mark.skip(
    reason='Dies with crds error no META.INSTRUMENT.NAME'
)
@runslow
@require_bigdata
def test_run_full(mk_tmp_dirs):
    """Test a full run"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    asn_path = update_asn_basedir(
        path.join(DATAPATH, 'mos_udf_g235M_spec3_asn.json'),
        root=path.join(DATAPATH, 'level2b')
    )
    args = [
        path.join(SCRIPT_DATA_PATH, 'calwebb_spec3_default.cfg'),
        asn_path,
    ]

    Step.from_cmdline(args)
    assert False
