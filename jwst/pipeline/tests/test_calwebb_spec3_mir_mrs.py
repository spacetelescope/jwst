"""Test calwebb_spec3 with MIRI MRS

Notes
-----
The test data has been arranged to match that of the pandokia jwst_test_data
file structure. As such the environmental variable TEST_BIGDATA points to
the top of the example data tree.
"""

from collections import defaultdict
from glob import glob
from os import path

from .helpers import (
    SCRIPT_DATA_PATH,
    abspath,
    require_bigdata,
    runslow,
    update_asn_basedir,
)

from ...associations import load_asn
from ...stpipe.step import (Step, remove_suffix)

DATAPATH = abspath(
    '$TEST_BIGDATA/miri/test_datasets/mrs/simulated'
)


@require_bigdata
def test_run_nothing(mk_tmp_dirs):
    """Run no steps. There should be no output."""

    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    asn_path = update_asn_basedir(
        path.join(DATAPATH, 'single_spec3_asn.json'),
        root=path.join(DATAPATH, 'level2b')
    )
    args = [
        path.join(SCRIPT_DATA_PATH, 'calwebb_spec3_default.cfg'),
        asn_path,
        '--steps.mrs_imatch.skip=true',
        '--steps.outlier_detection.skip=true',
        '--steps.resample_spec.skip=true',
        '--steps.cube_build.skip=true',
        '--steps.extract_1d.skip=true'
    ]

    Step.from_cmdline(args)

    assert len(glob('*')) == 0


@runslow
@require_bigdata
def test_run_extract_1d_only(mk_tmp_dirs):
    """Test only the extraction step.
    """
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    asn_path = update_asn_basedir(
        path.join(DATAPATH, 'single_spec3_asn.json'),
        root=path.join(DATAPATH, 'level2b')
    )
    args = [
        path.join(SCRIPT_DATA_PATH, 'calwebb_spec3_default.cfg'),
        asn_path,
        '--steps.mrs_imatch.skip=true',
        '--steps.outlier_detection.skip=true',
    ]

    Step.from_cmdline(args)

    with open(asn_path) as fd:
        asn = load_asn(fd)
    product_name_base = asn['products'][0]['name']
    product_name_glob = product_name_base + '_ch[34]-long_s3d.fits'
    assert len(glob(product_name_glob)) == 2
    product_name_glob = product_name_base + '_ch[34]-long_x1d.fits'
    assert len(glob(product_name_glob)) == 2

@runslow
@require_bigdata
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
        '--steps.mrs_imatch.skip=true',
        '--steps.outlier_detection.skip=true',
        '--steps.extract_1d.skip=true',
    ]

    Step.from_cmdline(args)

    with open(asn_path) as fd:
        asn = load_asn(fd)
    product_name_base = asn['products'][0]['name']
    product_name_glob = product_name_base + '_ch[34]-long_s3d.fits'
    assert len(glob(product_name_glob)) == 2


@runslow
@require_bigdata
def test_run_mrs_imatch_only(mk_tmp_dirs):
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
        '--steps.mrs_imatch.save_results=true',
    ]

    Step.from_cmdline(args)

    with open(asn_path) as fd:
        asn = load_asn(fd)
    product_name_base = asn['products'][0]['name']
    product_name = product_name_base + '_mrs_imatch.fits'
    assert path.isfile(product_name)


@runslow
@require_bigdata
def test_run_full(mk_tmp_dirs):
    """Test a basic run"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    asn_path = update_asn_basedir(
        path.join(DATAPATH, 'single_channel_spec3_asn.json'),
        root=path.join(DATAPATH, 'level2b')
    )
    args = [
        path.join(SCRIPT_DATA_PATH, 'calwebb_spec3_default.cfg'),
        asn_path,
    ]

    Step.from_cmdline(args)

    # Now test for file existence. Get the association
    with open(asn_path) as fh:
        asn = load_asn(fh)
    acid = asn['asn_id']
    product = asn['products'][0]
    product_name = product['name']
    members_by_type = defaultdict(list)
    for member in product['members']:
        expname = path.split(member['expname'])[1]
        members_by_type[member['exptype'].lower()].append(expname)

    output_files = [
        path.split(result_path)[1]
        for result_path in
        glob('*')
    ]
    print('Created files ares: {}'.format(output_files))

    # Check Level3 products
    level3_suffixes = [
        '_ch1-long_s3d.fits',
        '_ch2-long_s3d.fits',
        '_ch3-long_s3d.fits',
        '_ch4-long_s3d.fits',
        '_ch1-long_x1d.fits',
        '_ch2-long_x1d.fits',
        '_ch3-long_x1d.fits',
        '_ch4-long_x1d.fits',
    ]
    for level3_suffix in level3_suffixes:
        product_name_file = product_name + level3_suffix
        assert product_name_file in output_files
        output_files.remove(product_name_file)

    # Check Level2 products
    for member in members_by_type['science']:
        basename, ext = path.splitext(path.split(member)[1])
        basename, separator = remove_suffix(basename)

        name = basename + separator + acid + separator + 'crf' + ext
        assert name in output_files
        output_files.remove(name)

    # If there are files left, this is an error
    assert len(output_files) == 0


@runslow
@require_bigdata
def test_run_outlier_only(mk_tmp_dirs):
    """Test a basic run"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    asn_path = update_asn_basedir(
        path.join(DATAPATH, 'single_channel_spec3_asn.json'),
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
    # Now test for file existence. Get the association
    with open(asn_path) as fh:
        asn = load_asn(fh)
    acid = asn['asn_id']
    product = asn['products'][0]
    members_by_type = defaultdict(list)
    for member in product['members']:
        expname = path.split(member['expname'])[1]
        members_by_type[member['exptype'].lower()].append(expname)

    output_files = [
        path.split(result_path)[1]
        for result_path in
        glob('*')
    ]
    print('Created files ares: {}'.format(output_files))

    # Check Level2 products
    for member in members_by_type['science']:
        basename, ext = path.splitext(path.split(member)[1])
        basename, separator = remove_suffix(basename)

        name = basename + separator + acid + separator + 'crf' + ext
        assert name in output_files
        output_files.remove(name)

    # If there are files left, this is an error
    assert len(output_files) == 0
