"""Test calwebb_spec3 with MIRI MRS with 4-dither pattern

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
)

from ...associations import load_asn
from ...stpipe.step import Step

DATAPATH = abspath(
    '$TEST_BIGDATA/miri/test_datasets/mrs/cv3_ifu_dither'
)


@pytest.mark.slow
@pytest.mark.bigdata
def test_run_cube_build_only(mk_tmp_dirs):
    """Test only the extraction step.
    """
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    asn_path = path.join(
        DATAPATH, 'level2b', 'cube_build_4dither_495_asn.json'
    )

    args = [
        path.join(SCRIPT_DATA_PATH, 'calwebb_spec3_default.cfg'),
        asn_path,
        '--output_dir=' + tmp_data_path,
        '--steps.mrs_imatch.skip=true',
        '--steps.outlier_detection.skip=true',
        '--steps.resample_spec.skip=true',
        '--steps.extract_1d.skip=true',
    ]

    Step.from_cmdline(args)

    with open(asn_path) as fd:
        asn = load_asn(fd)
    product_name_base = path.join(tmp_data_path, asn['products'][0]['name'])
    product_name = product_name_base + '_ch1-short_s3d.fits'
    assert path.isfile(product_name)
    product_name = product_name_base + '_ch2-short_s3d.fits'
    assert path.isfile(product_name)


@pytest.mark.slow
@pytest.mark.bigdata
def test_run_extract_1d_only(mk_tmp_dirs):
    """Test only the extraction step.
    """
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    asn_path = path.join(
        DATAPATH, 'level2b', 'cube_build_4dither_495_asn.json'
    )

    args = [
        path.join(SCRIPT_DATA_PATH, 'calwebb_spec3_default.cfg'),
        asn_path,
        '--output_dir=' + tmp_data_path,
        '--steps.mrs_imatch.skip=true',
        '--steps.outlier_detection.skip=true',
        '--steps.resample_spec.skip=true',
    ]

    Step.from_cmdline(args)

    with open(asn_path) as fd:
        asn = load_asn(fd)
    product_name_base = path.join(tmp_data_path, asn['products'][0]['name'])
    product_name = product_name_base + '_ch1-short_s3d.fits'
    assert path.isfile(product_name)
    product_name = product_name_base + '_ch2-short_s3d.fits'
    assert path.isfile(product_name)
    product_name = product_name_base + '_ch1-short_x1d.fits'
    assert path.isfile(product_name)
    product_name = product_name_base + '_ch2-short_x1d.fits'
    assert path.isfile(product_name)
