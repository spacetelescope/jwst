"""Test calwebb_spec2 with NIRSpec Bright Object"""

from glob import glob
from os import path

import pytest

from .helpers import (
    SCRIPT_DATA_PATH,
    abspath,
)

from ...stpipe.step import (Step, remove_suffix)

DATAPATH = abspath(
    path.join('$TEST_BIGDATA', 'pipelines')
)


@pytest.mark.slow
@pytest.mark.bigdata
def test_run_full(mk_tmp_dirs):
    """Test a full run"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    input_file = 'jw84600042001_02101_00001_nrs2_rateints.fits'

    args = [
        path.join(SCRIPT_DATA_PATH, 'cfgs', 'calwebb_spec2.cfg'),
        path.join(DATAPATH, input_file),
        '--output_dir=' + tmp_data_path
    ]

    Step.from_cmdline(args)

    # Now test for file existence.
    output_files = [
        path.split(result_path)[1]
        for result_path in
        glob(path.join(tmp_data_path, '*'))
    ]
    print('Created files ares: {}'.format(output_files))

    # Check Level2b products
    product_name, ext = path.splitext(input_file)
    product_name, separator = remove_suffix(product_name)
    product_name_file = product_name + separator + 'calints' + ext
    assert product_name_file in output_files
    output_files.remove(product_name_file)

    product_name_file = product_name + separator + 'x1dints' + ext
    assert product_name_file in output_files
    output_files.remove(product_name_file)

    # If there are files left, this is an error
    assert len(output_files) == 0
