"""Test calwebb_detector1 with NIRCam"""

from collections import defaultdict
from copy import copy
from glob import glob
from os import path

from .helpers import (
    SCRIPT_PATH,
    SCRIPT_DATA_PATH,
    abspath,
    mk_tmp_dirs,
    require_bigdata,
    runslow,
    update_asn_basedir,
)

from ...stpipe.step import (Step, remove_suffix)

DATAPATH = abspath(
    path.join('$TEST_BIGDATA', 'pipelines')
)


@runslow
@require_bigdata
def test_run_full(mk_tmp_dirs):
    """Test a full run"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    args = [
        path.join(SCRIPT_DATA_PATH, 'cfgs', 'calwebb_detector1.cfg'),
        path.join(DATAPATH, 'jw82500001003_02101_00001_NRCALONG_uncal.fits')
    ]

    Step.from_cmdline(args)

    # Check for output
    output_files = [
        path.split(result_path)[1]
        for result_path in
        glob(path.join(tmp_data_path, '*'))
    ]
    print('Created files ares: {}'.format(output_files))

    # If there are files left, this is an error
    assert len(output_files) == 0


@runslow
@require_bigdata
def test_run_persistence(mk_tmp_dirs):
    """Test a full run"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    args = [
        path.join(SCRIPT_DATA_PATH, 'cfgs', 'calwebb_detector1.cfg'),
        path.join(DATAPATH, 'jw82500001003_02101_00001_NRCALONG_uncal.fits'),
        '--steps.group_scale.skip=true',
        '--steps.dq_init.skip=true',
        '--steps.saturation.skip=true',
        '--steps.ipc.skip=true',
        '--steps.superbias.skip=true',
        '--steps.refpix.skip=true',
        '--steps.rscd.skip=true',
        '--steps.firstframe.skip=true',
        '--steps.lastframe.skip=true',
        '--steps.linearity.skip=true',
        '--steps.dark_current.skip=true',
        '--steps.persistence.skip=false',
        '--steps.jump.skip=true',
        '--steps.ramp_fit.skip=true',
        '--steps.gain_scale.skip=true',
    ]

    Step.from_cmdline(args)

    # Check for output
    output_files = [
        path.split(result_path)[1]
        for result_path in
        glob(path.join(tmp_data_path, '*'))
    ]
    print('Created files ares: {}'.format(output_files))

    # If there are not files, this is an error
    assert len(output_files) != 0

    # If there are files left, this is an error
    assert len(output_files) == 0
