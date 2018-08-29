"""Test calwebb_detector1 with NIRCam"""

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

    input_file = 'jw82500001003_02101_00001_NRCALONG_uncal.fits'
    args = [
        path.join(SCRIPT_DATA_PATH, 'cfgs', 'calwebb_detector1.cfg'),
        path.join(DATAPATH, input_file)
    ]

    Step.from_cmdline(args)

    # Check for output
    output_files = [
        path.split(result_path)[1]
        for result_path in
        glob('*')
    ]
    print('Created files ares: {}'.format(output_files))

    # Check for primary pipeline result
    input_basename, ext = path.splitext(input_file)
    output_basename, separator = remove_suffix(input_basename)
    output_file = output_basename + separator + 'rate' + ext
    assert output_file in output_files
    output_files.remove(output_file)

    output_file = output_basename + separator + 'rateints' + ext
    assert output_file in output_files
    output_files.remove(output_file)

    # Check for trapsfilled file.
    output_file = output_basename + separator + 'trapsfilled' + ext
    assert output_file in output_files
    output_files.remove(output_file)

    # If there are files left, this is an error
    assert len(output_files) == 0


@runslow
@require_bigdata
def test_run_persistence_with_trapsfilled(mk_tmp_dirs):
    """Test a full run"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    input_file = 'jw82500001003_02101_00001_NRCALONG_uncal.fits'
    args = [
        path.join(SCRIPT_DATA_PATH, 'cfgs', 'calwebb_detector1.cfg'),
        path.join(DATAPATH, input_file),
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
        glob('*')
    ]
    print('Created files ares: {}'.format(output_files))

    # Check for primary pipeline result
    input_basename, ext = path.splitext(input_file)
    output_basename, separator = remove_suffix(input_basename)
    output_file = output_basename + separator + 'ramp' + ext
    assert output_file in output_files
    output_files.remove(output_file)

    # Check for trapsfilled file.
    output_file = output_basename + separator + 'trapsfilled' + ext
    assert output_file in output_files
    output_files.remove(output_file)

    # If there are files left, this is an error
    assert len(output_files) == 0


@runslow
@require_bigdata
def test_run_persistence_without_trapsfilled(mk_tmp_dirs):
    """Test a full run"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    input_file = 'jw82500001003_02101_00001_NRCALONG_uncal.fits'
    args = [
        path.join(SCRIPT_DATA_PATH, 'cfgs', 'calwebb_detector1.cfg'),
        path.join(DATAPATH, input_file),
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
        '--steps.persistence.save_trapsfilled=false',
        '--steps.jump.skip=true',
        '--steps.ramp_fit.skip=true',
        '--steps.gain_scale.skip=true',
    ]

    Step.from_cmdline(args)

    # Check for output
    output_files = [
        path.split(result_path)[1]
        for result_path in
        glob('*')
    ]
    print('Created files ares: {}'.format(output_files))

    # Check for primary pipeline result
    input_basename, ext = path.splitext(input_file)
    output_basename, separator = remove_suffix(input_basename)
    output_file = output_basename + separator + 'ramp' + ext
    assert output_file in output_files
    output_files.remove(output_file)

    # If there are files left, this is an error
    assert len(output_files) == 0
