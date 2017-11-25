"""Test step/pipeline saving"""
import os
from os.path import (
    dirname,
    isfile,
    join,
    splitext,
)
import shutil
import tempfile

import pytest

from ..step import Step

data_fn = 'flat.fits'
data_fn_path = join(dirname(__file__), 'data', data_fn)
data_name, data_ext = splitext(data_fn)


@pytest.fixture
def mk_tmp_dirs():
    tmp_current_path = tempfile.mkdtemp()
    tmp_data_path = tempfile.mkdtemp()
    tmp_config_path = tempfile.mkdtemp()

    old_path = os.getcwd()
    try:
        os.chdir(tmp_current_path)
        yield (tmp_current_path, tmp_data_path, tmp_config_path)
    finally:
        os.chdir(old_path)


def test_save_step_default(mk_tmp_dirs):
    """Default save should be current working directory"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    args = [
        'jwst.stpipe.tests.steps.StepWithModel',
        data_fn_path
    ]

    Step.from_cmdline(args)

    fname = 'flat_stepwithmodel.fits'
    assert isfile(fname)


def test_save_step_withoutput(mk_tmp_dirs):
    """Default save should be current working directory"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    output_file = 'junk.fits'

    args = [
        'jwst.stpipe.tests.steps.StepWithModel',
        data_fn_path,
        '--output_file=' + output_file
    ]

    Step.from_cmdline(args)

    assert isfile(output_file)


def test_save_step_withdir(mk_tmp_dirs):
    """Save to specified folder"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    args = [
        'jwst.stpipe.tests.steps.StepWithModel',
        data_fn_path,
        '--output_dir=' + tmp_data_path
    ]

    Step.from_cmdline(args)

    output_fn_path = join(
        tmp_data_path,
        data_name + '_stepwithmodel' + data_ext,
    )
    assert isfile(output_fn_path)


def test_save_step_withdir_withoutput(mk_tmp_dirs):
    """Save to specified folder"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    output_file = 'junk.fits'

    args = [
        'jwst.stpipe.tests.steps.StepWithModel',
        data_fn_path,
        '--output_dir=' + tmp_data_path,
        '--output_file=' + output_file
    ]

    Step.from_cmdline(args)

    output_fn_path = join(
        tmp_data_path,
        output_file
    )
    assert isfile(output_fn_path)


def test_save_pipeline_default(mk_tmp_dirs):
    """Default save should be current working directory"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs
    step_fn = 'save_pipeline.cfg'

    step_fn_path = join(dirname(__file__), 'steps', step_fn)

    tmp_step_fn_path = join(tmp_config_path, step_fn)
    tmp_data_fn_path = join(tmp_data_path, data_fn)
    shutil.copy(step_fn_path, tmp_step_fn_path)
    shutil.copy(data_fn_path, tmp_data_fn_path)

    args = [
        step_fn_path,
        tmp_data_fn_path,
        '--steps.savestep.skip=False'
    ]

    Step.from_cmdline(args)

    output_pipeline_fn_path = data_name + '_processed_savepipeline' + data_ext
    output_stepsave_fn_path = data_name + '_processed' + data_ext
    assert isfile(output_pipeline_fn_path)
    assert isfile(output_stepsave_fn_path)


def test_save_pipeline_withdir(mk_tmp_dirs):
    """Save to specified folder"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs
    step_fn = 'save_pipeline.cfg'
    step_fn_path = join(dirname(__file__), 'steps', step_fn)

    tmp_step_fn_path = join(tmp_config_path, step_fn)
    shutil.copy(step_fn_path, tmp_step_fn_path)

    args = [
        step_fn_path,
        data_fn_path,
        '--output_dir=' + tmp_data_path,
    ]

    Step.from_cmdline(args)

    output_pipeline_fn_path = join(
        tmp_data_path,
        data_name + '_savepipeline' + data_ext
    )
    assert isfile(output_pipeline_fn_path)


def test_save_substep_withdir(mk_tmp_dirs):
    """Save to specified folder"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs
    step_fn = 'save_pipeline.cfg'
    step_fn_path = join(dirname(__file__), 'steps', step_fn)

    tmp_step_fn_path = join(tmp_config_path, step_fn)
    shutil.copy(step_fn_path, tmp_step_fn_path)

    args = [
        step_fn_path,
        data_fn_path,
        '--steps.savestep.skip=False',
        '--steps.savestep.output_dir=' + tmp_data_path
    ]

    Step.from_cmdline(args)

    output_pipeline_fn_path = data_name + '_processed_savepipeline' + data_ext
    output_stepsave_fn_path = join(
        tmp_data_path,
        data_name + '_processed' + data_ext
    )
    assert isfile(output_pipeline_fn_path)
    assert isfile(output_stepsave_fn_path)


def test_save_proper_pipeline(mk_tmp_dirs):
    """Test how pipeline saving should work"""
    args = [
        'jwst.stpipe.tests.steps.ProperPipeline',
        data_fn_path,
    ]
    Step.from_cmdline(args)

    assert isfile('ppbase_pp.fits')


def test_save_proper_pipeline_withdir(mk_tmp_dirs):
    """Test how pipeline saving should work with output_dir"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    args = [
        'jwst.stpipe.tests.steps.ProperPipeline',
        data_fn_path,
        '--output_dir=' + tmp_data_path
    ]
    Step.from_cmdline(args)

    assert isfile(join(tmp_data_path, 'ppbase_pp.fits'))


def test_save_proper_pipeline_withdir_withoutput(mk_tmp_dirs):
    """Test how pipeline saving should work with output_dir"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    output_name = 'junk.fits'

    args = [
        'jwst.stpipe.tests.steps.ProperPipeline',
        data_fn_path,
        '--output_file=' + output_name,
        '--output_dir=' + tmp_data_path
    ]
    Step.from_cmdline(args)

    assert isfile(join(tmp_data_path, output_name))


def test_save_proper_pipeline_substeps(mk_tmp_dirs):
    """Test how pipeline saving should work"""
    args = [
        'jwst.stpipe.tests.steps.ProperPipeline',
        data_fn_path,
        '--steps.stepwithmodel.save_results=true',
        '--steps.another_stepwithmodel.save_results=true',
    ]
    Step.from_cmdline(args)

    assert isfile('ppbase_pp.fits')
    assert isfile('ppbase_swm.fits')
    assert isfile('ppbase_aswm.fits')


def test_save_proper_pipeline_substeps_skip(mk_tmp_dirs):
    """Test how pipeline saving should work"""
    args = [
        'jwst.stpipe.tests.steps.ProperPipeline',
        data_fn_path,
        '--steps.stepwithmodel.save_results=true',
        '--steps.another_stepwithmodel.save_results=true',
        '--steps.another_stepwithmodel.skip=true',
    ]
    Step.from_cmdline(args)

    assert isfile('ppbase_pp.fits')
    assert isfile('ppbase_swm.fits')
    assert not isfile('ppbase_aswm.fits')


def test_save_proper_pipeline_substeps_withdir(mk_tmp_dirs):
    """Test how pipeline saving should work"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    args = [
        'jwst.stpipe.tests.steps.ProperPipeline',
        data_fn_path,
        '--output_dir=' + tmp_data_path,
        '--steps.stepwithmodel.save_results=true',
        '--steps.another_stepwithmodel.save_results=true',
        '--steps.another_stepwithmodel.output_dir=' + tmp_config_path,
    ]
    Step.from_cmdline(args)

    assert isfile(join(tmp_data_path, 'ppbase_pp.fits'))
    assert isfile(join(tmp_data_path, 'ppbase_swm.fits'))
    assert isfile(join(tmp_config_path, 'ppbase_aswm.fits'))
