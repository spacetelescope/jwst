"""Test step/pipeline saving"""
from glob import glob
import os
from os import path
import shutil

from ..step import Step

data_fn = 'flat.fits'
data_fn_path = path.join(path.dirname(__file__), 'data', data_fn)
data_name, data_ext = path.splitext(data_fn)


def test_make_output_path():
    """Test the basic make_output_file method"""

    step = Step()
    output_path = step.make_output_path('junk_uncal.fits')
    assert output_path == 'junk_step.fits'

    output_path = step.make_output_path('junk_uncal.fits', idx=1)
    assert output_path == 'junk_1_step.fits'

    step.output_ext = '.asdf'
    output_path = step.make_output_path('junk_uncal')
    assert output_path == 'junk_step.asdf'

    output_path = step.make_output_path('junk_uncal.fits', ext='asdf')
    assert output_path == 'junk_step.asdf'

    step.output_dir = '/junk'
    step.output_ext = None
    output_path = step.make_output_path('junk_uncal.fits')
    assert output_path == path.join(step.output_dir, 'junk_step.fits')


def test_save_step_default(mk_tmp_dirs):
    """Default save should be current working directory"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    args = [
        'jwst.stpipe.tests.steps.StepWithModel',
        data_fn_path
    ]

    step = Step.from_cmdline(args)
    step.closeout()

    fname = 'flat_stepwithmodel.fits'
    assert path.isfile(fname)


def test_save_step_withoutput(mk_tmp_dirs):
    """Default save should be current working directory"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    output_file = 'junk.fits'

    args = [
        'jwst.stpipe.tests.steps.StepWithModel',
        data_fn_path,
        '--output_file=' + output_file
    ]

    step = Step.from_cmdline(args)
    step.closeout()

    output_path, output_ext = path.splitext(output_file)
    assert path.isfile(output_path + '_stepwithmodel' + output_ext)


def test_save_step_withoutputsuffix(mk_tmp_dirs):
    """Default save should be current working directory"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    output_file = 'junk_rate.fits'
    actual_output_file = 'junk_stepwithmodel.fits'

    args = [
        'jwst.stpipe.tests.steps.StepWithModel',
        data_fn_path,
        '--output_file=' + output_file
    ]

    step = Step.from_cmdline(args)
    step.closeout()

    assert path.isfile(actual_output_file)


def test_save_step_withdir(mk_tmp_dirs):
    """Save to specified folder"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    args = [
        'jwst.stpipe.tests.steps.StepWithModel',
        data_fn_path,
        '--output_dir=' + tmp_data_path
    ]

    step = Step.from_cmdline(args)
    step.closeout()

    output_fn_path = path.join(
        tmp_data_path,
        data_name + '_stepwithmodel' + data_ext,
    )
    assert path.isfile(output_fn_path)


def test_save_step_withdir_environment(mk_tmp_dirs):
    """Save to specified folder"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    os.environ['TSSWE_OUTPATH'] = tmp_data_path

    args = [
        'jwst.stpipe.tests.steps.StepWithModel',
        data_fn_path,
        '--output_dir=$TSSWE_OUTPATH'
    ]

    step = Step.from_cmdline(args)
    step.closeout()

    output_fn_path = path.join(
        tmp_data_path,
        data_name + '_stepwithmodel' + data_ext,
    )
    assert path.isfile(output_fn_path)


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

    step = Step.from_cmdline(args)
    step.closeout()

    output_path, output_ext = path.splitext(output_file)
    output_fn_path = path.join(
        tmp_data_path,
        output_path + '_stepwithmodel' + output_ext
    )
    assert path.isfile(output_fn_path)


def test_save_container(mk_tmp_dirs):
    """Step with output_use_model"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    args = [
        'jwst.stpipe.tests.steps.StepWithContainer',
        data_fn_path,
    ]

    Step.from_cmdline(args)

    assert path.isfile('flat_0_stepwithcontainer.fits')
    assert path.isfile('flat_1_stepwithcontainer.fits')


def test_save_container_usemodel(mk_tmp_dirs):
    """Step with output_use_model"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    args = [
        'jwst.stpipe.tests.steps.StepWithContainer',
        data_fn_path,
        '--output_use_model=true'
    ]

    Step.from_cmdline(args)

    assert path.isfile('swc_model1_stepwithcontainer.fits')
    assert path.isfile('swc_model2_stepwithcontainer.fits')


def test_save_container_withfile(mk_tmp_dirs):
    """Step with output_use_model"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    args = [
        'jwst.stpipe.tests.steps.StepWithContainer',
        data_fn_path,
        '--output_file=tscwf.fits',
    ]

    step = Step.from_cmdline(args)

    assert path.isfile('tscwf_0_stepwithcontainer.fits')
    assert path.isfile('tscwf_1_stepwithcontainer.fits')

    step.closeout()


def test_save_pipeline_default(mk_tmp_dirs):
    """Default save should be current working directory"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs
    step_fn = 'save_pipeline.cfg'

    step_fn_path = path.join(path.dirname(__file__), 'steps', step_fn)

    tmp_step_fn_path = path.join(tmp_config_path, step_fn)
    tmp_data_fn_path = path.join(tmp_data_path, data_fn)
    shutil.copy(step_fn_path, tmp_step_fn_path)
    shutil.copy(data_fn_path, tmp_data_fn_path)

    args = [
        step_fn_path,
        tmp_data_fn_path,
        '--steps.savestep.skip=False'
    ]

    Step.from_cmdline(args)

    # Output from the explicit `SaveStep.save_model`
    desired = data_name + '_processed' + data_ext
    assert path.isfile(desired)

    # Output from the Step's default saving of `SaveStep`
    desired = data_name + '_savestep' + data_ext
    assert path.isfile(desired)

    # Output from the Steps' default saving of `SavePipeline`
    desired = data_name + '_savepipeline' + data_ext
    assert path.isfile(desired)


def test_save_pipeline_withdir(mk_tmp_dirs):
    """Save to specified folder"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs
    step_fn = 'save_pipeline.cfg'
    step_fn_path = path.join(path.dirname(__file__), 'steps', step_fn)

    tmp_step_fn_path = path.join(tmp_config_path, step_fn)
    shutil.copy(step_fn_path, tmp_step_fn_path)

    args = [
        step_fn_path,
        data_fn_path,
        '--output_dir=' + tmp_data_path,
    ]

    Step.from_cmdline(args)

    output_pipeline_fn_path = path.join(
        tmp_data_path,
        data_name + '_savepipeline' + data_ext
    )
    assert path.isfile(output_pipeline_fn_path)


def test_save_substep_withdir(mk_tmp_dirs):
    """Save to specified folder"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs
    step_fn = 'save_pipeline.cfg'
    step_fn_path = path.join(path.dirname(__file__), 'steps', step_fn)

    tmp_step_fn_path = path.join(tmp_config_path, step_fn)
    shutil.copy(step_fn_path, tmp_step_fn_path)

    args = [
        step_fn_path,
        data_fn_path,
        '--steps.savestep.skip=False',
        '--steps.savestep.output_dir=' + tmp_data_path
    ]

    Step.from_cmdline(args)

    # Output from the explicit `SaveStep.save_model`
    desired = path.join(
        tmp_data_path,
        data_name + '_processed' + data_ext
    )
    assert path.isfile(desired)

    # Output from the Step's default saving of `SaveStep`
    desired = path.join(
        tmp_data_path,
        data_name + '_savestep' + data_ext
    )
    assert path.isfile(desired)

    # Output from the Steps' default saving of `SavePipeline`
    desired = data_name + '_savepipeline' + data_ext
    assert path.isfile(desired)


def test_save_proper_pipeline(mk_tmp_dirs):
    """Test how pipeline saving should work"""
    args = [
        'jwst.stpipe.tests.steps.ProperPipeline',
        data_fn_path,
        '--steps.stepwithcontainer.skip=true',
    ]

    Step.from_cmdline(args)

    assert path.isfile('flat_pp.fits')


def test_save_proper_pipeline_withdir(mk_tmp_dirs):
    """Test how pipeline saving should work with output_dir"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    args = [
        'jwst.stpipe.tests.steps.ProperPipeline',
        data_fn_path,
        '--output_dir=' + tmp_data_path,
        '--steps.stepwithcontainer.skip=true',
    ]

    Step.from_cmdline(args)

    assert path.isfile(path.join(tmp_data_path, 'flat_pp.fits'))


def test_save_proper_pipeline_withdir_withoutput(mk_tmp_dirs):
    """Test how pipeline saving should work with output_dir"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    output_name = 'junk.fits'

    args = [
        'jwst.stpipe.tests.steps.ProperPipeline',
        data_fn_path,
        '--output_file=' + output_name,
        '--output_dir=' + tmp_data_path,
        '--steps.stepwithcontainer.skip=true',
    ]

    Step.from_cmdline(args)

    output_path, output_ext = path.splitext(output_name)
    assert path.isfile(path.join(
        tmp_data_path, output_path + '_pp' + output_ext
    ))


def test_save_proper_pipeline_substeps(mk_tmp_dirs):
    """Test how pipeline saving should work"""
    args = [
        'jwst.stpipe.tests.steps.ProperPipeline',
        data_fn_path,
        '--steps.stepwithmodel.save_results=true',
        '--steps.another_stepwithmodel.save_results=true',
        '--steps.stepwithcontainer.skip=true',
    ]

    Step.from_cmdline(args)

    assert path.isfile('flat_pp.fits')
    assert path.isfile('flat_swm.fits')
    assert path.isfile('flat_aswm.fits')


def test_save_proper_pipeline_substeps_skip(mk_tmp_dirs):
    """Test how pipeline saving should work"""
    args = [
        'jwst.stpipe.tests.steps.ProperPipeline',
        data_fn_path,
        '--steps.stepwithmodel.save_results=true',
        '--steps.another_stepwithmodel.save_results=true',
        '--steps.another_stepwithmodel.skip=true',
        '--steps.stepwithcontainer.skip=true',
    ]

    Step.from_cmdline(args)

    assert path.isfile('flat_pp.fits')
    assert path.isfile('flat_swm.fits')
    assert not path.isfile('flat_aswm.fits')


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
        '--steps.stepwithcontainer.skip=true',
    ]

    Step.from_cmdline(args)

    assert path.isfile(path.join(tmp_data_path, 'flat_pp.fits'))
    assert path.isfile(path.join(tmp_data_path, 'flat_swm.fits'))
    assert path.isfile(path.join(tmp_config_path, 'flat_aswm.fits'))


def test_save_proper_pipeline_container(mk_tmp_dirs):
    """Test how pipeline saving should work"""
    args = [
        'jwst.stpipe.tests.steps.ProperPipeline',
        data_fn_path,
    ]

    Step.from_cmdline(args)

    assert path.isfile('flat_0_pp.fits')
    assert path.isfile('flat_1_pp.fits')


def test_save_proper_pipeline_container_withdir(mk_tmp_dirs):
    """Test how pipeline saving should work with output_dir"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    args = [
        'jwst.stpipe.tests.steps.ProperPipeline',
        data_fn_path,
        '--output_dir=' + tmp_data_path,
    ]

    Step.from_cmdline(args)

    assert path.isfile(path.join(tmp_data_path, 'flat_0_pp.fits'))
    assert path.isfile(path.join(tmp_data_path, 'flat_1_pp.fits'))


def test_save_proper_pipeline_container_withdir_withoutput(mk_tmp_dirs):
    """Test how pipeline saving should work with output_dir"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    output_name = 'junk.fits'

    args = [
        'jwst.stpipe.tests.steps.ProperPipeline',
        data_fn_path,
        '--output_file=' + output_name,
        '--output_dir=' + tmp_data_path,
    ]

    Step.from_cmdline(args)

    output_path, output_ext = path.splitext(output_name)
    assert path.isfile(path.join(
        tmp_data_path, output_path + '_0_pp' + output_ext
    ))
    assert path.isfile(path.join(
        tmp_data_path, output_path + '_1_pp' + output_ext
    ))


def test_save_proper_pipeline_container_substeps(mk_tmp_dirs):
    """Test how pipeline saving should work"""
    args = [
        'jwst.stpipe.tests.steps.ProperPipeline',
        data_fn_path,
        '--steps.stepwithmodel.save_results=true',
        '--steps.another_stepwithmodel.save_results=true',
        '--steps.stepwithcontainer.save_results=true',
    ]

    Step.from_cmdline(args)

    assert path.isfile('flat_0_pp.fits')
    assert path.isfile('flat_1_pp.fits')
    assert path.isfile('flat_swm.fits')
    assert path.isfile('flat_aswm.fits')
    assert path.isfile('flat_0_swc.fits')
    assert path.isfile('flat_1_swc.fits')


def test_save_proper_pipeline_container_substeps_skip(mk_tmp_dirs):
    """Test how pipeline saving should work"""
    args = [
        'jwst.stpipe.tests.steps.ProperPipeline',
        data_fn_path,
        '--steps.stepwithmodel.save_results=true',
        '--steps.another_stepwithmodel.save_results=true',
        '--steps.another_stepwithmodel.skip=true',
        '--steps.stepwithcontainer.save_results=true',
    ]

    Step.from_cmdline(args)

    assert path.isfile('flat_0_pp.fits')
    assert path.isfile('flat_1_pp.fits')
    assert path.isfile('flat_swm.fits')
    assert not path.isfile('flat_aswm.fits')
    assert path.isfile('flat_0_swc.fits')
    assert path.isfile('flat_1_swc.fits')


def test_save_proper_pipeline_container_substeps_withdir(mk_tmp_dirs):
    """Test how pipeline saving should work"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    args = [
        'jwst.stpipe.tests.steps.ProperPipeline',
        data_fn_path,
        '--output_dir=' + tmp_data_path,
        '--steps.stepwithmodel.save_results=true',
        '--steps.another_stepwithmodel.save_results=true',
        '--steps.another_stepwithmodel.output_dir=' + tmp_config_path,
        '--steps.stepwithcontainer.save_results=true',
    ]

    Step.from_cmdline(args)

    assert path.isfile(path.join(tmp_data_path, 'flat_0_pp.fits'))
    assert path.isfile(path.join(tmp_data_path, 'flat_1_pp.fits'))
    assert path.isfile(path.join(tmp_data_path, 'flat_swm.fits'))
    assert path.isfile(path.join(tmp_config_path, 'flat_aswm.fits'))
    assert path.isfile(path.join(tmp_data_path, 'flat_0_swc.fits'))
    assert path.isfile(path.join(tmp_data_path, 'flat_1_swc.fits'))


def test_save_proper_pipeline_container_usemodel(mk_tmp_dirs):
    """Test how pipeline saving should work"""
    args = [
        'jwst.stpipe.tests.steps.ProperPipeline',
        data_fn_path,
        '--output_use_model=true',
        '--steps.stepwithcontainer.output_use_model=true',
    ]

    Step.from_cmdline(args)

    output_files = [
        path.split(result_path)[1]
        for result_path in
        glob('*')
    ]
    print('Created files ares: {}'.format(output_files))

    valid_files = [
        'swc_model1_pp.fits',
        'swc_model2_pp.fits'
    ]
    for valid_file in valid_files:
        assert valid_file in output_files
        output_files.remove(valid_file)

    assert len(output_files) == 0


def test_save_proper_pipeline_container_nosearch(mk_tmp_dirs):
    """Test how pipeline saving should work"""
    args = [
        'jwst.stpipe.tests.steps.ProperPipeline',
        data_fn_path,
        '--steps.stepwithcontainer.save_results=true',
        '--steps.stepwithcontainer.search_output_file=false',
    ]

    Step.from_cmdline(args)

    output_files = [
        path.split(result_path)[1]
        for result_path in
        glob('*')
    ]
    print('Created files ares: {}'.format(output_files))

    valid_files = [
        'flat_0_pp.fits',
        'flat_1_pp.fits',
        'swc_model1_swc.fits',
        'swc_model2_swc.fits'
    ]
    for valid_file in valid_files:
        assert valid_file in output_files
        output_files.remove(valid_file)

    assert len(output_files) == 0
