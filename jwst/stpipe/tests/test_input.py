"""Test input directory usage and input defaults"""

from os import path
import pytest

from jwst.stpipe import Step
from jwst.datamodels import JwstDataModel, ModelContainer
from jwst import datamodels

from .steps import StepWithModel
from .util import t_path


def test_default_input_with_container(mk_tmp_dirs):
    """Test default input name from a ModelContainer"""

    model_path = t_path('data/flat.fits')
    with ModelContainer([model_path]) as container:
        step = StepWithModel()
        step.run(container)

        assert step._input_filename is None


def test_default_input_with_full_model():
    """Test default input name retrieval with actual model"""
    model_path = t_path('data/flat.fits')
    with datamodels.open(model_path) as model:
        step = StepWithModel()
        step.run(model)

        assert step._input_filename == model.meta.filename


def test_default_input_with_new_model():
    """Test getting input name with new model"""

    step = StepWithModel()

    model = JwstDataModel()
    step.run(model)

    assert step._input_filename is None


def test_default_input_dir(mk_tmp_dirs):
    """Test defaults"""
    input_file = t_path('data/flat.fits')

    step = Step.from_cmdline([
        'jwst.stpipe.tests.steps.StepWithModel',
        input_file
    ])

    # Check that `input_dir` is set.
    input_path = path.split(input_file)[0]
    assert step.input_dir == input_path


def test_set_input_dir(mk_tmp_dirs):
    """Simply set the path"""
    input_file = t_path('data/flat.fits')

    step = Step.from_cmdline([
        'jwst.stpipe.tests.steps.StepWithModel',
        input_file,
        '--input_dir', 'junkdir'
    ])

    # Check that `input_dir` is set.
    assert step.input_dir == 'junkdir'


def test_use_input_dir(mk_tmp_dirs):
    """Test with a specified path"""
    input_dir = t_path('data')
    input_file = 'flat.fits'

    step = Step.from_cmdline([
        'jwst.stpipe.tests.steps.StepWithModel',
        input_file,
        '--input_dir', input_dir
    ])

    # Check that `input_dir` is set.
    assert step.input_dir == input_dir


def test_fail_input_dir(mk_tmp_dirs):
    """Fail with a bad file path"""
    input_file = 'flat.fits'

    with pytest.raises(FileNotFoundError):
        Step.from_cmdline([
            'jwst.stpipe.tests.steps.StepWithModel',
            input_file,
        ])


def test_input_dir_with_model(mk_tmp_dirs):
    """Use with an already opened DataModel"""
    with datamodels.open(t_path('data/flat.fits')) as model:
        step = StepWithModel()
        step.run(model)

        assert step.input_dir == ''
