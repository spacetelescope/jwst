"""Test input directory usage"""
from os import path

from .util import (mk_tmp_dirs, t_path)
from ..step import Step


def test_default_input_dir(mk_tmp_dirs):

    input_file = t_path('data/flat.fits')

    step = Step.from_cmdline([
        'jwst.stpipe.tests.steps.StepWithModel',
        input_file
    ])

    # Check that `input_dir` is set.
    input_path = path.split(input_file)[0]
    assert step.input_dir == input_path


def test_set_input_dir(mk_tmp_dirs):

    input_file = t_path('data/flat.fits')

    step = Step.from_cmdline([
        'jwst.stpipe.tests.steps.StepWithModel',
        input_file,
        '--input_dir', 'junkdir'
    ])

    # Check that `input_dir` is set.
    assert step.input_dir == 'junkdir'
