"""Test initializing steps using ASDF and CRDS"""
from pathlib import Path
import pytest

from stpipe.config_parser import ValidationError
from stpipe import Step

from jwst.stpipe.tests.steps import MakeListStep, MakeListPipeline
from jwst.stpipe.tests.util import t_path

DEFAULT_PAR1 = 42.0
DEFAULT_PAR2 = 'Yes, a string'
DEFAULT_RESULT = [DEFAULT_PAR1, DEFAULT_PAR2, False]


def test_asdf_roundtrip_pipeline(_jail):
    """Save a Pipeline pars and re-instantiate with the save parameters"""

    # Save the parameters
    par_path = 'mkp_pars.asdf'
    args = [
        'jwst.stpipe.tests.steps.MakeListPipeline',
        'a.fits',
        'b',
        '--steps.make_list.par1', '10.',
        '--steps.make_list.par2', 'par2',
        '--save-parameters',
        par_path
    ]
    Step.from_cmdline(args)

    # Rerun with the parameter file
    # Initial condition is that `Step.from_cmdline`
    # succeeds.
    args = [
        par_path,
        'a.fits',
        'b'
    ]
    step = Step.from_cmdline(args)

    # As a secondary condition, ensure the required parameter
    # `par2` is set.
    assert step.make_list.par2 == 'par2'


def test_asdf_from_call():
    """Test using an ASDF file from call"""
    config_file = t_path(
        Path('steps') / 'jwst_generic_pars-makeliststep_0001.asdf'
    )
    results = MakeListStep.call(config_file=config_file)

    assert results == DEFAULT_RESULT


def test_from_command_line():
    """Test creating Step from command line using ASDF"""
    config_file = t_path(
        Path('steps') / 'jwst_generic_pars-makeliststep_0001.asdf'
    )
    args = [config_file]
    step = Step.from_cmdline(args)
    assert isinstance(step, MakeListStep)
    assert step.par1 == 42.0
    assert step.par2 == 'Yes, a string'

    results = step.run()
    assert results == DEFAULT_RESULT


def test_from_command_line_override():
    """Test creating Step from command line using ASDF"""
    config_file = t_path(
        Path('steps') / 'jwst_generic_pars-makeliststep_0001.asdf'
    )
    args = [config_file, '--par1=0.']
    step = Step.from_cmdline(args)
    assert isinstance(step, MakeListStep)
    assert step.par1 == 0.
    assert step.par2 == 'Yes, a string'

    results = step.run()
    assert results == [0., DEFAULT_PAR2, False]


def test_makeliststep_missingpars():
    """Test the testing step class when given insufficient information"""
    with pytest.raises(ValidationError):
        MakeListStep.call()


def test_makeliststep_test():
    """Test the testing step class for basic operation"""
    result = MakeListStep.call(par1=DEFAULT_PAR1, par2=DEFAULT_PAR2)

    assert result == DEFAULT_RESULT


def test_step_from_asdf():
    """Test initializing step completely from config"""
    config_file = t_path(
        Path('steps') / 'jwst_generic_pars-makeliststep_0001.asdf'
    )
    step = Step.from_config_file(config_file)
    assert isinstance(step, MakeListStep)
    assert step.name == 'make_list'

    results = step.run()
    assert results == DEFAULT_RESULT


def test_step_from_asdf_api_override():
    """Test initializing step completely from config"""
    config_file = t_path(
        Path('steps') / 'jwst_generic_pars-makeliststep_0001.asdf'
    )
    results = MakeListStep.call(config_file=config_file, par1=0.)
    assert results == [0., DEFAULT_PAR2, False]


def test_makeliststep_call_config_file():
    """Test override step asdf with .cfg"""
    config_file = t_path(
        Path('steps') / 'makelist.cfg'
    )
    results = MakeListStep.call(config_file=config_file)
    assert results == [43.0, 'My hovercraft is full of eels.', False]


def test_makeliststep_call_from_within_pipeline():
    """Test override step asdf with .cfg"""
    config_file = t_path(
        Path('steps') / 'makelist_pipeline.cfg'
    )
    results = MakeListPipeline.call(config_file=config_file)
    assert results == [43.0, 'My hovercraft is full of eels.', False]


def test_step_from_asdf_noname():
    """Test initializing step completely from config without a name specified"""
    root = 'jwst_generic_pars-makeliststep_0002'
    config_file = t_path(
        Path('steps') / (root + '.asdf')
    )
    step = Step.from_config_file(config_file)
    assert isinstance(step, MakeListStep)
    assert step.name == root

    results = step.run()
    assert results == DEFAULT_RESULT
