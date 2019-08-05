"""Test initializing steps using ASDF and CRDS"""
from pathlib import Path
import pytest

from jwst.stpipe.config_parser import ValidationError
from jwst.stpipe.step import Step

from .steps import MakeListStep
from .util import t_path

DEFAULT_PAR1 = 42.0
DEFAULT_PAR2 = 'Yes, a string'
DEFAULT_RESULT = [DEFAULT_PAR1, DEFAULT_PAR2, False]


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
