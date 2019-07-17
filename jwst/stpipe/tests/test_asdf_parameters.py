"""Test initializing steps using ASDF and CRDS"""
from pathlib import Path
import pytest

from jwst.stpipe.config_parser import ValidationError

from .steps import MakeListStep
from .util import t_path

DEFAULT_PAR1 = 42.0
DEFAULT_PAR2 = 'Yes, a string'


def test_asdf_from_call():
    """Test using an ASDF file from call"""
    config_file = t_path(Path('data') / 'step_parameters' / 'jwst_all_pars-makeliststep_0001.asdf')
    results = MakeListStep.call(config_file=config_file)

    assert results == [42.0, 'Yes, a string', False]


def test_makeliststep_missingpars():
    """Test the testing step class when given insufficient information"""
    with pytest.raises(ValidationError):
        MakeListStep.call()


def test_makeliststep_test():
    """Test the testing step class for basic operation"""
    result = MakeListStep.call(par1=DEFAULT_PAR1, par2=DEFAULT_PAR2)

    assert result == [42.0, 'Yes, a string', False]
