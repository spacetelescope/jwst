"""Test initializing steps using ASDF and CRDS"""
import pytest

from jwst.stpipe.config_parser import ValidationError

from .steps import MakeListStep

DEFAULT_PAR1 = 42.0
DEFAULT_PAR2 = 'Yes, a string'


def test_makeliststep_missingpars():
    """Test the testing step class when given insufficient information"""
    with pytest.raises(ValidationError):
        MakeListStep.call()


def test_makeliststep_test():
    """Test the testing step class for basic operation"""
    result = MakeListStep.call(par1=DEFAULT_PAR1, par2=DEFAULT_PAR2)

    assert result == [42.0, 'Yes, a string', False]
